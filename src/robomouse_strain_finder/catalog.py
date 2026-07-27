"""In-memory index over the MMRRC bulk catalog CSV.

The source CSV carries roughly 589K rows for roughly 69K stocks -- one row per
gene/allele annotation on a stock -- so rows are aggregated by
``STRAIN/STOCK_ID`` and indexed by MGI gene accession, which is the join key
back to genes found in the knowledge graph.

Two source quirks are handled here: the header for ``RESEARCH_AREAS`` has a
trailing space, and a stock's MGI allele accession and MGI gene accession
usually live on different rows.
"""

from __future__ import annotations

import csv
import datetime as dt
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import httpx

from .models import StrainAllele, StrainHit

# Some phenotype/name fields are long; the default field limit can trip on them.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

REQUIRED_COLUMNS = frozenset(
    {
        "STRAIN/STOCK_ID",
        "STRAIN/STOCK_DESIGNATION",
        "MGI_GENE_ACCESSION_ID",
        "GENE_SYMBOL",
    }
)


#: `STRAIN/STOCK_DESIGNATION` embeds presentational HTML:
#: `B6.129P2-<i>Esr2<sup>tm1Unc</sup></i>/Mmnc`.
#:
#: `ALLELE_SYMBOL` does NOT -- its angle brackets are MGI superscript
#: nomenclature (`Esr2<tm1Unc>`, `A<y>`). Stripping anything that looks like a
#: tag would destroy 95% of allele symbols, reducing them to bare gene symbols.
#: Only the designation is cleaned, and `<sup>x</sup>` is folded into MGI's own
#: `<x>` form rather than discarded.
_HTML_TAG_RE = re.compile(r"</?(?:i|b|em|strong|sup|sub|span|br|u|p|small)\s*/?>", re.I)
_SUP_RE = re.compile(r"<sup>(.*?)</sup>", re.I | re.S)
_OPEN, _CLOSE = "\x00", "\x01"


def _strip_html(value: str) -> str:
    """Render designation markup as canonical MGI text.

    Superscripts are protected with sentinels first, so an allele superscript
    that happens to spell a tag name (`<sup>b</sup>`) is not then eaten as HTML.
    """
    protected = _SUP_RE.sub(rf"{_OPEN}\1{_CLOSE}", value)
    cleaned = _HTML_TAG_RE.sub("", protected)
    return cleaned.replace(_OPEN, "<").replace(_CLOSE, ">").strip()


_PAREN_RE = re.compile(r"\(([^)]*)\)")
_TOKEN_RE = re.compile(r"[A-Za-z0-9._]+")

#: Transgene payloads that make a line a research *tool* rather than a model of
#: the gene's biology. `Tg(Mc3r-EGFP)BX153Gsat` is a GENSAT reporter: it drives
#: GFP off the Mc3r promoter and leaves Mc3r itself intact. It is not an obesity
#: model, but it is annotated with Mc3r and carries exactly one gene, so it
#: scores a perfect gene-match ratio of 1.00 and ranks alongside real knockouts.
#:
#: Deliberately restricted to `Tg(...)` transgenes. A targeted allele that
#: happens to insert lacZ (`Gabrq<tm1Hmo>`) is still a knockout.
_TOOL_PAYLOAD_RE = re.compile(
    r"^(EGFP|GFP|mCherry|tdTomato|YFP|CFP|RFP|DsRed|Venus|luciferase|"
    r"cre|creERT|creERT2|rtTA|tTA|FLPo|FLPe|Flp|Dre)$",
    re.I,
)
_TRANSGENE_RE = re.compile(r"^Tg\(([^)]*)\)", re.I)


def _is_tool_allele(allele_symbol: str) -> bool:
    """True when a transgene's payload is a marker or recombinase."""
    match = _TRANSGENE_RE.match(allele_symbol)
    if not match:
        return False
    parts = re.split(r"[-/,;]", match.group(1))
    # The first element is the promoter; the payload is what follows.
    return any(_TOOL_PAYLOAD_RE.match(part.strip()) for part in parts[1:])


def _genes_named_by_allele(allele_symbol: str, gene_by_lower: dict[str, str]) -> set[str]:
    """Infer which of a stock's genes an allele symbol refers to.

    The catalog puts the allele accession and the gene accession on *different*
    rows and never states which allele belongs to which gene, so the link is
    recovered from MGI naming: `Esr2<tm1Unc>` names `Esr2`, and a transgene
    `Tg(Myh6-Pln)11Egk` names both `Myh6` (promoter) and `Pln` (payload).

    This is a naming heuristic, not an authoritative MGI mapping -- see
    `StrainAllele.link` for how each result is labelled.
    """
    hits: set[str] = set()
    prefix = allele_symbol.split("<", 1)[0].strip().lower()
    if prefix in gene_by_lower:
        hits.add(gene_by_lower[prefix])
    for inner in _PAREN_RE.findall(allele_symbol):
        for token in _TOKEN_RE.findall(inner):
            gene = gene_by_lower.get(token.lower())
            if gene:
                hits.add(gene)
    return hits


def _split(value: str, separator: str = "|") -> list[str]:
    return [part.strip() for part in value.split(separator) if part.strip()]


@dataclass
class Stock:
    stock_id: str
    designation: str | None = None
    rrids: set[str] = field(default_factory=set)
    strain_types: set[str] = field(default_factory=set)
    states: set[str] = field(default_factory=set)
    mutation_types: set[str] = field(default_factory=set)
    gene_symbols: set[str] = field(default_factory=set)
    mgi_gene_ids: set[str] = field(default_factory=set)
    #: MGI gene accession -> symbol, so a matched accession can find its alleles.
    gene_symbol_by_mgi: dict[str, str] = field(default_factory=dict)
    alleles: dict[str, StrainAllele] = field(default_factory=dict)
    phenotypes: set[str] = field(default_factory=set)
    pubmed_ids: set[str] = field(default_factory=set)
    research_areas: set[str] = field(default_factory=set)
    sds_url: str | None = None
    accepted_date: str | None = None

    def matched_alleles(self, matched_mgi_gene_ids: list[str]) -> list[StrainAllele]:
        """Alleles belonging to the genes that matched the query."""
        wanted = {
            self.gene_symbol_by_mgi[mgi_id].lower()
            for mgi_id in matched_mgi_gene_ids
            if mgi_id in self.gene_symbol_by_mgi
        }
        if not wanted:
            return []
        return [
            allele
            for allele in self.alleles.values()
            if any(symbol.lower() in wanted for symbol in allele.gene_symbols)
        ]

    @property
    def is_tool_line(self) -> bool:
        """Every allele is a reporter/recombinase transgene, so nothing is disrupted."""
        symbols = [a.symbol for a in self.alleles.values() if a.symbol]
        return bool(symbols) and all(_is_tool_allele(symbol) for symbol in symbols)

    def to_hit(self, matched_mgi_gene_ids: list[str], matched_symbols: list[str]) -> StrainHit:
        annotated = len(self.mgi_gene_ids)
        return StrainHit(
            stock_id=self.stock_id,
            rrids=sorted(self.rrids),
            designation=self.designation,
            strain_types=sorted(self.strain_types),
            states=sorted(self.states),
            mutation_types=sorted(self.mutation_types),
            alleles=list(self.alleles.values()),
            matched_alleles=self.matched_alleles(matched_mgi_gene_ids),
            gene_symbols=sorted(self.gene_symbols),
            matched_mgi_gene_ids=matched_mgi_gene_ids,
            matched_gene_symbols=matched_symbols,
            annotated_gene_count=annotated,
            matched_fraction=(len(matched_mgi_gene_ids) / annotated) if annotated else 0.0,
            tool_line=self.is_tool_line,
            phenotypes=sorted(self.phenotypes),
            pubmed_ids=sorted(self.pubmed_ids),
            research_areas=sorted(self.research_areas),
            sds_url=self.sds_url,
            accepted_date=self.accepted_date,
        )


class MmrrcCatalog:
    """Loaded catalog with a MGI-gene-accession index."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.stocks: dict[str, Stock] = {}
        self.by_mgi_gene: dict[str, set[str]] = defaultdict(set)
        self.rows = 0
        self.modified_at: str | None = None

    @property
    def loaded(self) -> bool:
        return bool(self.stocks)

    def load(self) -> None:
        with self.path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            try:
                raw_header = next(reader)
            except StopIteration as exc:
                raise ValueError(f"Catalog {self.path} is empty") from exc

            header = [column.strip() for column in raw_header]
            missing = REQUIRED_COLUMNS - set(header)
            if missing:
                raise ValueError(f"Catalog {self.path} is missing columns: {sorted(missing)}")
            index = {name: position for position, name in enumerate(header)}

            def cell(row: list[str], name: str) -> str:
                position = index.get(name, -1)
                if position < 0 or position >= len(row):
                    return ""
                return row[position].strip()

            for row in reader:
                if not row:
                    continue
                stock_id = cell(row, "STRAIN/STOCK_ID")
                if not stock_id:
                    continue
                self.rows += 1

                stock = self.stocks.get(stock_id)
                if stock is None:
                    stock = Stock(stock_id=stock_id)
                    self.stocks[stock_id] = stock

                stock.designation = (
                    stock.designation or _strip_html(cell(row, "STRAIN/STOCK_DESIGNATION")) or None
                )
                stock.sds_url = stock.sds_url or cell(row, "SDS_URL") or None
                accepted = cell(row, "ACCEPTED_DATE")
                if accepted and accepted != "0000-00-00" and not stock.accepted_date:
                    stock.accepted_date = accepted

                for other in _split(cell(row, "OTHER_NAMES"), ","):
                    if other.startswith("RRID:"):
                        stock.rrids.add(other)

                for value, target in (
                    (cell(row, "STRAIN_TYPE"), stock.strain_types),
                    (cell(row, "MUTATION_TYPE"), stock.mutation_types),
                ):
                    if value:
                        target.add(value)
                for state in _split(cell(row, "STATE"), ","):
                    stock.states.add(state)

                gene_symbol = cell(row, "GENE_SYMBOL")
                if gene_symbol:
                    stock.gene_symbols.add(gene_symbol)
                mgi_gene = cell(row, "MGI_GENE_ACCESSION_ID")
                if mgi_gene:
                    stock.mgi_gene_ids.add(mgi_gene)
                    self.by_mgi_gene[mgi_gene].add(stock_id)
                    if gene_symbol:
                        stock.gene_symbol_by_mgi.setdefault(mgi_gene, gene_symbol)

                allele_id = cell(row, "MGI_ALLELE_ACCESSION_ID")
                # Not HTML-stripped: the angle brackets are MGI nomenclature.
                allele_symbol = cell(row, "ALLELE_SYMBOL")
                if allele_id or allele_symbol:
                    key = allele_id or allele_symbol
                    if key not in stock.alleles:
                        stock.alleles[key] = StrainAllele(
                            mgi_allele_id=allele_id or None,
                            symbol=allele_symbol or None,
                            name=_strip_html(cell(row, "ALLELE_NAME")) or None,
                        )

                for phenotype in _split(cell(row, "MPT_IDS")):
                    stock.phenotypes.add(phenotype)
                for pmid in _split(cell(row, "PUBMED_IDS")):
                    stock.pubmed_ids.add(pmid)
                for area in _split(cell(row, "RESEARCH_AREAS"), ","):
                    stock.research_areas.add(area)

        self._link_alleles_to_genes()

        stat = self.path.stat()
        self.modified_at = dt.datetime.fromtimestamp(stat.st_mtime, dt.UTC).isoformat()

    def _link_alleles_to_genes(self) -> None:
        """Attach genes to alleles once every row for a stock has been seen.

        This cannot run during the row loop: a stock's allele row and its gene
        rows are separate, and either order is possible.
        """
        for stock in self.stocks.values():
            if not stock.alleles:
                continue
            gene_by_lower = {symbol.lower(): symbol for symbol in stock.gene_symbols}
            sole_gene = next(iter(stock.gene_symbols)) if len(stock.gene_symbols) == 1 else None
            for allele in stock.alleles.values():
                named = (
                    _genes_named_by_allele(allele.symbol, gene_by_lower) if allele.symbol else set()
                )
                if named:
                    allele.gene_symbols = sorted(named)
                    allele.link = "symbol"
                elif sole_gene:
                    allele.gene_symbols = [sole_gene]
                    allele.link = "sole-gene"

    def stocks_for_genes(
        self,
        mgi_gene_ids: dict[str, str],
        *,
        exclude_mutation_types: frozenset[str] = frozenset(),
        exclude_tool_lines: bool = False,
    ) -> list[StrainHit]:
        """Find stocks annotated with any of the given MGI gene accessions.

        `mgi_gene_ids` maps an MGI accession to a display symbol. The result is
        a gene-level join: a stock is returned because it is annotated with the
        gene, not because its allele is known to model the queried phenotype.
        """
        matches: dict[str, tuple[set[str], set[str]]] = defaultdict(lambda: (set(), set()))
        for mgi_id, symbol in mgi_gene_ids.items():
            for stock_id in self.by_mgi_gene.get(mgi_id, ()):  # noqa: B038
                genes, symbols = matches[stock_id]
                genes.add(mgi_id)
                if symbol:
                    symbols.add(symbol)

        hits = []
        for stock_id, (genes, symbols) in matches.items():
            stock = self.stocks.get(stock_id)
            if stock is None:
                continue
            # `<=` alone would drop stocks with NO mutation type at all, since
            # the empty set is a subset of everything -- 6,151 gene-annotated
            # stocks silently vanished whenever any exclusion was active.
            if (
                exclude_mutation_types
                and stock.mutation_types
                and stock.mutation_types <= exclude_mutation_types
            ):
                continue
            hit = stock.to_hit(sorted(genes), sorted(symbols))
            if exclude_tool_lines and hit.tool_line:
                continue
            hits.append(hit)

        # Rank by the gene-match ratio, not raw match count. An ENU line with
        # 86 genes that happens to contain 3 query genes is a far weaker
        # candidate than a targeted mutant whose single annotated gene is a
        # query gene -- but raw count ranks the ENU line higher.
        # Reporter and cre-driver lines score a perfect gene-match ratio but do not
        # disrupt the gene, so they sort below real mutants at equal precision
        # rather than being hidden -- sometimes a reporter is what you want.
        hits.sort(
            key=lambda hit: (
                -hit.matched_fraction,
                hit.tool_line,
                -len(hit.matched_mgi_gene_ids),
                not hit.alleles,
                hit.stock_id,
            )
        )
        return hits


def download_catalog(url: str, destination: Path, timeout: float = 600.0) -> None:
    """Stream the bulk CSV to `destination`, replacing it atomically."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        dir=destination.parent, prefix=".mmrrc-", suffix=".part", delete=False
    )
    temporary = Path(handle.name)
    try:
        with handle, httpx.stream("GET", url, timeout=timeout, follow_redirects=True) as response:
            response.raise_for_status()
            for chunk in response.iter_bytes():
                handle.write(chunk)
        shutil.move(str(temporary), destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def load_catalog(path: Path, url: str, auto_download: bool) -> MmrrcCatalog:
    if not path.exists():
        if not auto_download:
            raise FileNotFoundError(
                f"MMRRC catalog not found at {path}. Download it or set "
                "RMSF_AUTO_DOWNLOAD_CATALOG=true."
            )
        download_catalog(url, path)
    catalog = MmrrcCatalog(path)
    catalog.load()
    return catalog
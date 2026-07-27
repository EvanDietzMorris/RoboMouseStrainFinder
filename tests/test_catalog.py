from __future__ import annotations

import csv
from pathlib import Path

import pytest

from robomouse_strain_finder.catalog import MmrrcCatalog

# Mirrors the real export, including the trailing space on RESEARCH_AREAS.
HEADER = [
    "STRAIN/STOCK_ID",
    "STRAIN/STOCK_DESIGNATION",
    "OTHER_NAMES",
    "STRAIN_TYPE",
    "STATE",
    "MGI_ALLELE_ACCESSION_ID",
    "ALLELE_SYMBOL",
    "ALLELE_NAME",
    "MUTATION_TYPE",
    "CHROMOSOME",
    "MGI_GENE_ACCESSION_ID",
    "GENE_SYMBOL",
    "GENE_NAME",
    "SDS_URL",
    "ACCEPTED_DATE",
    "MPT_IDS",
    "PUBMED_IDS",
    "RESEARCH_AREAS ",
]


def _row(**overrides: str) -> list[str]:
    values = dict.fromkeys(HEADER, "")
    for key, value in overrides.items():
        values[key.rstrip() if key.rstrip() in values else key] = value
    return [values[column] for column in HEADER]


@pytest.fixture
def catalog(tmp_path: Path) -> MmrrcCatalog:
    path = tmp_path / "catalog.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        # A stock split across two rows: the allele row carries no gene ID and
        # the gene row carries no allele -- the real catalog's usual shape.
        writer.writerow(
            _row(
                **{
                    "STRAIN/STOCK_ID": "MMRRC:011620-UNC",
                    "STRAIN/STOCK_DESIGNATION": "B6.129-Gabrq/Mmnc",
                    "OTHER_NAMES": "RRID:MMRRC_011620-UNC, some other name",
                    "STRAIN_TYPE": "MSR",
                    "STATE": "CA,SP",
                    "MGI_ALLELE_ACCESSION_ID": "MGI:3512216",
                    "ALLELE_SYMBOL": "Gabrq<tm1Hmo>",
                    "MUTATION_TYPE": "TM",
                    "SDS_URL": "https://example.org/11620",
                    "ACCEPTED_DATE": "2007-01-01",
                    "MPT_IDS": "decreased prepulse inhibition [MP:0009142]",
                    "PUBMED_IDS": "12345|67890",
                }
            )
        )
        writer.writerow(
            _row(
                **{
                    "STRAIN/STOCK_ID": "MMRRC:011620-UNC",
                    "STRAIN/STOCK_DESIGNATION": "B6.129-Gabrq/Mmnc",
                    "OTHER_NAMES": "RRID:MMRRC_011620-UNC",
                    "STRAIN_TYPE": "MSR",
                    "STATE": "CA,SP",
                    "MUTATION_TYPE": "TM",
                    "MGI_GENE_ACCESSION_ID": "MGI:1888498",
                    "GENE_SYMBOL": "Gabrq",
                    "ACCEPTED_DATE": "2007-01-01",
                    "RESEARCH_AREAS ": "Neurobiology, Research Tools",
                }
            )
        )
        # A second stock on a different gene.
        writer.writerow(
            _row(
                **{
                    "STRAIN/STOCK_ID": "MMRRC:000012-UNC",
                    "STRAIN/STOCK_DESIGNATION": "B6.129S4-Ccr2/Mmnc",
                    "STRAIN_TYPE": "MSR",
                    "STATE": "CA",
                    "MGI_GENE_ACCESSION_ID": "MGI:106185",
                    "GENE_SYMBOL": "Ccr2",
                    "MUTATION_TYPE": "TM",
                }
            )
        )
    loaded = MmrrcCatalog(path)
    loaded.load()
    return loaded


def test_rows_aggregate_into_stocks(catalog: MmrrcCatalog) -> None:
    assert catalog.rows == 3
    assert len(catalog.stocks) == 2


def test_allele_and_gene_rows_merge_into_one_stock(catalog: MmrrcCatalog) -> None:
    stock = catalog.stocks["MMRRC:011620-UNC"]
    assert stock.mgi_gene_ids == {"MGI:1888498"}
    assert [allele.mgi_allele_id for allele in stock.alleles.values()] == ["MGI:3512216"]
    assert stock.rrids == {"RRID:MMRRC_011620-UNC"}
    assert stock.pubmed_ids == {"12345", "67890"}


def _load(tmp_path: Path, rows: list[dict[str, str]], name: str = "t.csv") -> MmrrcCatalog:
    path = tmp_path / name
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        for values in rows:
            writer.writerow(_row(**values))
    loaded = MmrrcCatalog(path)
    loaded.load()
    return loaded


def test_designation_html_becomes_canonical_mgi_text(tmp_path: Path) -> None:
    catalog = _load(
        tmp_path,
        [
            {
                "STRAIN/STOCK_ID": "MMRRC:000211-UNC",
                "STRAIN/STOCK_DESIGNATION": "B6.129S7-<i>Chrnb2<sup>tm1Mdb</sup></i>/Mmnc",
            }
        ],
    )
    # <i> is dropped; <sup>x</sup> folds into MGI's own <x> superscript form.
    assert catalog.stocks["MMRRC:000211-UNC"].designation == "B6.129S7-Chrnb2<tm1Mdb>/Mmnc"


def test_allele_symbol_angle_brackets_are_mgi_notation_and_survive(tmp_path: Path) -> None:
    """Regression: stripping anything tag-shaped reduced 95% of allele symbols
    to bare gene symbols -- Mecp2<tm1.1Jae> became Mecp2."""
    catalog = _load(
        tmp_path,
        [
            {
                "STRAIN/STOCK_ID": "MMRRC:000011-UCD",
                "MGI_ALLELE_ACCESSION_ID": "MGI:1",
                "ALLELE_SYMBOL": "Mecp2<tm1.1Jae>",
                "ALLELE_NAME": "methyl CpG binding protein 2; targeted mutation 1.1",
            }
        ],
    )
    assert catalog.stocks["MMRRC:000011-UCD"].alleles["MGI:1"].symbol == "Mecp2<tm1.1Jae>"


def test_a_superscript_spelling_a_tag_name_is_not_eaten(tmp_path: Path) -> None:
    catalog = _load(
        tmp_path,
        [{"STRAIN/STOCK_ID": "MMRRC:1-MU", "STRAIN/STOCK_DESIGNATION": "X-<i>Lyst<sup>b</sup></i>/J"}],
    )
    assert catalog.stocks["MMRRC:1-MU"].designation == "X-Lyst<b>/J"


def test_trailing_space_research_areas_header_is_read(catalog: MmrrcCatalog) -> None:
    stock = catalog.stocks["MMRRC:011620-UNC"]
    assert stock.research_areas == {"Neurobiology", "Research Tools"}


def test_stocks_for_genes_reports_the_matching_accession(catalog: MmrrcCatalog) -> None:
    hits = catalog.stocks_for_genes({"MGI:1888498": "Gabrq"})
    assert [hit.stock_id for hit in hits] == ["MMRRC:011620-UNC"]
    assert hits[0].matched_mgi_gene_ids == ["MGI:1888498"]
    assert hits[0].matched_gene_symbols == ["Gabrq"]


def test_unknown_gene_matches_nothing(catalog: MmrrcCatalog) -> None:
    assert catalog.stocks_for_genes({"MGI:9999999": "Nope"}) == []


def test_multi_gene_matches_rank_first(catalog: MmrrcCatalog) -> None:
    hits = catalog.stocks_for_genes({"MGI:1888498": "Gabrq", "MGI:106185": "Ccr2"})
    assert {hit.stock_id for hit in hits} == {"MMRRC:011620-UNC", "MMRRC:000012-UNC"}


@pytest.fixture
def mixed_catalog(tmp_path: Path) -> MmrrcCatalog:
    """A targeted knockout of one gene, and an ENU line carrying 20 genes.

    The ENU line matches two query genes to the knockout's one, so raw match
    count would rank the noisy line first.
    """
    path = tmp_path / "mixed.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        writer.writerow(
            _row(
                **{
                    "STRAIN/STOCK_ID": "MMRRC:000001-UNC",
                    "MGI_GENE_ACCESSION_ID": "MGI:100",
                    "GENE_SYMBOL": "Targeted",
                    "MUTATION_TYPE": "TM",
                    "MGI_ALLELE_ACCESSION_ID": "MGI:900",
                    "ALLELE_SYMBOL": "Targeted<tm1>",
                }
            )
        )
        for index in range(20):
            writer.writerow(
                _row(
                    **{
                        "STRAIN/STOCK_ID": "MMRRC:000002-MU",
                        "MGI_GENE_ACCESSION_ID": f"MGI:{200 + index}",
                        "GENE_SYMBOL": f"Enu{index}",
                        "MUTATION_TYPE": "CI",
                    }
                )
            )
    loaded = MmrrcCatalog(path)
    loaded.load()
    return loaded


def test_targeted_stock_outranks_a_broader_enu_line(mixed_catalog: MmrrcCatalog) -> None:
    hits = mixed_catalog.stocks_for_genes(
        {"MGI:100": "Targeted", "MGI:200": "Enu0", "MGI:201": "Enu1"}
    )
    # The ENU line matches two genes to the knockout's one, but matches only
    # 2 of its own 20 annotations.
    assert [hit.stock_id for hit in hits] == ["MMRRC:000001-UNC", "MMRRC:000002-MU"]
    assert hits[0].matched_fraction == 1.0
    assert hits[0].annotated_gene_count == 1
    assert hits[1].annotated_gene_count == 20
    assert hits[1].matched_fraction == pytest.approx(0.1)


def test_excluding_a_mutation_type_drops_only_wholly_excluded_stocks(
    mixed_catalog: MmrrcCatalog,
) -> None:
    hits = mixed_catalog.stocks_for_genes(
        {"MGI:100": "Targeted", "MGI:200": "Enu0"},
        exclude_mutation_types=frozenset({"CI"}),
    )
    assert [hit.stock_id for hit in hits] == ["MMRRC:000001-UNC"]


def test_excluding_an_unrelated_mutation_type_keeps_everything(
    mixed_catalog: MmrrcCatalog,
) -> None:
    hits = mixed_catalog.stocks_for_genes(
        {"MGI:100": "Targeted", "MGI:200": "Enu0"},
        exclude_mutation_types=frozenset({"RAD"}),
    )
    assert len(hits) == 2


def test_allele_links_to_the_gene_its_symbol_names(tmp_path: Path) -> None:
    """The agouti case: the gene symbol is `a`, which alone looks like an error.

    The stock carries A<y> (agouti yellow), and the two live on different rows.
    """
    catalog = _load(
        tmp_path,
        [
            {
                "STRAIN/STOCK_ID": "MMRRC:000137-MU",
                "MGI_ALLELE_ACCESSION_ID": "MGI:1856798",
                "ALLELE_SYMBOL": "A<y>",
                "ALLELE_NAME": "nonagouti; agouti yellow",
            },
            {
                "STRAIN/STOCK_ID": "MMRRC:000137-MU",
                "MGI_GENE_ACCESSION_ID": "MGI:87853",
                "GENE_SYMBOL": "a",
            },
        ],
    )
    hits = catalog.stocks_for_genes({"MGI:87853": "a"})
    assert len(hits) == 1
    matched = hits[0].matched_alleles
    assert [a.symbol for a in matched] == ["A<y>"]
    assert matched[0].name == "nonagouti; agouti yellow"
    # Case-insensitive: allele `A<y>` names gene `a`.
    assert matched[0].link == "symbol"


def test_transgene_links_to_both_promoter_and_payload(tmp_path: Path) -> None:
    catalog = _load(
        tmp_path,
        [
            {
                "STRAIN/STOCK_ID": "MMRRC:000052-MU",
                "MGI_ALLELE_ACCESSION_ID": "MGI:2",
                "ALLELE_SYMBOL": "Tg(Myh6-Pln)11Egk",
            },
            {"STRAIN/STOCK_ID": "MMRRC:000052-MU", "MGI_GENE_ACCESSION_ID": "MGI:a", "GENE_SYMBOL": "Myh6"},
            {"STRAIN/STOCK_ID": "MMRRC:000052-MU", "MGI_GENE_ACCESSION_ID": "MGI:b", "GENE_SYMBOL": "Pln"},
        ],
    )
    allele = catalog.stocks["MMRRC:000052-MU"].alleles["MGI:2"]
    assert allele.gene_symbols == ["Myh6", "Pln"]
    assert allele.link == "symbol"
    # Querying either gene surfaces the same transgene.
    assert catalog.stocks_for_genes({"MGI:b": "Pln"})[0].matched_alleles[0].symbol == "Tg(Myh6-Pln)11Egk"


def test_a_lone_gene_links_by_fallback_and_is_labelled_as_such(tmp_path: Path) -> None:
    catalog = _load(
        tmp_path,
        [
            {
                "STRAIN/STOCK_ID": "MMRRC:3-MU",
                "MGI_ALLELE_ACCESSION_ID": "MGI:3",
                "ALLELE_SYMBOL": "Del(17Foo-Bar)1Xyz",
            },
            {"STRAIN/STOCK_ID": "MMRRC:3-MU", "MGI_GENE_ACCESSION_ID": "MGI:g", "GENE_SYMBOL": "Unrelated"},
        ],
    )
    allele = catalog.stocks["MMRRC:3-MU"].alleles["MGI:3"]
    assert allele.gene_symbols == ["Unrelated"]
    assert allele.link == "sole-gene"


def test_an_unlinkable_allele_on_a_multi_gene_stock_stays_unlinked(tmp_path: Path) -> None:
    catalog = _load(
        tmp_path,
        [
            {"STRAIN/STOCK_ID": "MMRRC:4-MU", "MGI_ALLELE_ACCESSION_ID": "MGI:4", "ALLELE_SYMBOL": "Tg(A94G6)1Rub"},
            {"STRAIN/STOCK_ID": "MMRRC:4-MU", "MGI_GENE_ACCESSION_ID": "MGI:x", "GENE_SYMBOL": "Csf2"},
            {"STRAIN/STOCK_ID": "MMRRC:4-MU", "MGI_GENE_ACCESSION_ID": "MGI:y", "GENE_SYMBOL": "Irf1"},
        ],
    )
    allele = catalog.stocks["MMRRC:4-MU"].alleles["MGI:4"]
    assert allele.gene_symbols == []
    assert allele.link == "none"
    # An unlinked allele must not be attributed to a matched gene.
    assert catalog.stocks_for_genes({"MGI:x": "Csf2"})[0].matched_alleles == []


def test_only_alleles_of_the_matched_gene_are_returned(tmp_path: Path) -> None:
    catalog = _load(
        tmp_path,
        [
            {"STRAIN/STOCK_ID": "MMRRC:5-MU", "MGI_ALLELE_ACCESSION_ID": "MGI:p", "ALLELE_SYMBOL": "Lep<ob>"},
            {"STRAIN/STOCK_ID": "MMRRC:5-MU", "MGI_ALLELE_ACCESSION_ID": "MGI:q", "ALLELE_SYMBOL": "Cd8a<tm1Mak>"},
            {"STRAIN/STOCK_ID": "MMRRC:5-MU", "MGI_GENE_ACCESSION_ID": "MGI:lep", "GENE_SYMBOL": "Lep"},
            {"STRAIN/STOCK_ID": "MMRRC:5-MU", "MGI_GENE_ACCESSION_ID": "MGI:cd8", "GENE_SYMBOL": "Cd8a"},
        ],
    )
    hit = catalog.stocks_for_genes({"MGI:lep": "Lep"})[0]
    assert [a.symbol for a in hit.matched_alleles] == ["Lep<ob>"]
    assert len(hit.alleles) == 2  # the full list is still available


def test_missing_required_column_is_an_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text("A,B\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing columns"):
        MmrrcCatalog(path).load()
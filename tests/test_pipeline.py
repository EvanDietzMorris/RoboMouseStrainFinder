"""Pipeline behaviour with the graph and name resolver stubbed.

The fixtures below encode the real shapes observed in RoboMouse KG: a human
MONDO term bridging to a mouse MP term through UPheno, human genes carrying
`NCBITaxon:9606`, and mouse orthologs carrying an `MGI:` equivalent identifier.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from robomouse_strain_finder.catalog import MmrrcCatalog
from robomouse_strain_finder.config import Settings
from robomouse_strain_finder.models import Candidate, Species
from robomouse_strain_finder.pipeline import StrainFinder

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


class FakeNameResolver:
    def __init__(self, candidates: list[Candidate]) -> None:
        self.candidates = candidates
        self.calls: list[tuple[str, Species, bool]] = []

    async def lookup(self, string, *, species=None, limit=10, autocomplete=False,
                     biolink_types=(), restrict_prefixes=True):
        self.calls.append((string, species, restrict_prefixes))
        return self.candidates


class FakeAutomat:
    """Dispatches on a distinctive fragment of each pipeline query."""

    def __init__(self, responses: dict[str, list[dict]]) -> None:
        self.responses = responses
        self.queries: list[str] = []

    async def cypher(self, query: str) -> list[dict]:
        self.queries.append(query)
        for marker, rows in self.responses.items():
            if marker in query:
                return rows
        return []

    async def node(self, curie: str):
        return {"id": curie, "name": "stub", "category": ["biolink:Disease"]}

    async def graph_version(self):
        return "RoboMouse KG 1.0.0"


@pytest.fixture
def catalog(tmp_path: Path) -> MmrrcCatalog:
    path = tmp_path / "catalog.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        row = dict.fromkeys(HEADER, "")
        row["STRAIN/STOCK_ID"] = "MMRRC:000012-UNC"
        row["STRAIN/STOCK_DESIGNATION"] = "B6.129S4-Ccr2/Mmnc"
        row["MGI_GENE_ACCESSION_ID"] = "MGI:106185"
        row["GENE_SYMBOL"] = "Ccr2"
        row["STATE"] = "CA"
        writer.writerow([row[column] for column in HEADER])
    loaded = MmrrcCatalog(path)
    loaded.load()
    return loaded


def build_finder(catalog: MmrrcCatalog, nameres, automat) -> StrainFinder:
    return StrainFinder(
        settings=Settings(mmrrc_catalog_path=catalog.path, auto_download_catalog=False),
        nameres=nameres,
        automat=automat,
        catalog=catalog,
    )


@pytest.mark.asyncio
async def test_human_query_reaches_a_mouse_strain_through_ortholog(catalog: MmrrcCatalog) -> None:
    nameres = FakeNameResolver(
        [Candidate(curie="MONDO:0005027", label="epilepsy", score=6386.0,
                   types=["biolink:Disease"])]
    )
    automat = FakeAutomat(
        {
            "homologous_to": [
                {"curie": "MP:0002064", "label": "seizures",
                 "labels": ["biolink:PhenotypicFeature"], "source": "MONDO:0005027"}
            ],
            "biolink:Gene`)\n        WHERE t.id IN": [
                {
                    "curie": "NCBIGene:729230", "symbol": "CCR2", "taxon": "NCBITaxon:9606",
                    "predicates": ["biolink:genetically_associated_with"],
                    "seed_curies": ["MONDO:0005027"], "sources": ["infores:disgenet"],
                    "mgi_ids": [], "seed_count": 1, "predicate_count": 1,
                }
            ],
            "orthologous_to": [
                {
                    "curie": "NCBIGene:12772", "symbol": "Ccr2", "taxon": "NCBITaxon:10090",
                    "mgi_ids": ["MGI:106185"], "ortholog_of": ["NCBIGene:729230"],
                }
            ],
            "biolink:Pathway": [
                {
                    "curie": "GO:0006955", "name": "immune response",
                    "gene_curies": ["NCBIGene:729230"], "gene_symbols": ["CCR2"],
                    "predicates": ["biolink:actively_involved_in"], "gene_count": 1,
                }
            ],
        }
    )
    finder = build_finder(catalog, nameres, automat)

    result = await finder.search(term="epilepsy", species=Species.HUMAN)

    assert result.selected is not None
    assert result.selected.curie == "MONDO:0005027"
    # UPheno bridge added the mouse phenotype.
    assert [seed.curie for seed in result.seed_terms] == ["MONDO:0005027", "MP:0002064"]
    assert result.counts.human_genes == 1
    assert result.counts.mouse_genes == 1

    mouse_gene = next(gene for gene in result.genes if gene.species is Species.MOUSE)
    assert mouse_gene.via == "ortholog"
    assert mouse_gene.ortholog_of == "NCBIGene:729230"
    # Seed provenance carries across the orthology hop.
    assert mouse_gene.seed_curies == ["MONDO:0005027"]

    assert [strain.stock_id for strain in result.strains] == ["MMRRC:000012-UNC"]
    assert result.strains[0].matched_mgi_gene_ids == ["MGI:106185"]
    assert result.pathways[0].curie == "GO:0006955"


@pytest.mark.asyncio
async def test_orthologs_can_be_disabled(catalog: MmrrcCatalog) -> None:
    nameres = FakeNameResolver(
        [Candidate(curie="MONDO:0005027", label="epilepsy", score=1.0)]
    )
    automat = FakeAutomat(
        {
            "biolink:Gene`)\n        WHERE t.id IN": [
                {
                    "curie": "NCBIGene:729230", "symbol": "CCR2", "taxon": "NCBITaxon:9606",
                    "predicates": ["biolink:genetically_associated_with"],
                    "seed_curies": ["MONDO:0005027"], "sources": [], "mgi_ids": [],
                    "seed_count": 1, "predicate_count": 1,
                }
            ],
            "orthologous_to": [
                {"curie": "NCBIGene:12772", "symbol": "Ccr2", "taxon": "NCBITaxon:10090",
                 "mgi_ids": ["MGI:106185"], "ortholog_of": ["NCBIGene:729230"]}
            ],
        }
    )
    finder = build_finder(catalog, nameres, automat)

    result = await finder.search(term="epilepsy", species=Species.HUMAN, include_orthologs=False)

    assert result.counts.mouse_genes == 0
    assert result.strains == []
    assert not any("orthologous_to" in query for query in automat.queries)


@pytest.mark.asyncio
async def test_prefix_filter_falls_back_when_species_vocabulary_misses(
    catalog: MmrrcCatalog,
) -> None:
    class PickyResolver(FakeNameResolver):
        async def lookup(self, string, *, species=None, limit=10, autocomplete=False,
                         biolink_types=(), restrict_prefixes=True):
            self.calls.append((string, species, restrict_prefixes))
            if restrict_prefixes:
                return []
            return [Candidate(curie="MONDO:0005027", label="epilepsy", score=1.0)]

    nameres = PickyResolver([])
    finder = build_finder(catalog, nameres, FakeAutomat({}))

    result = await finder.search(term="epilepsy", species=Species.MOUSE)

    assert result.selected is not None
    assert [call[2] for call in nameres.calls] == [True, False]
    assert any("fell back" in note for note in result.notes)


@pytest.mark.asyncio
async def test_unresolvable_term_returns_an_empty_result_with_a_note(
    catalog: MmrrcCatalog,
) -> None:
    finder = build_finder(catalog, FakeNameResolver([]), FakeAutomat({}))

    result = await finder.search(term="not a real phenotype", species=Species.MOUSE)

    assert result.selected is None
    assert result.genes == []
    assert result.strains == []
    assert result.counts.strains == 0
    assert any("no mouse disease/phenotype match" in note for note in result.notes)


@pytest.mark.asyncio
async def test_mouse_genes_without_mgi_accession_are_reported(catalog: MmrrcCatalog) -> None:
    nameres = FakeNameResolver([Candidate(curie="MP:0002064", label="seizures", score=1.0)])
    automat = FakeAutomat(
        {
            "biolink:Gene`)\n        WHERE t.id IN": [
                {
                    "curie": "NCBIGene:99999", "symbol": "Xyz", "taxon": "NCBITaxon:10090",
                    "predicates": ["biolink:has_phenotype"], "seed_curies": ["MP:0002064"],
                    "sources": [], "mgi_ids": [], "seed_count": 1, "predicate_count": 1,
                }
            ]
        }
    )
    finder = build_finder(catalog, nameres, automat)

    result = await finder.search(term="seizures", species=Species.MOUSE)

    assert result.counts.mouse_genes == 1
    assert result.strains == []
    assert any("none carried an MGI accession" in note for note in result.notes)
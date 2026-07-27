"""Response schemas.

The shapes here deliberately keep repository stock, allele, gene, and species
as separate fields. A gene-level link between a phenotype and an MMRRC stock is
a *candidate*, not evidence that the stock's allele models the phenotype, and
the schema should not let those collapse into each other.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, Field

MOUSE_TAXON = "NCBITaxon:10090"
HUMAN_TAXON = "NCBITaxon:9606"


class Species(StrEnum):
    MOUSE = "mouse"
    HUMAN = "human"

    @property
    def taxon(self) -> str:
        return MOUSE_TAXON if self is Species.MOUSE else HUMAN_TAXON

    @property
    def other(self) -> "Species":
        return Species.HUMAN if self is Species.MOUSE else Species.MOUSE

    @property
    def prefixes(self) -> tuple[str, ...]:
        """Ontology prefixes that carry this species' phenotype vocabulary."""
        if self is Species.MOUSE:
            return ("MP",)
        return ("MONDO", "HP")


class Candidate(BaseModel):
    """One Name Resolver hit for a user-supplied phenotype or disease string."""

    curie: str
    label: str
    score: float
    types: list[str] = Field(default_factory=list)
    synonyms: list[str] = Field(default_factory=list)
    taxa: list[str] = Field(default_factory=list)


class ResolveResponse(BaseModel):
    query: str
    species: Species
    candidates: list[Candidate]


class SeedTerm(BaseModel):
    """A phenotype/disease node used as a traversal starting point.

    `via` records how the term entered the seed set: `query` for the term the
    user selected, `upheno` for one reached through UPheno cross-species
    phenotype homology.
    """

    curie: str
    label: str | None = None
    category: list[str] = Field(default_factory=list)
    via: str = "query"
    homologous_to: str | None = None


class GeneHit(BaseModel):
    curie: str
    symbol: str | None = None
    taxon: str | None = None
    species: Species | None = None
    mgi_ids: list[str] = Field(default_factory=list)
    # How this gene was reached: direct edge from a seed term, or orthology.
    via: str = "direct"
    predicates: list[str] = Field(default_factory=list)
    seed_curies: list[str] = Field(default_factory=list)
    ortholog_of: str | None = None
    knowledge_sources: list[str] = Field(default_factory=list)


class PathwayHit(BaseModel):
    curie: str
    name: str | None = None
    predicates: list[str] = Field(default_factory=list)
    gene_curies: list[str] = Field(default_factory=list)
    gene_symbols: list[str] = Field(default_factory=list)


class StrainAllele(BaseModel):
    """An allele carried by a stock.

    The catalog never states which gene an allele belongs to, so `gene_symbols`
    is recovered from MGI naming and `link` records how confident that is:

    - `symbol` -- the allele symbol names the gene (`Esr2<tm1Unc>` -> `Esr2`,
      `Tg(Myh6-Pln)11Egk` -> `Myh6`, `Pln`). ~72% of alleles.
    - `sole-gene` -- inferred only because the stock annotates a single gene.
      ~28% of alleles.
    - `none` -- could not be linked (~0.2%).
    """

    mgi_allele_id: str | None = None
    symbol: str | None = None
    name: str | None = None
    gene_symbols: list[str] = Field(default_factory=list)
    link: str = "none"


class GeneEvidence(BaseModel):
    """Why a gene on a stock was implicated in the queried phenotype.

    Carried on the strain so a hit can be judged without cross-referencing the
    gene table. `via` matters for reading it honestly: on an `ortholog` row the
    predicates and sources describe the **human** gene named in `ortholog_of`,
    not the mouse gene the stock actually carries.
    """

    gene_symbol: str | None = None
    gene_curie: str | None = None
    via: str = "direct"
    ortholog_of: str | None = None
    ortholog_of_symbol: str | None = None
    predicates: list[str] = Field(default_factory=list)
    knowledge_sources: list[str] = Field(default_factory=list)
    seed_curies: list[str] = Field(default_factory=list)


class StrainHit(BaseModel):
    """An MMRRC stock reached from the gene set.

    `matched_mgi_gene_ids` is the join key that produced this stock. It is a
    gene-level link: it does not assert that the stock's allele reproduces the
    queried phenotype.

    `annotated_gene_count` is what separates a targeted mutant from an ENU
    mutagenesis line. A knockout carries one or two annotated genes; a
    chemically induced Mutagenetix stock can carry dozens of incidental
    variants, so a raw match count would rank the noisy stock first.
    `matched_fraction` is that ratio, and is what the ranking sorts on.
    """

    stock_id: str
    rrids: list[str] = Field(default_factory=list)
    designation: str | None = None
    strain_types: list[str] = Field(default_factory=list)
    states: list[str] = Field(default_factory=list)
    mutation_types: list[str] = Field(default_factory=list)
    alleles: list[StrainAllele] = Field(default_factory=list)
    #: The subset of `alleles` linked to a gene that matched the query. This is
    #: what makes a hit legible: the gene may be `a`, but the allele is `A<y>`.
    matched_alleles: list[StrainAllele] = Field(default_factory=list)
    #: The gene->phenotype evidence behind each matched gene, joined on from the
    #: gene stage so the strain row is self-contained.
    gene_evidence: list[GeneEvidence] = Field(default_factory=list)
    gene_symbols: list[str] = Field(default_factory=list)
    matched_mgi_gene_ids: list[str] = Field(default_factory=list)
    matched_gene_symbols: list[str] = Field(default_factory=list)
    annotated_gene_count: int = 0
    matched_fraction: float = 0.0
    #: True when every allele on the stock is a reporter or recombinase
    #: transgene -- a GENSAT `Tg(gene-EGFP)` line or a cre driver. These carry
    #: one gene and so score a perfect gene-match ratio, but they leave the gene
    #: intact and are research tools rather than models of its disease.
    tool_line: bool = False
    phenotypes: list[str] = Field(default_factory=list)
    pubmed_ids: list[str] = Field(default_factory=list)
    research_areas: list[str] = Field(default_factory=list)
    sds_url: str | None = None
    accepted_date: str | None = None


class Counts(BaseModel):
    seed_terms: int = 0
    genes: int = 0
    mouse_genes: int = 0
    human_genes: int = 0
    genes_with_mgi: int = 0
    pathways: int = 0
    strains: int = 0


class Provenance(BaseModel):
    name_resolver_url: str
    graph_url: str
    graph_version: str | None = None
    catalog_path: str
    catalog_modified_at: str | None = None
    catalog_rows: int | None = None
    retrieved_at: str


class SearchResponse(BaseModel):
    query: str
    species: Species
    selected: Candidate | None = None
    seed_terms: list[SeedTerm] = Field(default_factory=list)
    genes: list[GeneHit] = Field(default_factory=list)
    pathways: list[PathwayHit] = Field(default_factory=list)
    strains: list[StrainHit] = Field(default_factory=list)
    counts: Counts = Field(default_factory=Counts)
    provenance: Provenance
    notes: list[str] = Field(default_factory=list)
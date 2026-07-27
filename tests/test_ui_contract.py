"""Pins the JSON field names the browser client reads.

`static/app.js` addresses response fields by name. If a model field is renamed,
FastAPI keeps returning valid JSON and the UI silently renders blanks -- nothing
else in the suite would catch it. These lists mirror the accesses in app.js.
"""

from __future__ import annotations

import pytest

from robomouse_strain_finder.models import (
    Candidate,
    Counts,
    GeneEvidence,
    GeneHit,
    PathwayHit,
    Provenance,
    SearchResponse,
    SeedTerm,
    StrainAllele,
    StrainHit,
)

EXPECTED = {
    SearchResponse: {
        "query", "species", "selected", "seed_terms", "genes",
        "pathways", "strains", "counts", "provenance", "notes",
    },
    Counts: {
        "seed_terms", "genes", "mouse_genes", "human_genes",
        "genes_with_mgi", "pathways", "strains",
    },
    Candidate: {"curie", "label"},
    SeedTerm: {"curie", "label", "via", "homologous_to"},
    GeneHit: {
        "curie", "symbol", "species", "taxon", "via", "ortholog_of",
        "predicates", "mgi_ids", "knowledge_sources",
    },
    PathwayHit: {"curie", "name", "gene_curies", "gene_symbols"},
    StrainHit: {
        "stock_id", "rrids", "designation", "matched_gene_symbols",
        "matched_mgi_gene_ids", "matched_fraction", "annotated_gene_count", "tool_line",
        "alleles", "matched_alleles", "gene_evidence", "mutation_types", "states", "strain_types",
        "phenotypes", "pubmed_ids", "research_areas", "accepted_date", "sds_url",
    },
    StrainAllele: {"symbol", "mgi_allele_id", "name", "link", "gene_symbols"},
    GeneEvidence: {
        "gene_symbol", "gene_curie", "via", "ortholog_of", "ortholog_of_symbol",
        "predicates", "knowledge_sources", "seed_curies",
    },
    Provenance: {
        "name_resolver_url", "graph_url", "graph_version", "catalog_path",
        "catalog_rows", "catalog_modified_at", "retrieved_at",
    },
}


@pytest.mark.parametrize("model,fields", EXPECTED.items(), ids=lambda v: getattr(v, "__name__", ""))
def test_ui_fields_exist_on_the_model(model, fields) -> None:
    missing = fields - set(model.model_fields)
    assert not missing, f"{model.__name__} no longer provides {sorted(missing)} — static/app.js reads it"
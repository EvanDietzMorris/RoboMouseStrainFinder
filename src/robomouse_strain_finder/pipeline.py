"""Phenotype/disease -> genes -> orthologs -> pathways -> MMRRC strains.

The traversal is deliberately staged so each hop stays inspectable and bounded:

1. **Seed terms** -- the resolved CURIE, plus cross-species phenotype homologs
   from UPheno (``biolink:homologous_to``: HP<->MP, and MONDO->MP). This is what
   makes a human disease query able to reach mouse biology at all.
2. **Genes** -- gene nodes adjacent to any seed via a curated set of
   association predicates. Ranked by how many distinct seeds and predicates
   support them.
3. **Orthologs** -- Alliance orthology (``biolink:orthologous_to``) used to
   reach mouse genes from human genes and vice versa.
4. **Pathways** -- pathway membership for the assembled gene set.
5. **Strains** -- mouse genes carry MGI accessions in their equivalent
   identifiers; those join to the MMRRC catalog.

Gene hub nodes can have thousands of edges, so every stage caps its own result
set rather than relying on a single global limit.
"""

from __future__ import annotations

import datetime as dt

from .catalog import MmrrcCatalog
from .clients import AutomatClient, NameResolverClient
from .clients.automat import curie_list_literal, validate_curie
from .config import Settings
from .models import (
    Candidate,
    Counts,
    GeneEvidence,
    GeneHit,
    PathwayHit,
    Provenance,
    SearchResponse,
    SeedTerm,
    Species,
    StrainHit,
)

#: Predicates treated as gene<->condition association evidence.
#:
#: Excludes `biolink:related_to` and `biolink:target_for`, which in ROBOKOP
#: carry Hetionet co-occurrence and drug-target semantics respectively -- neither
#: means the gene is implicated in the condition.
ASSOCIATION_PREDICATES = (
    "biolink:has_phenotype",
    "biolink:genetically_associated_with",
    "biolink:gene_associated_with_condition",
    "biolink:causes",
    "biolink:contributes_to",
    "biolink:model_of",
    "biolink:biomarker_for",
)

PATHWAY_PREDICATES = (
    "biolink:actively_involved_in",
    "biolink:affects",
    "biolink:acts_upstream_of",
    "biolink:acts_upstream_of_positive_effect",
    "biolink:acts_upstream_of_negative_effect",
    "biolink:regulates",
)

#: Node labels worth surfacing; Neo4j labels include the whole Biolink ancestry.
INTERESTING_LABELS = frozenset(
    {
        "biolink:Disease",
        "biolink:PhenotypicFeature",
        "biolink:BehavioralFeature",
        "biolink:ClinicalFinding",
    }
)

TAXON_TO_SPECIES = {
    Species.MOUSE.taxon: Species.MOUSE,
    Species.HUMAN.taxon: Species.HUMAN,
}


def _predicate_list_literal(predicates: tuple[str, ...]) -> str:
    return "[" + ", ".join(f'"{predicate}"' for predicate in predicates) + "]"


def _tidy_labels(labels: list[str] | None) -> list[str]:
    return sorted(set(labels or []) & INTERESTING_LABELS)


class StrainFinder:
    def __init__(
        self,
        *,
        settings: Settings,
        nameres: NameResolverClient,
        automat: AutomatClient,
        catalog: MmrrcCatalog,
    ) -> None:
        self._settings = settings
        self._nameres = nameres
        self._automat = automat
        self._catalog = catalog

    async def resolve(
        self, term: str, species: Species, *, limit: int = 10, restrict_prefixes: bool = True
    ) -> list[Candidate]:
        return await self._nameres.lookup(
            term, species=species, limit=limit, restrict_prefixes=restrict_prefixes
        )

    async def search(
        self,
        *,
        term: str,
        species: Species,
        curie: str | None = None,
        max_genes: int = 100,
        max_pathways: int = 40,
        max_strains: int = 200,
        include_orthologs: bool = True,
        bridge_species: bool = True,
        exclude_mutation_types: frozenset[str] = frozenset(),
        exclude_tool_lines: bool = False,
    ) -> SearchResponse:
        notes: list[str] = []

        selected, candidates = await self._select_term(term, species, curie, notes)
        if selected is None:
            return self._empty_response(term, species, notes)

        seeds = await self._seed_terms(selected, species, bridge_species, notes)
        genes = await self._genes_for_seeds(seeds, max_genes)
        if not genes:
            notes.append(
                "No genes were associated with the seed terms through the curated "
                "association predicates."
            )

        if include_orthologs:
            genes = await self._add_orthologs(genes, max_genes)

        mouse_genes = [gene for gene in genes if gene.species is Species.MOUSE]
        human_genes = [gene for gene in genes if gene.species is Species.HUMAN]

        pathways = await self._pathways(genes, max_pathways)
        strains, gene_index = self._strains(
            mouse_genes,
            max_strains,
            notes,
            exclude_mutation_types,
            exclude_tool_lines,
            gene_by_curie={gene.curie: gene for gene in genes},
        )

        counts = Counts(
            seed_terms=len(seeds),
            genes=len(genes),
            mouse_genes=len(mouse_genes),
            human_genes=len(human_genes),
            genes_with_mgi=len(gene_index),
            pathways=len(pathways),
            strains=len(strains),
        )
        return SearchResponse(
            query=term,
            species=species,
            selected=selected,
            seed_terms=seeds,
            genes=genes,
            pathways=pathways,
            strains=strains,
            counts=counts,
            provenance=await self._provenance(),
            notes=notes + self._candidate_note(selected, candidates),
        )

    # -- stages ---------------------------------------------------------------

    async def _select_term(
        self, term: str, species: Species, curie: str | None, notes: list[str]
    ) -> tuple[Candidate | None, list[Candidate]]:
        """Pick the term to traverse from: an explicit CURIE, or the top hit."""
        if curie:
            validate_curie(curie)
            node = await self._automat.node(curie)
            if node is None:
                notes.append(
                    f"{curie} was not found in {self._settings.graph}; it may be valid "
                    "in another graph build or resolve to a different canonical identifier."
                )
                return None, []
            return (
                Candidate(
                    curie=node.get("id", curie),
                    label=node.get("name") or curie,
                    score=0.0,
                    types=_tidy_labels(node.get("category")),
                ),
                [],
            )

        candidates = await self.resolve(term, species)
        if not candidates:
            # Species prefix filters can be too tight for an unusual term.
            candidates = await self.resolve(term, species, restrict_prefixes=False)
            if candidates:
                notes.append(
                    f"No {'/'.join(species.prefixes)} term matched {term!r}; fell back to "
                    "the unrestricted vocabulary."
                )
        if not candidates:
            notes.append(f"Name Resolver returned no {species} disease/phenotype match for {term!r}.")
            return None, []
        return candidates[0], candidates

    async def _seed_terms(
        self, selected: Candidate, species: Species, bridge: bool, notes: list[str]
    ) -> list[SeedTerm]:
        seeds = [
            SeedTerm(
                curie=selected.curie,
                label=selected.label,
                category=selected.types,
                via="query",
            )
        ]
        if not bridge:
            return seeds

        query = f"""
        MATCH (t)-[:`biolink:homologous_to`]-(h)
        WHERE t.id IN {curie_list_literal([selected.curie])}
          AND (h:`biolink:PhenotypicFeature` OR h:`biolink:Disease`)
        RETURN DISTINCT h.id AS curie, h.name AS label, labels(h) AS labels, t.id AS source
        ORDER BY curie
        LIMIT 50
        """
        rows = await self._automat.cypher(query)
        for row in rows:
            seeds.append(
                SeedTerm(
                    curie=row["curie"],
                    label=row.get("label"),
                    category=_tidy_labels(row.get("labels")),
                    via="upheno",
                    homologous_to=row.get("source"),
                )
            )
        if len(seeds) == 1:
            notes.append(
                f"{selected.curie} has no UPheno cross-species phenotype homolog in "
                f"{self._settings.graph}; only direct gene associations were used."
            )
        return seeds

    async def _genes_for_seeds(self, seeds: list[SeedTerm], max_genes: int) -> list[GeneHit]:
        query = f"""
        MATCH (t)-[e]-(g:`biolink:Gene`)
        WHERE t.id IN {curie_list_literal([seed.curie for seed in seeds])}
          AND type(e) IN {_predicate_list_literal(ASSOCIATION_PREDICATES)}
        WITH g,
             collect(DISTINCT type(e)) AS predicates,
             collect(DISTINCT t.id) AS seed_curies,
             collect(DISTINCT e.primary_knowledge_source) AS sources
        RETURN g.id AS curie, g.name AS symbol, g.taxon AS taxon,
               predicates, seed_curies, sources,
               [x IN g.equivalent_identifiers WHERE x STARTS WITH 'MGI:'] AS mgi_ids,
               size(seed_curies) AS seed_count, size(predicates) AS predicate_count
        ORDER BY seed_count DESC, predicate_count DESC, curie
        LIMIT {int(max_genes)}
        """
        rows = await self._automat.cypher(query)
        return [
            GeneHit(
                curie=row["curie"],
                symbol=row.get("symbol"),
                taxon=row.get("taxon"),
                species=TAXON_TO_SPECIES.get(row.get("taxon") or ""),
                mgi_ids=row.get("mgi_ids") or [],
                via="direct",
                predicates=sorted(row.get("predicates") or []),
                seed_curies=sorted(row.get("seed_curies") or []),
                knowledge_sources=sorted(x for x in (row.get("sources") or []) if x),
            )
            for row in rows
        ]

    async def _add_orthologs(self, genes: list[GeneHit], max_genes: int) -> list[GeneHit]:
        """Expand the gene set across species using Alliance orthology.

        Mouse orthologs of human hits are what ultimately reach MMRRC stocks;
        human orthologs of mouse hits are kept so the human-facing context
        stays visible.
        """
        if not genes:
            return genes
        known = {gene.curie for gene in genes}
        query = f"""
        MATCH (g)-[:`biolink:orthologous_to`]-(o:`biolink:Gene`)
        WHERE g.id IN {curie_list_literal([gene.curie for gene in genes])}
          AND o.taxon IN ["{Species.MOUSE.taxon}", "{Species.HUMAN.taxon}"]
        RETURN DISTINCT o.id AS curie, o.name AS symbol, o.taxon AS taxon,
               [x IN o.equivalent_identifiers WHERE x STARTS WITH 'MGI:'] AS mgi_ids,
               collect(DISTINCT g.id) AS ortholog_of
        // Mouse first (only mouse reaches a stock), then a stable tiebreak.
        // A bare LIMIT would truncate arbitrarily and give different orthologs
        // run to run -- a promiscuous family like SYCP3 has 76 of them.
        ORDER BY (o.taxon = "{Species.MOUSE.taxon}") DESC, size(ortholog_of) DESC, curie
        LIMIT {int(max_genes) * 4}
        """
        rows = await self._automat.cypher(query)

        seeds_by_gene = {gene.curie: gene.seed_curies for gene in genes}
        added: list[GeneHit] = []
        for row in rows:
            curie = row["curie"]
            if curie in known:
                continue
            known.add(curie)
            sources = row.get("ortholog_of") or []
            added.append(
                GeneHit(
                    curie=curie,
                    symbol=row.get("symbol"),
                    taxon=row.get("taxon"),
                    species=TAXON_TO_SPECIES.get(row.get("taxon") or ""),
                    mgi_ids=row.get("mgi_ids") or [],
                    via="ortholog",
                    ortholog_of=sources[0] if sources else None,
                    seed_curies=sorted(
                        {
                            seed
                            for source in sources
                            for seed in seeds_by_gene.get(source, ())
                        }
                    ),
                )
            )
        # Mouse orthologs first: they are the ones that can reach a stock.
        added.sort(key=lambda gene: (gene.species is not Species.MOUSE, gene.curie))
        return genes + added

    async def _pathways(self, genes: list[GeneHit], max_pathways: int) -> list[PathwayHit]:
        if not genes:
            return []
        query = f"""
        MATCH (g)-[e]-(p:`biolink:Pathway`)
        WHERE g.id IN {curie_list_literal([gene.curie for gene in genes])}
          AND type(e) IN {_predicate_list_literal(PATHWAY_PREDICATES)}
        WITH p,
             collect(DISTINCT g.id) AS gene_curies,
             collect(DISTINCT g.name) AS gene_symbols,
             collect(DISTINCT type(e)) AS predicates
        RETURN p.id AS curie, p.name AS name, gene_curies, gene_symbols, predicates,
               size(gene_curies) AS gene_count
        ORDER BY gene_count DESC, curie
        LIMIT {int(max_pathways)}
        """
        rows = await self._automat.cypher(query)
        return [
            PathwayHit(
                curie=row["curie"],
                name=row.get("name"),
                predicates=sorted(row.get("predicates") or []),
                gene_curies=sorted(row.get("gene_curies") or []),
                gene_symbols=sorted(x for x in (row.get("gene_symbols") or []) if x),
            )
            for row in rows
        ]

    def _strains(
        self,
        mouse_genes: list[GeneHit],
        max_strains: int,
        notes: list[str],
        exclude_mutation_types: frozenset[str] = frozenset(),
        exclude_tool_lines: bool = False,
        gene_by_curie: dict[str, GeneHit] | None = None,
    ) -> tuple[list[StrainHit], dict[str, str]]:
        gene_by_curie = gene_by_curie or {}
        gene_index: dict[str, str] = {}
        gene_by_mgi: dict[str, GeneHit] = {}
        for gene in mouse_genes:
            for mgi_id in gene.mgi_ids:
                gene_index[mgi_id] = gene.symbol or gene.curie
                gene_by_mgi[mgi_id] = gene

        if not gene_index:
            if mouse_genes:
                notes.append(
                    "Mouse genes were found but none carried an MGI accession in their "
                    "equivalent identifiers, so no MMRRC join was possible."
                )
            return [], gene_index

        hits = self._catalog.stocks_for_genes(
            gene_index,
            exclude_mutation_types=exclude_mutation_types,
            exclude_tool_lines=exclude_tool_lines,
        )
        if hits and hits[0].matched_fraction < 0.5:
            notes.append(
                "No stock matched on a majority of its annotated genes. Top hits are "
                "likely mutagenesis lines carrying the query gene as one of many "
                "incidental variants rather than targeted models."
            )
        if len(hits) > max_strains:
            notes.append(
                f"{len(hits)} MMRRC stocks matched; showing the first {max_strains}. "
                "Narrow the query or raise max_strains to see the rest."
            )
            hits = hits[:max_strains]

        # Attach the gene->phenotype evidence to each hit. No extra queries --
        # this is the gene stage's own output, joined on the MGI accession that
        # produced the match.
        for hit in hits:
            hit.gene_evidence = [
                self._evidence_for(gene_by_mgi[mgi_id], gene_by_curie)
                for mgi_id in hit.matched_mgi_gene_ids
                if mgi_id in gene_by_mgi
            ]
        return hits, gene_index

    def _evidence_for(self, gene: GeneHit, gene_by_curie: dict[str, GeneHit]) -> GeneEvidence:
        """Build the evidence for one matched gene.

        A gene reached by orthology has no association of its own -- the
        predicates and sources belong to its human partner, so they are read
        from that gene and reported with `via="ortholog"` so the attribution
        stays explicit rather than being silently transferred to the mouse gene.
        """
        source = gene
        partner = gene_by_curie.get(gene.ortholog_of) if gene.ortholog_of else None
        if gene.via == "ortholog" and partner is not None:
            source = partner
        return GeneEvidence(
            gene_symbol=gene.symbol,
            gene_curie=gene.curie,
            via=gene.via,
            ortholog_of=gene.ortholog_of,
            # None when the partner fell outside the returned gene set; the UI
            # then shows the CURIE rather than inventing a symbol.
            ortholog_of_symbol=partner.symbol if partner else None,
            predicates=source.predicates,
            knowledge_sources=source.knowledge_sources,
            seed_curies=source.seed_curies or gene.seed_curies,
        )

    # -- helpers --------------------------------------------------------------

    def _candidate_note(self, selected: Candidate, candidates: list[Candidate]) -> list[str]:
        others = [c for c in candidates if c.curie != selected.curie][:4]
        if not others:
            return []
        alternatives = ", ".join(f"{c.label} ({c.curie})" for c in others)
        return [
            f"Traversed {selected.label} ({selected.curie}). Other name matches: {alternatives}."
        ]

    async def _provenance(self) -> Provenance:
        return Provenance(
            name_resolver_url=self._settings.name_resolver_url,
            graph_url=self._settings.graph_url,
            graph_version=await self._automat.graph_version(),
            catalog_path=str(self._catalog.path),
            catalog_modified_at=self._catalog.modified_at,
            catalog_rows=self._catalog.rows,
            retrieved_at=dt.datetime.now(dt.UTC).isoformat(),
        )

    def _empty_response(self, term: str, species: Species, notes: list[str]) -> SearchResponse:
        return SearchResponse(
            query=term,
            species=species,
            counts=Counts(),
            provenance=Provenance(
                name_resolver_url=self._settings.name_resolver_url,
                graph_url=self._settings.graph_url,
                catalog_path=str(self._catalog.path),
                catalog_modified_at=self._catalog.modified_at,
                catalog_rows=self._catalog.rows,
                retrieved_at=dt.datetime.now(dt.UTC).isoformat(),
            ),
            notes=notes,
        )
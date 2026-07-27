# RoboMouse Strain Finder — how it works

**The question it answers:** *"I study X. Is there a mouse I can order?"*

A researcher types a disease or phenotype in plain English. The app returns
MMRRC mouse strains that are genetically relevant to it — with the gene, the
allele, and the evidence trail that connects them.

Nothing about this is a lookup. The catalog has no "disease" column. Getting
from *"Rett syndrome"* to a shippable stock number requires crossing four
independent resources, two species, and five identifier systems.

```
  "Rett syndrome"
        │
        │  1. Name Resolver            text ──▶ CURIE
        ▼
   MONDO:0010726
        │
        │  2. UPheno                   human phenotype ◀──▶ mouse phenotype
        │  3. RoboMouse KG             condition ──▶ genes
        ▼
   MECP2 (human)  ──── 4. Alliance orthology ────▶  Mecp2 (mouse)
        │                                                │
        │  5. GO / Reactome                              │  MGI:99918
        ▼                                                ▼
    pathways                                    6. MMRRC catalog join
                                                         │
                                                         ▼
                                          MMRRC:000011-UCD  Mecp2<tm1.1Jae>
```

---

## The resources

| Resource | Role | Endpoint |
|---|---|---|
| **NameResolution** | free text → normalized CURIE | `name-resolution-exp.apps.renci.org` |
| **RoboMouse KG** (Automat) | the knowledge graph | `robokop-automat.apps.renci.org/robomousekg` |
| **MMRRC bulk catalog** | the orderable inventory | `mmrrc.org/about/mmrrc_catalog_data.csv` |

RoboMouse KG is ROBOKOP plus the mouse layer that makes this work: mouse GO
annotations, **Alliance of Genome Resources orthology**, **UPheno human/mouse
phenotype homology**, and MGI phenotype/disease associations.

The catalog is downloaded once and held in memory:

```
588,980 rows  →  69,388 stocks  →  23,738 distinct MGI gene accessions indexed
```

Rows are not stocks. A single stock is spread across many rows — its allele
accession and its gene accession are on *different* rows. Counting rows as
stocks is the first mistake available in this dataset.

---

## Worked example: "Rett syndrome"

Every query and every result below is a real call against the live services.

### Step 1 — Resolve the text to an identifier

`GET /lookup` on NameResolution. Disease and phenotype terms are not
taxon-tagged the way genes are, so the mouse/human selector is applied as an
**ontology prefix filter**, not a taxon filter:

```http
POST https://name-resolution-exp.apps.renci.org/lookup
  ?string=Rett syndrome
  &biolink_type=biolink:Disease
  &biolink_type=biolink:PhenotypicFeature
  &only_prefixes=MONDO|HP          # "human"; "mouse" sends only_prefixes=MP
  &limit=5
```

```
MONDO:0010726   Rett syndrome                    score=5993
MONDO:0017746   atypical Rett syndrome           score=1060
MONDO:0100040   FOXG1 disorder                   score=462
MONDO:7770183   Rett syndrome, rhesus macaque    score=354
```

The top hit is used, and the alternatives are reported back so the user can see
what else matched. In the UI, picking from the type-ahead **pins** that CURIE so
the traversal uses exactly the term chosen rather than the top-scoring guess.

### Step 2 — Bridge species through UPheno

Human and mouse describe phenotypes in different vocabularies (HP vs MP).
UPheno's `biolink:homologous_to` edges connect them, so a human query can reach
mouse biology.

```cypher
MATCH (t)-[:`biolink:homologous_to`]-(h)
WHERE t.id IN ["MONDO:0010726"]
  AND (h:`biolink:PhenotypicFeature` OR h:`biolink:Disease`)
RETURN DISTINCT h.id AS curie, h.name AS label
```

```
(no homolog for this term)
```

**This one returns nothing, and that's worth showing.** Rett syndrome has no
UPheno phenotype homolog in this build. The app says so in a note rather than
failing silently — and the traversal still succeeds, because MGI's `model_of`
edges reach mouse genes directly in the next step.

For a term where it *does* fire, `hearing loss disorder` (`MONDO:0005365`)
bridges to `MP:0006325 impaired hearing`, pulling in mouse genes that the human
disease term alone would never reach.

### Step 3 — Condition → genes

The traversal is restricted to a **curated set of association predicates**:

```cypher
MATCH (t)-[e]-(g:`biolink:Gene`)
WHERE t.id IN ["MONDO:0010726"]
  AND type(e) IN ["biolink:has_phenotype", "biolink:genetically_associated_with",
                  "biolink:gene_associated_with_condition", "biolink:causes",
                  "biolink:contributes_to", "biolink:model_of", "biolink:biomarker_for"]
WITH g, collect(DISTINCT type(e)) AS predicates,
        collect(DISTINCT t.id) AS seed_curies,
        collect(DISTINCT e.primary_knowledge_source) AS sources
RETURN g.id, g.name, g.taxon, predicates, sources
ORDER BY size(seed_curies) DESC, size(predicates) DESC
```

```
NCBIGene:4204    MECP2   NCBITaxon:9606    causes, genetically_associated_with,
                                           gene_associated_with_condition
                 sources: clingen, omim, diseases, uniprot, eram, disgenet,
                          ctd, monarchinitiative, ubergraph

NCBIGene:15228   Foxg1   NCBITaxon:10090   model_of          sources: mgi
NCBIGene:17257   Mecp2   NCBITaxon:10090   model_of          sources: mgi
NCBIGene:1959    EGR2    NCBITaxon:9606    genetically_associated_with
NCBIGene:2290    FOXG1   NCBITaxon:9606    genetically_associated_with
NCBIGene:2332    FMR1    NCBITaxon:9606    genetically_associated_with
```

*MECP2* ranks first because nine independent sources and three different
predicates support it — and it is in fact the causal Rett gene. Ranking is by
how many seed terms and distinct predicates back a gene, not by an opaque score.

**Two predicates are deliberately excluded.** `biolink:related_to` and
`biolink:target_for` exist between conditions and genes in ROBOKOP, but they
carry *Hetionet co-occurrence* and *drug-target* semantics. Neither means the
gene is implicated in the disease. Including them roughly doubles the gene count
and quietly fills it with noise.

### Step 4 — Cross to mouse via orthology

Human genes dominate the association data; mouse strains need mouse genes.
Alliance orthology closes the gap.

```cypher
MATCH (g)-[:`biolink:orthologous_to`]-(o:`biolink:Gene`)
WHERE g.id IN ["NCBIGene:4204"]
  AND o.taxon IN ["NCBITaxon:10090", "NCBITaxon:9606"]
RETURN DISTINCT o.id, o.name, o.taxon,
       [x IN o.equivalent_identifiers WHERE x STARTS WITH "MGI:"] AS mgi
```

```
NCBIGene:17257   Mecp2   NCBITaxon:10090   mgi=['MGI:99918']
```

That `MGI:99918` is the join key for the entire final step. It is not a property
of the node — it is pulled out of the gene's **equivalent identifiers**, the
clique of IDs that Babel considers the same entity. The graph knows this gene as
`NCBIGene:17257`; the MMRRC catalog has never heard of that identifier.

### Step 5 — Pathways

```cypher
MATCH (g)-[e]-(p:`biolink:Pathway`)
WHERE g.id IN ["NCBIGene:4204", "NCBIGene:17257"]
  AND type(e) IN ["biolink:actively_involved_in", "biolink:affects"]
WITH p, collect(DISTINCT g.name) AS gene_symbols
RETURN p.id, p.name, gene_symbols
```

```
GO:0099565   chemical synaptic transmission, postsynaptic   ['MECP2']
GO:0007219   Notch signaling pathway                        ['Mecp2', 'MECP2']
```

### Step 6 — Join to the MMRRC catalog

This step is local — no network. Every mouse gene's MGI accessions are looked up
in the in-memory index:

```
MGI:99918 (Mecp2)  →
  MMRRC:000011-UCD   STOCK Mecp2<tm1.1Jae>/Mmucd          1/1   Mecp2<tm1.1Jae>
  MMRRC:000415-UCD   B6.Cg-Mecp2<tm1.1Jae>/Mmucd          1/1   Mecp2<tm1.1Jae>
  MMRRC:011918-UCD   B6.129S4(C)-Mecp2<tm1Jae>/Mmucd      1/1   Mecp2<tm1Jae>
  MMRRC:029888-UCD   B6.Cg-Mecp2<tm1Jae>/UtaMmucd         1/1   Mecp2<tm1Jae>
  MMRRC:071393-UCD   B6J.B6N-Mecp2<tm1.1Dhy>/Mmucd        1/1   Mecp2<tm1.1Dhy>
```

### The result

```
seed terms 1 · human genes 8 · mouse genes 8 · pathways 22 · strains 103
graph: RoboMouse KG 1.0.0
```

Top of the strain table:

| Stock | Designation | Gene | Allele |
|---|---|---|---|
| `MMRRC:000011-UCD` | STOCK Mecp2\<tm1.1Jae\>/Mmucd | Mecp2 | Mecp2\<tm1.1Jae\> |
| `MMRRC:000190-UCD` | STOCK Bdnf\<tm1Lfr\>/Mmucd | Bdnf | Bdnf\<tm1Lfr\> |
| `MMRRC:000415-UCD` | B6.Cg-Mecp2\<tm1.1Jae\>/Mmucd | Mecp2 | Mecp2\<tm1.1Jae\> |
| `MMRRC:011918-UCD` | B6.129S4(C)-Mecp2\<tm1Jae\>/Mmucd | Mecp2 | Mecp2\<tm1Jae\> |

From four words of English to four orderable *Mecp2* knockouts.

---

## Two decisions that determine whether the output is usable

### 1. Rank by how much of the stock the match accounts for, not by match count

The obvious ranking — "how many of my genes does this stock carry?" — produces
garbage, and it took a live run to see it.

The MMRRC catalog is dominated by the Missouri ENU mutagenesis archive, where a
single mutagenized line carries **69–86 incidental gene annotations**. Such a
line will match 3 of your query genes where a targeted knockout matches 1. Raw
count therefore ranks every ENU line above every purpose-built model.

Searching `epilepsy` under raw-count ranking returned this:

```
MMRRC:043098-MU   C57BL/6J-MtgxR5540Btlr/Mmmh   Celf2, Dyrk1a, Phactr1
MMRRC:045729-MU   C57BL/6J-MtgxR7652Btlr/Mmmh   Alb, Rnf13, Stambp
MMRRC:038328-MU   C57BL/6J-MtgxR0034Btlr/Mmmh   Chrna7, Rhobtb2
```

Ranking by **matched ÷ annotated genes** instead — 3/86 = 0.03 for the ENU line,
1/1 = 1.00 for a knockout — returns this:

```
MMRRC:000211-UNC   B6.129S7-Chrnb2<tm1Mdb>    1/1   Chrnb2
MMRRC:000420-UNC   B6.129S7-Chrna7<tm2Bay>    1/1   Chrna7
MMRRC:029171-UCD   B6.129S1-Chrna4<tm2Lst>    1/1   Chrna4
```

*Chrnb2*, *Chrna7*, *Chrna4* — nicotinic receptor subunits, the classic ADNFLE
epilepsy genes. Both numbers are shown in the UI so the ranking is auditable,
and ENU lines can be filtered out entirely.

A second trap sits behind the same number. A GENSAT reporter line —
`Tg(Mc3r-EGFP)BX153Gsat` — is annotated with exactly one gene, so it also scores
a perfect **1.00** and ranked 6th for `obesity`. But it drives GFP off the *Mc3r*
promoter and leaves *Mc3r* itself completely intact: it is a tool for visualising
where the gene is expressed, not a model of what happens when it fails. There are
1,371 such lines in the catalog, and every one scores the maximum.

They are now flagged and sorted below real mutants at an equal ratio. The same
reporter moved from rank 6 to **rank 97** — behind all 96 genuine mutants that
also score 1.00 — and can be filtered out entirely.

The lesson generalises: the ratio measures whether the gene is the *point* of the
stock. It says nothing about what the allele actually *does*. That needs a
separate signal.

### 2. Show the allele, not just the gene

A search for `obesity` returns stocks whose gene symbol is **`a`**. That looks
like a data error. It isn't — `a` is the official MGI symbol for *agouti*
(`MGI:87853`, "nonagouti"). Mouse genetics uses single-letter symbols for
classic coat-color loci.

But the gene symbol alone is unusable. What the researcher needs is the allele:

| Stock | Gene | Allele | |
|---|---|---|---|
| `MMRRC:000137-MU` | a | `A<y>` | nonagouti; **agouti yellow** |
| `MMRRC:000375-MU` | a | `A<vy>` | nonagouti; **agouti viable yellow** |
| `MMRRC:000154-MU` | a | `a<m>` | nonagouti; mottled |

The first two are the *lethal yellow* and *viable yellow* mice — among the
best-known genetic obesity models there are. The third is a pigmentation
allele that is **not** an obesity model, despite sitting on the same gene and
therefore matching the same query.

That distinction is invisible if you only show `a`.

**The catalog does not record which allele belongs to which gene** — they are on
separate rows. The link is recovered from MGI naming conventions, measured
across all 32,163 stocks that carry alleles:

| How the link was made | Share |
|---|---|
| Allele symbol names the gene (`Esr2<tm1Unc>`→*Esr2*; `Tg(Myh6-Pln)`→*Myh6*+*Pln*) | 72.3% |
| Inferred — the stock annotates only one gene | 27.5% |
| Not linkable; left unlinked rather than guessed | 0.2% |

Because it is inference, every allele is labelled with *how* it was linked, and
the weaker single-gene inference is marked `?` in the UI.

---

## What a result means — and what it does not

Every strain returned is a **gene-level candidate**:

> This stock is annotated with a gene that the knowledge graph associates with
> your phenotype.

It is **not**:

> This stock's allele reproduces your phenotype.

`a<m>` (mottled) above is the concrete counterexample: right gene, wrong
biology. Allele, zygosity, and genetic background must be confirmed against the
MMRRC strain detail sheet and the originating publication before ordering.

Related limits worth stating out loud:

- **Absence is not evidence of absence.** No match means this catalog snapshot
  has no stock annotated with that gene — not that no such mouse exists.
- **Gene-level similarity is not allele equivalence.** Two stocks on the same
  gene may model completely different things.
- **The graph is a build.** Results are versioned against RoboMouse KG 1.0.0;
  edges change between builds. Every response carries the graph version, the
  catalog path, its row count, and a retrieval timestamp.

---

## In one slide

1. **Text → CURIE** via NameResolution, with the species selector choosing the
   vocabulary (MP for mouse, MONDO/HP for human).
2. **Cross species** via UPheno phenotype homology — when it exists.
3. **Condition → genes** via seven curated association predicates, excluding
   co-occurrence and drug-target edges that would masquerade as evidence.
4. **Human → mouse** via Alliance orthology.
5. **Genes → pathways** via GO/Reactome membership.
6. **Mouse gene → MMRRC stock** via MGI accessions pulled from the gene's
   equivalent-identifier clique.

Ranked by the gene-match ratio, annotated with the actual allele, and labelled with the
evidence that produced each hop.
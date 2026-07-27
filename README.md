# RoboMouse Strain Finder

A FastAPI service that goes from a free-text phenotype or disease to the MMRRC
mouse strains you could actually order.

```
"epilepsy"  ──NameRes──▶  MONDO:0005027
                              │
                              ├──UPheno──▶ MP:xxxxxxx  (cross-species phenotype homology)
                              │
                         ──RoboMouse KG──▶ human genes ──orthology──▶ mouse genes
                                                │                          │
                                             pathways              MGI accession
                                                                           │
                                                              ──▶ MMRRC catalog stocks
```

## Quick start

```bash
uv run robomouse-strain-finder          # http://127.0.0.1:8000
# or
uv run uvicorn robomouse_strain_finder.main:app --reload
```

On first start the app loads the MMRRC bulk catalog (~147 MB, ~589K rows,
~69K stocks) into memory, downloading it if it is not already cached at
`~/.cache/query-mmrrc-catalog/mmrrc_catalog_data.csv`. Expect a few seconds of
startup time. Interactive API docs are at `/docs`.

## Endpoints

| Route | Purpose |
| --- | --- |
| `GET /` | Landing page |
| `GET /app` | Single-page search UI |
| `GET /api/resolve` | Free text → ranked disease/phenotype CURIE candidates |
| `GET /api/search` | Full traversal → seed terms, genes, pathways, strains |
| `GET /healthz` | Catalog and upstream configuration |

```bash
curl 'http://127.0.0.1:8000/api/search?term=epilepsy&species=human&max_strains=10'
curl 'http://127.0.0.1:8000/api/search?term=obesity&species=mouse'
curl 'http://127.0.0.1:8000/api/search?term=hearing%20loss&species=human&exclude_mutation_types=CI'
```

## The UI

`GET /app` serves a single page (`src/robomouse_strain_finder/static/`) that calls
the same public API from the browser — vanilla JS, no framework, no build step.

- **Type-ahead** against `/api/resolve`. Picking a suggestion *pins* that CURIE,
  so the search skips name resolution and traverses exactly the term you chose
  rather than the top-scoring guess. The pin clears when you edit the text.
- **KPI row** of stat tiles, then collapsible Strains / Genes / Pathways /
  Seed terms / Provenance sections. Strains open by default.
- **Client-side strain filter** and **CSV export** of the full strain table.
- **Genes-matched meter** per strain, showing matched ÷ annotated genes at a glance.
- **Reporter / driver badge** on stocks that are research tools rather than models.
- Shareable URLs: `term`, `species`, and `curie` are kept in the query string,
  so `/app?term=epilepsy&species=human` runs that search on load.
- Light and dark, following the system theme.

### The landing page

`GET /` serves `static/landing.html` — the front door, and the only page that is
not the tool itself. It ships no JavaScript: the hero is a hand-authored inline
SVG, and every colour it paints with is a token in `landing.css` with a light and
a dark value, so the illustration follows the system theme along with the rest of
the page. `landing.css` loads *after* `app.css` and reuses its tokens, radii and
`.meter` / `.chip` / table styles, so the two pages stay one product.

All rendering builds DOM nodes and sets `textContent`; the client never assigns
`innerHTML`. That is deliberate — MMRRC designations and allele symbols contain
markup (`B6.129S7-<i>Chrnb2<sup>tm1Mdb</sup></i>/J`) in the source catalog. It is
stripped server-side on load *and* never parsed as HTML in the browser, and a
test asserts `app.js` contains no `innerHTML`.

### The mouse/human selector

Disease and phenotype terms are not taxon-tagged the way genes are, so the
species option selects the **vocabulary** used for name resolution rather than a
taxon filter:

- **mouse** → `MP` (Mammalian Phenotype)
- **human** → `MONDO` / `HP`

Either way the traversal reaches both species: seed terms are expanded through
UPheno `biolink:homologous_to` (HP↔MP and MONDO→MP), and the gene set is
expanded through Alliance `biolink:orthologous_to`. Picking "human" for
*epilepsy* still lands on mouse strains — via human genes and their mouse
orthologs. Set `bridge_species=false` or `include_orthologs=false` to see the
unexpanded result.

## How results are ranked, and what they mean

**Genes** are collected only through a curated set of association predicates
(`has_phenotype`, `genetically_associated_with`, `gene_associated_with_condition`,
`causes`, `contributes_to`, `model_of`, `biomarker_for`). `related_to` and
`target_for` are deliberately excluded — in ROBOKOP they carry Hetionet
co-occurrence and drug-target semantics, neither of which means the gene is
implicated in the condition. Every returned gene reports the predicates and
primary knowledge sources that produced it.

**Strains** are ranked by *how much of the stock the match accounts for*, not by
raw match count. This matters more
than it sounds: the MMRRC catalog is dominated by the Missouri ENU
mutagenesis archive, where a single stock can carry 80+ incidental gene
annotations. Ranking by "number of query genes matched" puts those lines above
every targeted knockout. Ranking by matched ÷ annotated genes puts the targeted
models first. `annotated_gene_count` and `matched_fraction` are returned so the
ranking is auditable, and `exclude_mutation_types=CI` drops the ENU lines
entirely.

**Reporter and cre-driver lines are demoted.** A GENSAT `Tg(Mc3r-EGFP)` line
carries exactly one gene and so scores a perfect 1.00 — but it drives GFP off
the promoter and leaves the gene intact. It is a research tool, not a model of
that gene's disease. Such stocks are flagged `tool_line`, sorted below real
mutants at an equal gene-match ratio, and can be removed with
`exclude_tool_lines=true`. Detection is narrow: only `Tg(...)` transgenes whose
payload is a marker or recombinase. A targeted allele that happens to insert
lacZ is still a knockout.

### Evidence caveat

A returned strain is a **gene-level candidate**. It is annotated with a gene
that the knowledge graph associates with the queried phenotype — it is *not* a
claim that the stock's allele reproduces that phenotype. Before using a stock,
confirm the allele, zygosity, and genetic background against the MMRRC strain
detail sheet and the originating publication. Gene-level and phenotype-level
similarity are not allele equivalence.

## Configuration

Environment variables, prefixed `RMSF_` (or a `.env` file):

| Variable | Default |
| --- | --- |
| `RMSF_NAME_RESOLVER_URL` | `https://name-resolution-exp.apps.renci.org` |
| `RMSF_AUTOMAT_URL` | `https://robokop-automat.apps.renci.org` |
| `RMSF_GRAPH` | `robomousekg` |
| `RMSF_MMRRC_CATALOG_PATH` | `~/.cache/query-mmrrc-catalog/mmrrc_catalog_data.csv` |
| `RMSF_AUTO_DOWNLOAD_CATALOG` | `true` |
| `RMSF_HTTP_TIMEOUT` | `120` |
| `RMSF_USE_SYSTEM_TRUST_STORE` | `true` |

Set `RMSF_GRAPH=robokopkg` to run the same traversal against baseline ROBOKOP —
though without RoboMouse's UPheno and Alliance orthology edges the mouse side
will be much thinner.

### TLS note

On networks that terminate TLS with their own CA, Python's bundled certificate
bundle rejects these hosts while `curl` succeeds. The app calls
`truststore.inject_into_ssl()` at startup so it validates against the OS trust
store instead. Disable with `RMSF_USE_SYSTEM_TRUST_STORE=false`.

## Notes on the graph

`robomousekg` has no parameter binding on its `/cypher` endpoint — it accepts a
raw query string only. Every CURIE interpolated into a query is therefore
validated against a strict pattern (`clients/automat.py`) that excludes quotes,
backslashes, and whitespace, then emitted as a JSON string literal. Anything
injection-shaped is rejected before it reaches the database.

Useful references:
- <https://robokop-automat.apps.renci.org/robomousekg/graph-metadata>
- <https://robokop-automat.apps.renci.org/robomousekg/schema>

## Further reading

[HOW_IT_WORKS.md](HOW_IT_WORKS.md) — a step-by-step walkthrough of the whole
pipeline with the actual live queries and their real responses, written for
presenting. Includes the worked *Rett syndrome* example end to end and the two
design decisions that determine whether the output is usable.

## Tests

```bash
uv run pytest
```

The suite stubs the graph and name resolver, so it needs no network and no
catalog download.
// RoboMouse Strain Finder -- browser client for /api/resolve and /api/search.
//
// No framework and no build step. All rendering goes through el()/text nodes
// rather than innerHTML: MMRRC designations and allele symbols contain markup
// in the source catalog, so data must never be parsed as HTML.

const $ = (id) => document.getElementById(id);

const form = $("search-form");
const termInput = $("term");
const suggestionList = $("suggestions");
const statusBox = $("status");
const resultsBox = $("results");
const pinnedWrap = $("pinned-wrap");
const submitButton = $("submit");

/** CURIE pinned from the suggestion list; when set, the server skips name resolution. */
let pinnedCurie = null;
let suggestions = [];
let activeSuggestion = -1;
let searchToken = 0;

// ---------------------------------------------------------------- helpers

function el(tag, props = {}, children = []) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(props)) {
    if (value === null || value === undefined || value === false) continue;
    if (key === "class") node.className = value;
    else if (key === "text") node.textContent = value;
    else if (key === "html") throw new Error("innerHTML is not allowed here");
    else if (key.startsWith("on")) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? "" : String(value));
  }
  for (const child of [].concat(children)) {
    if (child === null || child === undefined || child === false) continue;
    node.append(child.nodeType ? child : document.createTextNode(String(child)));
  }
  return node;
}

const chip = (value, cls = "chip") => el("span", { class: cls, text: value });

function chips(values, limit = Infinity) {
  const list = (values || []).filter(Boolean);
  const shown = list.slice(0, limit).map((v) => chip(v));
  if (list.length > limit) shown.push(chip(`+${list.length - limit}`));
  return shown.length ? shown : [el("span", { class: "muted", text: "—" })];
}

/** Stat-tile values are compacted; they use proportional figures, not tabular. */
function compact(n) {
  if (n === null || n === undefined) return "—";
  if (n < 1000) return String(n);
  if (n < 1e6) return `${(n / 1000).toFixed(n < 10000 ? 1 : 0)}K`;
  return `${(n / 1e6).toFixed(1)}M`;
}

const debounce = (fn, ms) => {
  let timer;
  return (...args) => {
    clearTimeout(timer);
    timer = setTimeout(() => fn(...args), ms);
  };
};

function setStatus(nodes) {
  statusBox.replaceChildren(...[].concat(nodes).filter(Boolean));
}

function speciesValue() {
  return form.querySelector('input[name="species"]:checked').value;
}

// ---------------------------------------------------------------- pinned term

function renderPinned() {
  pinnedWrap.replaceChildren();
  if (!pinnedCurie) return;
  pinnedWrap.append(
    el("span", { class: "pinned" }, [
      el("span", { text: `${pinnedCurie.label} ` }),
      el("code", { text: pinnedCurie.curie }),
      el("button", {
        type: "button",
        "aria-label": "Clear pinned term",
        title: "Clear pinned term",
        onclick: () => {
          pinnedCurie = null;
          renderPinned();
        },
        text: "×",
      }),
    ]),
  );
}

// ---------------------------------------------------------------- autocomplete

function closeSuggestions() {
  suggestions = [];
  activeSuggestion = -1;
  suggestionList.hidden = true;
  suggestionList.replaceChildren();
  termInput.setAttribute("aria-expanded", "false");
}

function renderSuggestions() {
  suggestionList.replaceChildren(
    ...suggestions.map((candidate, index) =>
      el(
        "li",
        {
          role: "option",
          id: `suggestion-${index}`,
          "aria-selected": index === activeSuggestion ? "true" : "false",
          onmousedown: (event) => {
            // mousedown, not click: blur would close the list first.
            event.preventDefault();
            choose(index);
          },
        },
        [
          el("span", { text: candidate.label }),
          el("span", { class: "curie", text: candidate.curie }),
        ],
      ),
    ),
  );
  suggestionList.hidden = suggestions.length === 0;
  termInput.setAttribute("aria-expanded", String(suggestions.length > 0));
}

function choose(index) {
  const candidate = suggestions[index];
  if (!candidate) return;
  pinnedCurie = { curie: candidate.curie, label: candidate.label };
  termInput.value = candidate.label;
  renderPinned();
  closeSuggestions();
  runSearch();
}

const fetchSuggestions = debounce(async (value) => {
  if (value.trim().length < 3) return closeSuggestions();
  try {
    const params = new URLSearchParams({ term: value.trim(), species: speciesValue(), limit: "8" });
    const response = await fetch(`/api/resolve?${params}`);
    if (!response.ok) return closeSuggestions();
    const body = await response.json();
    suggestions = body.candidates || [];
    activeSuggestion = -1;
    renderSuggestions();
  } catch {
    closeSuggestions(); // Autocomplete is a convenience; failures stay silent.
  }
}, 220);

termInput.addEventListener("input", () => {
  pinnedCurie = null;
  renderPinned();
  fetchSuggestions(termInput.value);
});

termInput.addEventListener("keydown", (event) => {
  if (suggestionList.hidden) return;
  if (event.key === "ArrowDown" || event.key === "ArrowUp") {
    event.preventDefault();
    const step = event.key === "ArrowDown" ? 1 : -1;
    activeSuggestion = (activeSuggestion + step + suggestions.length) % suggestions.length;
    renderSuggestions();
    termInput.setAttribute("aria-activedescendant", `suggestion-${activeSuggestion}`);
  } else if (event.key === "Enter" && activeSuggestion >= 0) {
    event.preventDefault();
    choose(activeSuggestion);
  } else if (event.key === "Escape") {
    closeSuggestions();
  }
});

termInput.addEventListener("blur", () => setTimeout(closeSuggestions, 120));

// ---------------------------------------------------------------- search

function buildParams() {
  const params = new URLSearchParams({
    term: termInput.value.trim(),
    species: speciesValue(),
    max_genes: $("max_genes").value,
    max_strains: $("max_strains").value,
    max_pathways: $("max_pathways").value,
    include_orthologs: String($("include_orthologs").checked),
    bridge_species: String($("bridge_species").checked),
  });
  if (pinnedCurie) params.set("curie", pinnedCurie.curie);
  if ($("exclude_enu").checked) params.append("exclude_mutation_types", "CI");
  if ($("exclude_tools").checked) params.set("exclude_tool_lines", "true");
  return params;
}

async function runSearch() {
  const term = termInput.value.trim();
  if (!term) return;
  closeSuggestions();

  const token = ++searchToken;
  submitButton.disabled = true;
  resultsBox.replaceChildren();
  setStatus(el("p", { class: "muted" }, [el("span", { class: "spinner" }), `Searching for “${term}”…`]));

  try {
    const params = buildParams();
    const response = await fetch(`/api/search?${params}`);
    const body = await response.json().catch(() => null);
    if (token !== searchToken) return; // A newer search superseded this one.

    if (!response.ok) {
      const detail = body && body.detail ? JSON.stringify(body.detail) : response.statusText;
      setStatus(el("div", { class: "error", text: `Search failed (${response.status}): ${detail}` }));
      return;
    }
    syncUrl(params);
    render(body);
  } catch (error) {
    if (token !== searchToken) return;
    setStatus(el("div", { class: "error", text: `Could not reach the API: ${error.message}` }));
  } finally {
    if (token === searchToken) submitButton.disabled = false;
  }
}

form.addEventListener("submit", (event) => {
  event.preventDefault();
  runSearch();
});

form.querySelectorAll('input[name="species"]').forEach((radio) =>
  radio.addEventListener("change", () => {
    pinnedCurie = null;
    renderPinned();
    if (termInput.value.trim()) runSearch();
  }),
);

// ---------------------------------------------------------------- rendering

function render(result) {
  const notes = (result.notes || []).map((note) => el("div", { class: "note", text: note }));
  setStatus(notes);

  resultsBox.replaceChildren();
  if (result.selected) resultsBox.append(selectedLine(result.selected));
  resultsBox.append(kpiRow(result.counts || {}));
  resultsBox.append(strainPanel(result));
  resultsBox.append(genePanel(result.genes || []));
  resultsBox.append(pathwayPanel(result.pathways || []));
  resultsBox.append(seedPanel(result.seed_terms || []));
  resultsBox.append(provenancePanel(result.provenance || {}));
}

function selectedLine(selected) {
  return el("p", { class: "tagline", style: "margin:1.1rem 0 0" }, [
    "Traversed ",
    el("strong", { text: selected.label }),
    " ",
    el("code", { text: selected.curie }),
  ]);
}

function kpiRow(counts) {
  const tiles = [
    ["Mouse genes", counts.mouse_genes, `${counts.genes_with_mgi ?? 0} with MGI`],
    ["Human genes", counts.human_genes, null],
    ["Pathways", counts.pathways, null],
    ["Seed terms", counts.seed_terms, null],
    ["MMRRC strains", counts.strains, null],
  ];
  return el(
    "ul",
    { class: "kpis" },
    tiles.map(([label, value, sub]) =>
      el("li", {}, [
        el("span", { class: "value", text: compact(value ?? 0) }),
        el("span", { class: "label", text: label }),
        sub ? el("span", { class: "sub", text: sub }) : null,
      ]),
    ),
  );
}

function panel(title, count, bodyNodes, { open = false, tools = null, footer = null } = {}) {
  const summary = el("summary", {}, [
    el("span", { text: title }),
    count !== null && count !== undefined ? el("span", { class: "count-badge", text: String(count) }) : null,
  ]);
  if (tools) {
    const holder = el("div", { class: "panel-tools" }, tools);
    // Keep control clicks from toggling the <details>.
    holder.addEventListener("click", (event) => event.stopPropagation());
    summary.append(holder);
  }
  return el("details", { class: "panel", open: open || null }, [
    summary,
    el("div", { class: "panel-body" }, bodyNodes),
    footer ? el("div", { class: "panel-note", text: footer }) : null,
  ]);
}

function table(headers, rows) {
  if (!rows.length) return el("p", { class: "empty", text: "Nothing to show." });
  return el("div", { class: "scroll" }, [
    el("table", {}, [
      el("thead", {}, [el("tr", {}, headers.map((h) => el("th", { text: h })))]),
      el("tbody", {}, rows),
    ]),
  ]);
}

function meter(fraction) {
  const pct = Math.max(0, Math.min(1, fraction || 0)) * 100;
  return el("span", { class: "meter" }, [
    el("span", { class: "meter-track" }, [
      el("span", { class: "meter-fill", style: `width:${pct.toFixed(0)}%` }),
    ]),
  ]);
}

function strainPanel(result) {
  const strains = result.strains || [];
  const filterInput = el("input", {
    type: "text",
    placeholder: "Filter strains…",
    "aria-label": "Filter strains",
  });
  const exportButton = el("button", {
    type: "button",
    class: "ghost",
    text: "Export CSV",
    onclick: () => exportCsv(result),
  });
  const body = el("div", {});

  const draw = () => {
    const needle = filterInput.value.trim().toLowerCase();
    const visible = needle
      ? strains.filter((s) =>
          [s.stock_id, s.designation, ...(s.matched_gene_symbols || []), ...(s.rrids || [])]
            .filter(Boolean)
            .join(" ")
            .toLowerCase()
            .includes(needle),
        )
      : strains;

    body.replaceChildren(
      visible.length
        ? table(["Stock", "Designation", "Matched gene / allele", "Gene → phenotype evidence", "Genes matched", "Mutation", "State", "Phenotypes"],
            visible.map((s) => strainRow(s)),
          )
        : el("p", { class: "empty", text: "No strains match this filter." }),
    );
  };

  filterInput.addEventListener("input", draw);
  draw();

  return panel("Strains", strains.length, body, {
    open: true,
    tools: [filterInput, exportButton],
    footer:
      "Matched by MGI gene accession. A gene-level link — not a claim that the stock's " +
      "allele reproduces the queried phenotype. The allele shown is inferred from MGI " +
      "naming, since the catalog does not state which allele belongs to which gene; “?” " +
      "marks one inferred only because the stock annotates a single gene. Confirm allele, " +
      "zygosity, and background against the strain detail sheet before ordering.",
  });
}

function strainRow(s) {
  const stockCell = el("td", {}, [
    s.sds_url
      ? el("a", { href: s.sds_url, target: "_blank", rel: "noopener", text: s.stock_id })
      : el("span", { text: s.stock_id }),
    ...(s.rrids || []).map((r) => el("div", { class: "chip", text: r })),
    s.tool_line
      ? el("div", {
          class: "chip tool",
          title: "Reporter or cre-driver transgene: the gene is not disrupted, so this is a research tool rather than a disease model.",
          text: "reporter / driver",
        })
      : null,
  ]);

  // Gene symbol alone can be cryptic -- mouse uses single-letter symbols for
  // classic loci (agouti is `a`). Pair each matched gene with the allele the
  // stock actually carries, so `a` reads as `a  A<y> (agouti yellow)`.
  const matchedAlleles = s.matched_alleles || [];
  const geneCell = el("td", {}, [
    ...chips(s.matched_gene_symbols),
    ...matchedAlleles.slice(0, 3).map((a) =>
      el("div", { class: "allele", title: a.name || "" }, [
        el("span", { class: "mono", text: a.symbol || a.mgi_allele_id || "—" }),
        a.link === "sole-gene" ? el("span", { class: "muted", text: " ?" }) : null,
        a.name ? el("div", { class: "muted allele-name", text: a.name }) : null,
      ]),
    ),
    matchedAlleles.length > 3 ? chip(`+${matchedAlleles.length - 3} alleles`) : null,
  ]);

  return el("tr", {}, [
    stockCell,
    el("td", { text: s.designation || "—" }),
    geneCell,
    evidenceCell(s.gene_evidence || []),
    el("td", { class: "num" }, [
      meter(s.matched_fraction),
      el("span", {
        class: "muted",
        text: ` ${(s.matched_mgi_gene_ids || []).length}/${s.annotated_gene_count ?? 0}`,
      }),
    ]),
    el("td", {}, chips(s.mutation_types)),
    el("td", {}, chips(s.states)),
    el("td", {}, (s.phenotypes || []).slice(0, 3).map((p) => el("div", { class: "muted", text: p }))
      .concat(
        (s.phenotypes || []).length > 3
          ? [chip(`+${s.phenotypes.length - 3} more`)]
          : (s.phenotypes || []).length
            ? []
            : [el("span", { class: "muted", text: "—" })],
      )),
  ]);
}

/** Why each matched gene is linked to the queried phenotype.
 *  On an ortholog row the evidence belongs to the human partner, not the mouse
 *  gene the stock carries -- so that is stated rather than implied. */
function evidenceCell(evidence) {
  if (!evidence.length) return el("td", {}, [el("span", { class: "muted", text: "—" })]);
  return el(
    "td",
    {},
    evidence.slice(0, 3).map((e) =>
      el("div", { class: "evidence" }, [
        ...chips((e.predicates || []).map((p) => p.replace("biolink:", "").replace(/_/g, " "))),
        e.via === "ortholog"
          ? el("div", { class: "muted evidence-note" }, [
              "via ortholog of ",
              el("span", { class: "mono", text: e.ortholog_of_symbol || e.ortholog_of || "?" }),
              " (human)",
            ])
          : null,
        (e.knowledge_sources || []).length
          ? el("div", { class: "muted evidence-note" }, [
              (e.knowledge_sources || [])
                .slice(0, 4)
                .map((k) => k.replace("infores:", ""))
                .join(", "),
              (e.knowledge_sources || []).length > 4
                ? ` +${e.knowledge_sources.length - 4}`
                : "",
            ])
          : null,
      ]),
    ),
  );
}

function genePanel(genes) {
  const rows = genes.map((g) =>
    el("tr", {}, [
      el("td", {}, [el("code", { text: g.curie })]),
      el("td", { text: g.symbol || "—" }),
      el("td", { text: g.species || g.taxon || "—" }),
      el("td", {}, [
        el("span", { text: g.via }),
        g.ortholog_of ? el("div", { class: "muted", text: `of ${g.ortholog_of}` }) : null,
      ]),
      el("td", {}, chips((g.predicates || []).map((p) => p.replace("biolink:", "")))),
      el("td", {}, chips(g.mgi_ids)),
      el("td", {}, chips((g.knowledge_sources || []).map((s) => s.replace("infores:", "")), 4)),
    ]),
  );
  return panel(
    "Genes",
    genes.length,
    table(["Gene", "Symbol", "Species", "Via", "Predicates", "MGI", "Sources"], rows),
  );
}

function pathwayPanel(pathways) {
  const rows = pathways.map((p) =>
    el("tr", {}, [
      el("td", {}, [el("code", { text: p.curie })]),
      el("td", { text: p.name || "—" }),
      el("td", { class: "num", text: String((p.gene_curies || []).length) }),
      el("td", {}, chips(p.gene_symbols, 10)),
    ]),
  );
  return panel("Pathways", pathways.length, table(["Pathway", "Name", "Genes", "Members"], rows));
}

function seedPanel(seeds) {
  const rows = seeds.map((s) =>
    el("tr", {}, [
      el("td", {}, [el("code", { text: s.curie })]),
      el("td", { text: s.label || "—" }),
      el("td", { text: s.via }),
      el("td", {}, s.homologous_to ? [el("code", { text: s.homologous_to })] : ["—"]),
    ]),
  );
  return panel("Seed terms", seeds.length, table(["CURIE", "Label", "Via", "Homologous to"], rows));
}

function provenancePanel(p) {
  const line = (label, value) =>
    value ? el("div", {}, [el("span", { class: "muted", text: `${label}: ` }), el("code", { text: value })]) : null;
  return panel("Provenance", null, [
    el("div", { class: "prov" }, [
      line("Name Resolver", p.name_resolver_url),
      line("Graph", p.graph_url),
      line("Graph version", p.graph_version),
      line("Catalog", p.catalog_path),
      line("Catalog rows", p.catalog_rows ? String(p.catalog_rows) : null),
      line("Catalog modified", p.catalog_modified_at),
      line("Retrieved", p.retrieved_at),
    ]),
  ]);
}

// ---------------------------------------------------------------- CSV export

function exportCsv(result) {
  const headers = [
    "stock_id", "rrids", "designation", "matched_gene_symbols", "matched_mgi_gene_ids",
    "matched_fraction", "annotated_gene_count", "matched_allele_symbols",
    "matched_allele_names", "matched_allele_link", "evidence_predicates",
    "evidence_sources", "evidence_via", "tool_line", "allele_symbols", "mgi_allele_ids",
    "mutation_types", "states", "strain_types", "phenotypes", "pubmed_ids",
    "research_areas", "accepted_date", "sds_url",
  ];
  const cell = (value) => {
    const text = Array.isArray(value) ? value.join("; ") : value === null || value === undefined ? "" : String(value);
    return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
  };
  const lines = [headers.join(",")];
  for (const s of result.strains || []) {
    lines.push([
      s.stock_id, s.rrids, s.designation, s.matched_gene_symbols, s.matched_mgi_gene_ids,
      (s.matched_fraction ?? 0).toFixed(3), s.annotated_gene_count,
      (s.matched_alleles || []).map((a) => a.symbol).filter(Boolean),
      (s.matched_alleles || []).map((a) => a.name).filter(Boolean),
      [...new Set((s.matched_alleles || []).map((a) => a.link))],
      [...new Set((s.gene_evidence || []).flatMap((e) => e.predicates || []))]
        .map((p) => p.replace("biolink:", "")),
      [...new Set((s.gene_evidence || []).flatMap((e) => e.knowledge_sources || []))]
        .map((k) => k.replace("infores:", "")),
      [...new Set((s.gene_evidence || []).map((e) => e.via))],
      String(s.tool_line ?? false),
      (s.alleles || []).map((a) => a.symbol).filter(Boolean),
      (s.alleles || []).map((a) => a.mgi_allele_id).filter(Boolean),
      s.mutation_types, s.states, s.strain_types, s.phenotypes, s.pubmed_ids,
      s.research_areas, s.accepted_date, s.sds_url,
    ].map(cell).join(","));
  }

  const slug = (result.selected?.curie || result.query || "strains").replace(/[^A-Za-z0-9]+/g, "-");
  const url = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8" }));
  const link = el("a", { href: url, download: `mmrrc-${slug}.csv` });
  document.body.append(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

// ---------------------------------------------------------------- URL state

function syncUrl(params) {
  const keep = new URLSearchParams();
  for (const key of ["term", "species", "curie"]) {
    if (params.get(key)) keep.set(key, params.get(key));
  }
  history.replaceState(null, "", keep.toString() ? `/app?${keep}` : "/app");
}

function restoreFromUrl() {
  const params = new URLSearchParams(location.search);
  const term = params.get("term");
  if (!term) return;
  termInput.value = term;
  const species = params.get("species");
  if (species) {
    const radio = form.querySelector(`input[name="species"][value="${CSS.escape(species)}"]`);
    if (radio) radio.checked = true;
  }
  const curie = params.get("curie");
  if (curie) {
    pinnedCurie = { curie, label: term };
    renderPinned();
  }
  runSearch();
}

restoreFromUrl();
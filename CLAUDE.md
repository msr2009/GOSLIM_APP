# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A data analysis project mapping *C. elegans* disease-associated genes to their subcellular
localization (GO Cellular Component slim terms), using human-worm orthology and disease
annotation data from the Alliance of Genome Resources and WormBase. It has two parts: a batch
pipeline script (`map_genes_to_go_cc.py`) and an interactive Shiny app (`goslim_app/`) for
building/curating the GO slim definition the pipeline uses.

Not a git repository — no VCS history to reference.

## Environment

- Python env: `go_terms` (conda/mamba), defined in `environment.yml`
  (deps: `pandas`, `goatools`, `matplotlib`, `shiny`)
- Run with `mamba run -n go_terms python <script>.py ...` or activate the env first

## Running the pipeline

```
python map_genes_to_go_cc.py \
    --gaf gene_association.WS298.wb \
    --obo gene_ontology.WS298.obo \
    --slim goslim_generic.obo goslim_agr.obo \
    --genes data/worm_disease_genes.tsv \
    --out-prefix worm_disease_go
```

Writes `<out-prefix>_per_term.tsv` (one row per gene × slim term) and `<out-prefix>_per_gene.tsv`
(one row per gene, multi-valued fields joined with `;`).

`map_genes_to_go_cc.py` is the productionized version of the exploratory logic in
`tag_orthology_disease_analysis.ipynb` — prefer editing the script over the notebook for
anything beyond one-off exploration.

## Running the GO slim explorer app

```
shiny run --reload goslim_app/app.py
```

An interactive curation tool for the GO slim definition consumed by the pipeline: lets a
curator promote uncovered terms into the slim, remove terms, drill a term down into its
more-specific subterms, flag/exclude uninformative or obsolete terms, and inspect
coverage/group-size stats — all recomputed live against the same GAF/OBO/slim files the
pipeline uses. Finished slim definitions save as JSON under `goslim_app/saved_slims/`.
Per-gene gene→group tables and gene lists can be exported as TSV/text from the UI.

Deploying a hosted copy to shinyapps.io: see `goslim_app/DEPLOY.md` — run
`goslim_app/prepare_deploy.sh` first to stage repo-root data files into
`goslim_app/vendor/`, since a shinyapps deploy only bundles the app directory.

App architecture (`goslim_app/`):
- `config.py` — all data-file paths, resolved relative to the repo root (`APP_DIR.parent`), or
  to `goslim_app/vendor/` instead when that directory exists (see `DEPLOY.md`)
- `startup.py` — builds the in-memory `DATA` bundle once at import (GAF, GO DAG, ancestor
  closure, base slims, gene lists) — this is the slow step; everything downstream recomputes
  in milliseconds against it
- `go_data.py` — raw file parsers (GAF, obsolete-term map, slim ids, gene lists)
- `ancestors.py` — ancestor-closure builder and memoized slim mapping (a faster replacement
  for goatools' `mapslim`, used by the app instead of the pipeline script's approach)
- `annotations.py` — reduces the GAF to `{aspect: {gene: frozenset(go_ids)}}`
- `aspects.py` — per-aspect qualifier/root/namespace tables and `EXTRA_SLIM_TERMS`
- `slimming.py` — per-request recompute of gene→slim-group coverage
- `partition.py` — forces a single gene→one-group assignment (rarest-group-wins tie-break)
- `drilldown.py` / `excluded.py` / `saved_slims.py` / `exports.py` — term drilldown, the
  exclude-list audit/persistence, saved-slim JSON I/O, and TSV/text export builders
- `ui_layout.py` / `server_main.py` — UI layout (5 tabs: Coverage, Slim groups, Add terms,
  Gene detail, Slim hygiene) and all reactive wiring

Gotchas specific to the app:
- `excluded_terms.json` at the repo root is shared, global state the app mutates directly
  on every exclude/un-exclude — it is not versioned per-slim.
- The "use part_of" checkbox in the UI is saved into exported slim JSON but is **not** wired
  to the live ancestor closure (only `config.USE_PART_OF_DEFAULT`, read once at startup,
  controls that) — toggling it in the UI has no runtime effect.
- `EXTRA_SLIM_TERMS` only applies for the Cellular Component aspect.

## Building the proteome gene list

`build_gene_list.py` downloads the WormBase geneIDs file (via the EBI FTP mirror, since
downloads.wormbase.org sits behind a Cloudflare challenge) and derives
`data/celegans_proteome.tsv` — one protein-coding gene symbol per line, used by the app as an
alternative gene list to `data/worm_disease_genes.tsv`.

```
mamba run -n go_terms python build_gene_list.py --release WS298
```

## Data flow / architecture

1. **Orthology**: `ORTHOLOGY-ALLIANCE_COMBINED.tsv` (Alliance orthology calls across species)
   is filtered to `Gene1SpeciesName == Caenorhabditis elegans` /
   `Gene2SpeciesName == Homo sapiens` to get worm genes with human orthologs.
2. **Disease genes**: `DISEASE-ALLIANCE_WB.tsv` filtered to `DBobjectType == 'gene'` gives the
   set of worm genes with disease associations (`data/worm_disease_genes.tsv` is this gene list,
   one symbol per line — the `--genes` input to the script).
3. **GO annotation**: `gene_association.WS298.wb` (GAF format, see `GAF_COLS` in the script for
   column layout) is filtered to CC-aspect rows with `Qualifier` in `{located_in, part_of}` to
   get each gene's cellular-component GO terms. Other qualifiers (e.g. `enables`, `involved_in`)
   are intentionally excluded as not indicating physical localization.
4. **GO DAG + slimming**: `gene_ontology.WS298.obo` is the full ontology, release-matched to the
   GAF (do not mix with `go-basic.obo`, which is a different/generic release kept around from
   earlier exploration). `goslim_generic.obo` and `goslim_agr.obo` are merged into one custom
   slim DAG via `build_custom_slim_dag`, then patched with `EXTRA_SLIM_TERMS` — CC terms observed
   to be common in the worm GAF but missing from both stock slims.
5. **Obsolete term handling**: `build_replaced_by_map` parses `replaced_by:` stanzas directly out
   of the OBO file (goatools' `GODag` doesn't surface this) so obsolete GO IDs in the GAF get
   substituted with their replacement before slim mapping.
6. **Slim mapping**: `goatools.mapslim.mapslim` maps each specific GO term to its direct
   ancestor(s) in the custom slim DAG. `GO_CC_ROOT` (`GO:0005575`, the bare "cellular_component"
   root) is treated as uninformative and folded into `None` unless `--keep-root` is passed.
7. **Output**: `per_term`/`per_gene` TSVs (see "Running the pipeline" above) — one row per unique
   `(gene, go_slim_term)` pair, or one row per gene with multi-valued fields joined by `;`.

## Gotchas

- GAF and OBO releases must match (`gene_association.WS298.wb` ↔ `gene_ontology.WS298.obo`) —
  mismatched releases cause obsolete-ID misses even with the replaced_by patch.
- `worm_disease_go_per_term.tsv` / `worm_disease_go_per_gene.tsv` in the repo root are generated
  outputs, not inputs — regenerate by rerunning `map_genes_to_go_cc.py` rather than hand-editing.
  `worm_disease_go.tsv` is an older single-file output from before the per-term/per-gene split.
- `c_elegans.PRJNA13758.WS279.go_annotations.gaf` is a separate, older-release GAF from earlier
  exploration — the pipeline and app both use `gene_association.WS298.wb`, not this file.
- The notebook (`tag_orthology_disease_analysis.ipynb`) has known rough edges the script fixed:
  it uses `go-basic.obo` for the replaced_by map while annotating against
  `gene_ontology.WS298.obo` (release mismatch), and it hardcodes a fix for one obsolete term
  (`GO:0005615`) rather than using the general replaced_by mechanism.
- `symbol_not_found_*.txt` files are exports from the app's coverage tab (gene symbols present
  in a gene list but absent from the GAF), not hand-maintained inputs.

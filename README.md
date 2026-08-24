# GO Slim Exporer - a Shiny app for exploring grouping *C. elegans* genes by GO terms

Proteome-scale gene tagging in animals will require multi-lab collaboration, so
how do you split up the genes? This app allows users to interactively group
genes by Gene Ontology terms, starting with GO slim terms lists provided by the
Alliance of Genome Resources and Wormbase. 

Two parts:

- **`map_genes_to_go_cc.py`** — batch pipeline script
- **`goslim_app/`** — interactive Shiny app for building/curating the GO slim definition the
  pipeline uses

## Environment

Python env `go_terms` (conda/mamba), defined in `environment.yml`:

- `pandas`
- `goatools`
- `matplotlib`
- `shiny`

```
mamba env create -f environment.yml
mamba activate go_terms
```

## Running the pipeline

```
python map_genes_to_go_cc.py \
    --gaf gene_association.WS298.wb \
    --obo gene_ontology.WS298.obo \
    --slim goslim_generic.obo goslim_agr.obo \
    --genes data/worm_disease_genes.tsv \
    --out-prefix worm_disease_go
```

Writes:

- `<out-prefix>_per_term.tsv` — one row per gene x slim term
- `<out-prefix>_per_gene.tsv` — one row per gene, multi-valued fields joined with `;`

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

The default gene list on load is `celegans_proteome` (falls back to whichever gene list is
first alphabetically if that file is missing).

The Coverage tab's group-size bar chart also overlays two diamond-marker traces sourced from
[WormTagDB](https://wormtagdb.rc.duke.edu/) (`wormtagdb_tagged_genes.tsv`, at the repo root):
per-group counts of genes with an endogenous fluorescent-protein tag ("WormTagDB fluorescent
tags"), and the subset of those with a CGC-orderable strain ("CGC-available fluorescent
tags"). Each marker's tooltip shows the count and % of that group's genes. These are ordinary
plotly traces, so click their legend entries to toggle visibility — no separate control needed.

Deploying a hosted copy to shinyapps.io: see `goslim_app/DEPLOY.md` — run
`goslim_app/prepare_deploy.sh` first to stage repo-root data files into
`goslim_app/vendor/`, since a shinyapps deploy only bundles the app directory.

### App architecture (`goslim_app/`)

| File | Role |
|---|---|
| `config.py` | Data-file paths, resolved relative to the repo root, or to `goslim_app/vendor/` when that directory exists (see `DEPLOY.md`) |
| `startup.py` | Builds the in-memory `DATA` bundle once at import (GAF, GO DAG, ancestor closure, base slims, gene lists, WormTagDB tagged-gene sets) — the slow step; everything downstream recomputes in milliseconds against it |
| `go_data.py` | Raw file parsers (GAF, obsolete-term map, slim ids, gene lists, WormTagDB tagged-gene TSV) |
| `ancestors.py` | Ancestor-closure builder and memoized slim mapping (a faster replacement for goatools' `mapslim`) |
| `annotations.py` | Reduces the GAF to `{aspect: {gene: frozenset(go_ids)}}` |
| `aspects.py` | Per-aspect qualifier/root/namespace tables and `EXTRA_SLIM_TERMS` |
| `slimming.py` | Per-request recompute of gene→slim-group coverage |
| `partition.py` | Forces a single gene→one-group assignment (rarest-group-wins tie-break) |
| `drilldown.py` / `excluded.py` / `saved_slims.py` / `exports.py` | Term drilldown, exclude-list audit/persistence, saved-slim JSON I/O, TSV/text export builders |
| `ui_layout.py` / `server_main.py` | UI layout (5 tabs: Coverage, Slim groups, Add terms, Gene detail, Slim hygiene) and reactive wiring |

Gotchas specific to the app:

- `excluded_terms.json` at the repo root is shared, global state the app mutates directly on
  every exclude/un-exclude — it is not versioned per-slim.
- The "use part_of" checkbox in the UI is saved into exported slim JSON but is **not** wired to
  the live ancestor closure (only `config.USE_PART_OF_DEFAULT`, read once at startup, controls
  that) — toggling it in the UI has no runtime effect.
- `EXTRA_SLIM_TERMS` only applies for the Cellular Component aspect.
- `wormtagdb_tagged_genes.tsv` at the repo root is a snapshot pulled from WormTagDB's search
  export (not auto-updating) — one row per tagged gene, with `fluor_tags`/`other_tags`/
  `cgc_strains` columns (semicolon-joined; empty when not applicable). It's kept out of
  `data/` deliberately: `GENE_LIST_DIR` auto-loads every `*.tsv`/`*.txt` there as a plain
  gene list, and this file's extra columns would corrupt that. `prepare_deploy.sh` vendors it
  alongside the other repo-root data files for shinyapps.io deploys.

## Data sources (not committed to this repo — download before running)

These are large, externally versioned files (`.gitignore`'d here) — fetch them into the repo
root before running the pipeline or app. Get GAF/OBO releases that match (see Gotchas below).

```bash
# WormBase GAF + matching GO ontology, via the EBI mirror (downloads.wormbase.org sits
# behind a Cloudflare challenge that blocks plain curl/requests)
curl -O https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS298/ONTOLOGY/gene_association.WS298.wb
curl -O https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS298/ONTOLOGY/gene_ontology.WS298.obo

# generic/current GO ontology (a different, non-release-pinned tree - see Gotchas)
curl -o go-basic.obo http://purl.obolibrary.org/obo/go/go-basic.obo

# Alliance of Genome Resources: combined cross-species orthology calls
curl -o ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz \
    https://download.alliancegenome.org/9.0.0/downloads/ORTHOLOGY-ALLIANCE_TSV_COMBINED.tsv.gz
gunzip ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz
```

`DISEASE-ALLIANCE_WB.tsv` (Alliance disease associations, C. elegans only) isn't reliably
available as a direct TSV URL in the current Alliance release schema — download it from the
[Alliance downloads page](https://www.alliancegenome.org/downloads) (Disease Associations,
WB) instead; if only JSON is offered, convert the WB JSON export to TSV.

The Alliance release version (`9.0.0` above) changes over time — check
https://www.alliancegenome.org/api/downloads for the current release's file list/URLs if the
link above 404s.

## Building the proteome gene list

`build_gene_list.py` downloads the WormBase geneIDs file (via the EBI FTP mirror, since
downloads.wormbase.org sits behind a Cloudflare challenge) and derives
`data/celegans_proteome.tsv` — one protein-coding gene symbol per line, used by the app as an
alternative gene list to `data/worm_disease_genes.tsv`.

```
mamba run -n go_terms python build_gene_list.py --release WS298
```

## Data flow / architecture

1. **Orthology** — `ORTHOLOGY-ALLIANCE_COMBINED.tsv` (Alliance orthology calls across species)
   is filtered to `Gene1SpeciesName == Caenorhabditis elegans` /
   `Gene2SpeciesName == Homo sapiens` to get worm genes with human orthologs.
2. **Disease genes** — `DISEASE-ALLIANCE_WB.tsv` filtered to `DBobjectType == 'gene'` gives the
   set of worm genes with disease associations (`data/worm_disease_genes.tsv` is this gene
   list, one symbol per line — the `--genes` input to the script).
3. **GO annotation** — `gene_association.WS298.wb` (GAF format) is filtered to CC-aspect rows
   with `Qualifier` in `{located_in, part_of}` to get each gene's cellular-component GO terms.
   Other qualifiers (e.g. `enables`, `involved_in`) are intentionally excluded as not indicating
   physical localization.
4. **GO DAG + slimming** — `gene_ontology.WS298.obo` is the full ontology, release-matched to
   the GAF (do not mix with `go-basic.obo`, a different/generic release kept around from earlier
   exploration). `goslim_generic.obo` and `goslim_agr.obo` are merged into one custom slim DAG,
   then patched with `EXTRA_SLIM_TERMS` — CC terms observed to be common in the worm GAF but
   missing from both stock slims.
5. **Obsolete term handling** — `replaced_by:` stanzas are parsed directly out of the OBO file
   (goatools' `GODag` doesn't surface this) so obsolete GO IDs in the GAF get substituted with
   their replacement before slim mapping.
6. **Slim mapping** — `goatools.mapslim.mapslim` maps each specific GO term to its direct
   ancestor(s) in the custom slim DAG. `GO:0005575` (the bare "cellular_component" root) is
   treated as uninformative and folded into `None` unless `--keep-root` is passed.
7. **Output** — `per_term`/`per_gene` TSVs, one row per unique `(gene, go_slim_term)` pair, or
   one row per gene with multi-valued fields joined by `;`.

## Gotchas

- GAF and OBO releases must match (`gene_association.WS298.wb` <-> `gene_ontology.WS298.obo`) —
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

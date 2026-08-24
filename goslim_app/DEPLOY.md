# Deploying to shinyapps.io

The app reads several large data files (`gene_association.WS298.wb`,
`gene_ontology.WS298.obo`, the two `goslim_*.obo` files, `data/*.tsv`,
`excluded_terms.json`) from the repo root, one directory above `goslim_app/`.
shinyapps.io only uploads the app directory you point `rsconnect` at, so those
files need to be staged *inside* `goslim_app/` first.

## One-time setup

```
pip install rsconnect-python
rsconnect add --account <your-shinyapps-account> --name <nickname> --token <token> --secret <secret>
```

(Token/secret come from shinyapps.io -> Account -> Tokens.)

## Every deploy

```
bash goslim_app/prepare_deploy.sh      # stages gene_association.WS298.wb, the
                                        # OBOs, data/*.tsv, excluded_terms.json
                                        # into goslim_app/vendor/
rsconnect deploy shiny goslim_app \
    --name <nickname> \
    --title goslim-explorer
```

`config.py` prefers `goslim_app/vendor/` over the repo-root copies whenever
that directory exists, so nothing else needs to change for local development
- `vendor/` just doesn't exist there.

`requirements.txt` in this directory pins the Python deps (`pandas`,
`goatools`, `plotly`, `shiny`, `shinywidgets`, `openpyxl`) that
`rsconnect deploy` installs into the app's environment; keep it in sync with
`environment.yml` if you upgrade a dependency.

## Known limitations of a shinyapps.io deployment

- **`excluded_terms.json` and `saved_slims/*.json` are mutated at runtime.**
  shinyapps.io's filesystem is writable while an instance is running but is
  **not** persisted across redeploys or instance restarts, and multiple
  concurrent worker processes each get their own copy. Curation changes made
  in a deployed instance can silently disappear or fail to show up for other
  users. This is fine for a single-curator demo/review session; it is not a
  substitute for the local workflow if you need durable, shared curation
  state - export slims (the JSON/xlsx downloads) instead of relying on the
  server retaining them.
- **Bundle size.** The vendored data adds ~80 MB to the deploy bundle. Well
  within shinyapps.io's limits, but re-run `prepare_deploy.sh` (don't hand-edit
  `vendor/`) whenever the source GAF/OBO files change so the bundle doesn't go
  stale.
- **Startup time.** `startup.py` builds the full in-memory `DATA` bundle
  (GO DAG, ancestor closure, GAF index) once per process at import; expect a
  ~1-2s cold start per instance, same as running locally.

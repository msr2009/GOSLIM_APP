#!/usr/bin/env bash
# Deploy the GO slim explorer (goslim_app/) to shinyapps.io, or dry-run the bundle.
#
#   ./deploy.sh              stage vendor data, then deploy (updates the existing app in place)
#   ./deploy.sh --manifest   stage vendor data, print what WOULD be uploaded, then clean up
#   ./deploy.sh --name foo   any extra args are passed through to `rsconnect deploy shiny`
#
# WHY THERE IS NO .rscignore
# --------------------------
# .rscignore is an R rsconnect feature. rsconnect-python does not read it: in 1.30.0
# (tested 2026-09-30) the string "rscignore" does not appear anywhere in the package, and a
# test bundle with a .rscignore naming a file and a directory included both of them, plus
# the .rscignore file itself. Do not "fix" bundle bloat by editing a .rscignore -- add a
# pattern to EXCLUDES below and check the result with --manifest.
#
# rsconnect-python's only built-in ignores (bundle.py, directories_ignore_list) are
# .git/ .svn/ __pycache__/ node_modules/ renv/ packrat/ rsconnect/ rsconnect-python/
# .Rproj.user/ -- and __pycache__/ only at the TOP level of the deploy dir. Everything
# else ships.
#
# Exclude-glob notes (patterns are relative to goslim_app/):
#   - '**/*.pyc' is the only form that catches nested __pycache__; '**/__pycache__/' and
#     '*/__pycache__/' silently match nothing.
#   - '*~' catches both vim backups (foo~) and vim undo files (.foo.un~).
#   - a bare directory name ('somedir/') is expanded by rsconnect to somedir/**/*.
#
# WHAT MUST SHIP (the app hard-fails at startup without these; see goslim_app/config.py):
#   app.py and every *.py module, requirements.txt (the dependency manifest on shinyapps.io),
#   vendor/gene_association.WS298.wb, vendor/gene_ontology.WS298.obo,
#   vendor/goslim_{agr,generic}.obo, vendor/data/*.tsv gene lists,
#   vendor/excluded_terms.json, vendor/wormtagdb_tagged_genes.tsv, saved_slims/*.json.
# Do not exclude *.py files based on an import scan -- keep all code.
#
# The existing app's id lives in goslim_app/rsconnect-python/goslim_app.json, so a plain
# `rsconnect deploy shiny goslim_app` updates that app in place. Never pass --app-id and
# never deploy the repo root (it holds ~400 MB of source data the app does not read).
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="$REPO_DIR/goslim_app"
ENTRYPOINT="app:app"

EXCLUDES=(
    # compiled bytecode in nested __pycache__ dirs (top-level one is ignored by rsconnect)
    -x '**/*.pyc'
    # editor backups and vim undo files
    -x '*~'
    -x '**/.DS_Store'
    # notebooks, tool caches, dev-only env spec (requirements.txt is the deploy manifest)
    -x '**/*.ipynb'
    -x '.ipynb_checkpoints/'
    -x '.ruff_cache/'
    -x '.claude/'
    -x 'environment.yml'
    # docs and deploy tooling, not read at runtime
    -x '**/MODELS.md'
    -x 'DEPLOY.md'
    -x 'prepare_deploy.sh'
    # multi-column derived tables that prepare_deploy.sh copies from data/; config.py would
    # auto-load them as (garbage) gene lists since it globs every *.tsv in vendor/data/
    -x 'vendor/data/human_disease_worm_orthology_*'
)

# refresh vendor/ from the repo-root data so the bundle is never stale
bash "$APP_DIR/prepare_deploy.sh" >/dev/null

# dry-run mode: build a manifest, summarize it, and remove it
if [[ "${1:-}" == "--manifest" ]]; then
    MANIFEST="$APP_DIR/manifest.json"
    # a stray manifest.json changes how later deploys resolve dependencies - always remove it
    trap 'rm -f "$MANIFEST"' EXIT
    rsconnect write-manifest shiny "$APP_DIR" --entrypoint "$ENTRYPOINT" --overwrite \
        "${EXCLUDES[@]}" >/dev/null
    python3 - "$APP_DIR" "$MANIFEST" <<'EOF'
import json, os, sys

app_dir, manifest = sys.argv[1], sys.argv[2]
files = sorted(json.load(open(manifest))["files"])
sizes = {f: os.path.getsize(os.path.join(app_dir, f)) for f in files}
print(f"{len(files)} files, {sum(sizes.values()) / 1048576:.1f} MB")
print("files over 500 KB:")
for f in files:
    if sizes[f] > 500 * 1024:
        print(f"  {sizes[f] / 1048576:8.1f} MB  {f}")
print("all files:")
for f in files:
    print(f"  {f}")
EOF
    exit 0
fi

rsconnect deploy shiny "$APP_DIR" --entrypoint "$ENTRYPOINT" "${EXCLUDES[@]}" "$@"

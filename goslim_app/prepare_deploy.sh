#!/usr/bin/env bash
# Stage a self-contained copy of the data files this app needs under
# goslim_app/vendor/, so a shinyapps.io deploy (which only bundles the app
# directory) has everything config.py looks for. Run this before every
# `rsconnect deploy shiny goslim_app`; re-run whenever the source data
# files change, since it copies rather than links.
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$APP_DIR")"
VENDOR_DIR="$APP_DIR/vendor"

mkdir -p "$VENDOR_DIR/data"

cp "$REPO_DIR/gene_association.WS298.wb" "$VENDOR_DIR/"
cp "$REPO_DIR/gene_ontology.WS298.obo" "$VENDOR_DIR/"
cp "$REPO_DIR/goslim_agr.obo" "$VENDOR_DIR/"
cp "$REPO_DIR/goslim_generic.obo" "$VENDOR_DIR/"
cp "$REPO_DIR/data/"*.tsv "$VENDOR_DIR/data/"

if [ -f "$REPO_DIR/wormtagdb_tagged_genes.tsv" ]; then
    cp "$REPO_DIR/wormtagdb_tagged_genes.tsv" "$VENDOR_DIR/"
fi

if [ -f "$REPO_DIR/excluded_terms.json" ]; then
    cp "$REPO_DIR/excluded_terms.json" "$VENDOR_DIR/"
fi

echo "Staged vendor data under $VENDOR_DIR:"
du -sh "$VENDOR_DIR"/* "$VENDOR_DIR"/data/* 2>/dev/null

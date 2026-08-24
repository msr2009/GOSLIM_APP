"""Persistent exclude list plus an audit of base slims for junk terms
(obsolete, replaced, missing, uninformative, or empty for the current genes)."""

import json
from pathlib import Path

import pandas as pd

# seeded on first run; each is verified against the release-matched OBO
# rather than trusted blindly, since GO ids can be split/merged over time
SEED_EXCLUDES = {
    "global": [
        {"id": "GO:0005575", "name": "cellular_component", "reason": "aspect root, uninformative"},
        {"id": "GO:0003674", "name": "molecular_function", "reason": "aspect root, uninformative"},
        {"id": "GO:0008150", "name": "biological_process", "reason": "aspect root, uninformative"},
    ],
    "C": [
        {"id": "GO:0110165", "name": "cellular anatomical entity", "reason": "near-root, groups everything"},
    ],
    "F": [],
    "P": [],
}

UNINFORMATIVE_DEPTH = 1
UNINFORMATIVE_COVERAGE_FRACTION = 0.6


def load_excluded(path):
    """Read excluded_terms.json, seeding it with SEED_EXCLUDES if absent."""
    path = Path(path)
    if not path.exists():
        save_excluded(path, SEED_EXCLUDES)
        return dict(SEED_EXCLUDES)
    with open(path) as f:
        return json.load(f)


def save_excluded(path, data):
    """Write the exclude list back out, human-readable and diff-friendly."""
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def excluded_ids_for_aspect(excluded, aspect):
    """Flatten global + per-aspect excludes into a single id set."""
    ids = {e["id"] for e in excluded.get("global", [])}
    ids |= {e["id"] for e in excluded.get(aspect, [])}
    return frozenset(ids)


def add_exclude(excluded, aspect, go_id, name, reason):
    """Append one term to the per-aspect exclude list, if not already present."""
    bucket = excluded.setdefault(aspect, [])
    if not any(e["id"] == go_id for e in bucket):
        bucket.append({"id": go_id, "name": name, "reason": reason})
    return excluded


def remove_exclude(excluded, aspect, go_id):
    """Remove one term from the per-aspect exclude list (or global)."""
    for key in ("global", aspect):
        excluded[key] = [e for e in excluded.get(key, []) if e["id"] != go_id]
    return excluded


def audit_slim(slim_ids, godag, replaced_by, names, term_gene_counts, n_genes_total, depth):
    """Flag stock-slim terms that are obsolete, replaced, missing from the
    release OBO, uninformative (near-root or covering most of the gene
    list), or empty (map zero genes in the current gene list)."""
    rows = []
    for go_id in sorted(slim_ids):
        term = godag.get(go_id)
        flags = []
        replacement = None

        if go_id in replaced_by:
            flags.append("replaced")
            replacement = replaced_by[go_id]
        if term is None:
            flags.append("missing")
        elif term.is_obsolete:
            flags.append("obsolete")

        n_genes = term_gene_counts.get(go_id, 0)
        term_depth = depth.get(go_id, 0)
        if term_depth <= UNINFORMATIVE_DEPTH or n_genes_total and n_genes / n_genes_total > UNINFORMATIVE_COVERAGE_FRACTION:
            flags.append("uninformative")
        if n_genes == 0:
            flags.append("empty")

        if flags:
            rows.append({
                "go_id": go_id,
                "name": names.get(go_id, (term.name if term else None)),
                "flags": ";".join(flags),
                "replacement": replacement,
                "n_genes": n_genes,
                "depth": term_depth,
            })

    return pd.DataFrame(rows, columns=["go_id", "name", "flags", "replacement", "n_genes", "depth"])

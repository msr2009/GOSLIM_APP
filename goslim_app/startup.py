"""Build the in-memory DATA bundle once, at app import time."""

import time

import config
from ancestors import build_ancestor_closure, build_depth_map
from annotations import build_annotation_index
from aspects import ASPECT_NAMESPACE
from go_data import (
    build_replaced_by_map,
    load_gaf,
    load_genes,
    load_godag,
    load_slim_ids,
    load_tagged_genes,
)


def build_data(progress=print):
    """Load the GAF/OBO/slims/gene-lists and derive everything the app needs
    to recompute a grouping in milliseconds per iteration."""
    t0 = time.time()

    def phase(n, total, msg):
        progress(f"[{n}/{total}] {msg}... ({time.time() - t0:.1f}s elapsed)")

    phase(1, 6, f"Loading GO DAG ({config.OBO_PATH.name})")
    godag = load_godag(config.OBO_PATH)

    phase(2, 6, "Building obsolete-term replacement map")
    replaced_by = build_replaced_by_map(config.OBO_PATH)

    phase(3, 6, f"Loading GAF ({config.GAF_PATH.name})")
    gaf_df = load_gaf(config.GAF_PATH)

    phase(4, 6, "Building per-aspect annotation index")
    ann = build_annotation_index(gaf_df, replaced_by, progress=progress)
    del gaf_df  # the raw frame is ~300k rows; only the per-aspect dicts are kept

    phase(5, 6, "Loading base slims")
    base_slims = {}
    for slim_key, slim_path in config.SLIM_PATHS.items():
        for aspect, namespace in ASPECT_NAMESPACE.items():
            base_slims[(slim_key, aspect)] = frozenset(load_slim_ids(slim_path, namespace))

    phase(6, 6, "Building ancestor closure and gene lists")
    seed_ids = set()
    for aspect_ann in ann.values():
        for terms in aspect_ann.values():
            seed_ids |= terms
    for slim_ids in base_slims.values():
        seed_ids |= slim_ids
    anc = build_ancestor_closure(godag, seed_ids, use_part_of=config.USE_PART_OF_DEFAULT)
    depth = build_depth_map(godag, seed_ids)
    names = {go_id: term.name for go_id, term in godag.items()}

    gene_lists = {}
    for key, path in config.GENE_LIST_PATHS.items():
        if path.exists():
            gene_lists[key] = load_genes(str(path))
        else:
            progress(f"WARNING: gene list '{key}' not found at {path}, skipping")

    symbol_lookup = _build_symbol_lookup(ann)
    tagged_genes, tagged_genes_cgc = load_tagged_genes(config.TAGGED_GENES_PATH)

    progress(f"Done. Total load time {time.time() - t0:.1f}s")

    return {
        "godag": godag,
        "replaced_by": replaced_by,
        "ann": ann,
        "anc": anc,
        "depth": depth,
        "names": names,
        "base_slims": base_slims,
        "gene_lists": gene_lists,
        "symbol_lookup": symbol_lookup,
        "tagged_genes": tagged_genes,
        "tagged_genes_cgc": tagged_genes_cgc,
    }


def _build_symbol_lookup(ann):
    """Casefolded symbol -> canonical symbol, across all three aspects, so
    gene-list lookups tolerate case differences."""
    lookup = {}
    for aspect_ann in ann.values():
        for gene in aspect_ann:
            lookup.setdefault(gene.casefold(), gene)
    return lookup


# built once at import; shared across all browser sessions
DATA = build_data()

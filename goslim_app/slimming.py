"""Core per-iteration recompute: map a gene list to slim groups."""

import pandas as pd
from ancestors import map_term_to_slim
from partition import assign_single_group


def compute_grouping(genes, aspect, slim_ids, data):
    """Map every gene in `genes` to its slim groups for one aspect.

    Returns a dict of DataFrames: per_pair, term_counts, partition,
    partition_counts, uncovered_terms, and a coverage summary dict.
    """
    ann = data["ann"][aspect]
    names = data["names"]
    anc = data["anc"]
    depth = data["depth"]
    slim_ids = frozenset(slim_ids)

    canonical_genes, unmatched = _resolve_symbols(genes, ann, data["symbol_lookup"])

    pair_rows = []
    uncovered_counter = {}  # go_id -> [n_uncovered_genes, n_total_genes]
    no_annotation = []
    annotated_unmapped = []

    for gene in canonical_genes:
        raw_terms = ann.get(gene)
        if not raw_terms:
            no_annotation.append(gene)
            continue

        hit_ids = set()
        for term in raw_terms:
            hits = map_term_to_slim(term, slim_ids, anc)
            hit_ids |= hits
            # tally every raw term toward the uncovered-terms panel, whether
            # or not it happened to map, so the panel always reflects the
            # true frequency of terms among genes not fully covered
            uncovered_counter.setdefault(term, [0, 0])[1] += 1

        if hit_ids:
            for slim_id in hit_ids:
                pair_rows.append({"gene": gene, "slim_id": slim_id, "slim_name": names.get(slim_id)})
        else:
            annotated_unmapped.append(gene)
            for term in raw_terms:
                uncovered_counter[term][0] += 1

    per_pair = pd.DataFrame(pair_rows, columns=["gene", "slim_id", "slim_name"])

    term_counts = _build_term_counts(per_pair, depth)
    partition, partition_counts = assign_single_group(
        per_pair, term_counts, no_annotation, annotated_unmapped, unmatched
    )
    partition_counts = _add_leaf_breakdown(partition, partition_counts, ann, anc)
    uncovered_terms = _build_uncovered_terms(uncovered_counter, names, depth)

    coverage = {
        "n_input": len(genes),
        "n_unmatched_symbol": len(unmatched),
        "n_no_annotation": len(no_annotation),
        "n_annotated_but_unmapped": len(annotated_unmapped),
        "n_covered": len(canonical_genes) - len(no_annotation) - len(annotated_unmapped),
    }

    return {
        "per_pair": per_pair,
        "term_counts": term_counts,
        "partition": partition,
        "partition_counts": partition_counts,
        "uncovered_terms": uncovered_terms,
        "coverage": coverage,
        "unmatched_symbols": unmatched,
    }


def _resolve_symbols(genes, ann, symbol_lookup):
    """Map user-provided gene symbols onto canonical GAF symbols, tolerating
    case differences. Returns (canonical_genes, unmatched_originals)."""
    canonical = []
    unmatched = []
    seen = set()
    for g in genes:
        if g in ann:
            match = g
        else:
            match = symbol_lookup.get(g.casefold())
        if match is None:
            unmatched.append(g)
        elif match not in seen:
            seen.add(match)
            canonical.append(match)
    return canonical, unmatched


def _add_leaf_breakdown(partition, partition_counts, ann, anc):
    """Split each group's n_genes into n_leaf (the assigned slim term is the
    gene's most specific annotation - can't be split any further) and
    n_splittable (the gene also carries a more specific raw term that isn't
    in the slim, so drilling down could pull it into a smaller group)."""
    if partition.empty:
        partition_counts["n_leaf"] = 0
        partition_counts["n_splittable"] = 0
        return partition_counts

    def is_leaf(row):
        for term in ann.get(row.gene, ()):
            if term != row.group_id and row.group_id in anc.get(term, ()):
                return False
        return True

    is_leaf_col = partition.apply(is_leaf, axis=1)
    leaf_counts = partition.loc[is_leaf_col].groupby("group_id")["gene"].nunique()

    partition_counts["n_leaf"] = partition_counts["group_id"].map(leaf_counts).fillna(0).astype(int)
    # leftover pseudo-rows (unannotated/unmapped/symbol-not-found) have no group_id and
    # aren't real slim terms, so "splitting" doesn't apply - count them all as leaf
    no_group = partition_counts["group_id"].isna()
    partition_counts.loc[no_group, "n_leaf"] = partition_counts.loc[no_group, "n_genes"]
    partition_counts["n_splittable"] = partition_counts["n_genes"] - partition_counts["n_leaf"]
    return partition_counts


def _build_term_counts(per_pair, depth):
    """One row per slim term actually used, with multi-membership gene
    counts and how many genes are exclusive to that term."""
    if per_pair.empty:
        return pd.DataFrame(columns=["slim_id", "slim_name", "n_genes", "n_genes_exclusive", "depth"])

    counts = per_pair.groupby(["slim_id", "slim_name"], as_index=False)["gene"].nunique()
    counts = counts.rename(columns={"gene": "n_genes"})

    genes_per_term = per_pair.groupby("gene")["slim_id"].nunique()
    exclusive_genes = genes_per_term[genes_per_term == 1].index
    exclusive_pairs = per_pair[per_pair["gene"].isin(exclusive_genes)]
    exclusive_counts = exclusive_pairs.groupby("slim_id")["gene"].nunique()
    counts["n_genes_exclusive"] = counts["slim_id"].map(exclusive_counts).fillna(0).astype(int)

    counts["depth"] = counts["slim_id"].map(depth).fillna(0).astype(int)
    return counts.sort_values("n_genes", ascending=False).reset_index(drop=True)


def _build_uncovered_terms(uncovered_counter, names, depth):
    """Frequency table of raw GO terms seen among genes that are not (fully)
    covered by the current slim, ordered by uncovered-gene count desc."""
    rows = [
        {
            "go_id": go_id,
            "go_name": names.get(go_id),
            "n_uncovered_genes": n_uncovered,
            "n_total_genes": n_total,
            "depth": depth.get(go_id, 0),
        }
        for go_id, (n_uncovered, n_total) in uncovered_counter.items()
        if n_uncovered > 0
    ]
    df = pd.DataFrame(rows, columns=["go_id", "go_name", "n_uncovered_genes", "n_total_genes", "depth"])
    return df.sort_values("n_uncovered_genes", ascending=False).reset_index(drop=True)

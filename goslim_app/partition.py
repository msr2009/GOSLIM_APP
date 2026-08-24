"""Forced single-group assignment: partition genes into exactly one slim
group each, for even lab-splitting, rather than the natural multi-membership."""

import pandas as pd

UNANNOTATED_LABEL = "— unannotated —"
UNMAPPED_LABEL = "— unmapped —"
UNMATCHED_LABEL = "— symbol not found —"


def assign_single_group(per_pair, term_counts, no_annotation, annotated_unmapped, unmatched_symbols=()):
    """Pick one slim group per covered gene: rarest candidate group wins,
    ties broken by deeper (more specific) term, then GO id lexically. Both
    orderings are computed from term_counts BEFORE any assignment, so the
    result does not depend on iteration order.

    Returns (partition_df, partition_counts_df); partition_counts sums to
    the full input gene count once the leftover pseudo-rows (unannotated,
    unmapped, symbol-not-found) are included.
    """
    if per_pair.empty:
        partition = pd.DataFrame(columns=["gene", "group_id", "group_name", "n_candidate_groups"])
    else:
        rank = {
            row.slim_id: (row.n_genes, -row.depth, row.slim_id)
            for row in term_counts.itertuples()
        }
        ranked = per_pair.assign(_rank_key=per_pair["slim_id"].map(rank))
        ranked = ranked.sort_values(["gene", "_rank_key"])

        n_candidates = ranked.groupby("gene")["slim_id"].transform("nunique")
        ranked = ranked.assign(n_candidate_groups=n_candidates)

        partition = (
            ranked.drop_duplicates("gene")
            .rename(columns={"slim_id": "group_id", "slim_name": "group_name"})
            [["gene", "group_id", "group_name", "n_candidate_groups"]]
            .reset_index(drop=True)
        )

    partition_counts = (
        partition.groupby(["group_id", "group_name"], as_index=False)["gene"]
        .nunique()
        .rename(columns={"gene": "n_genes"})
    )

    # append visible leftover rows so the column sums to the full gene total
    leftover_rows = []
    if no_annotation:
        leftover_rows.append({"group_id": None, "group_name": UNANNOTATED_LABEL, "n_genes": len(no_annotation)})
    if annotated_unmapped:
        leftover_rows.append({"group_id": None, "group_name": UNMAPPED_LABEL, "n_genes": len(annotated_unmapped)})
    if unmatched_symbols:
        leftover_rows.append({"group_id": None, "group_name": UNMATCHED_LABEL, "n_genes": len(unmatched_symbols)})
    if leftover_rows:
        partition_counts = pd.concat(
            [partition_counts, pd.DataFrame(leftover_rows)], ignore_index=True
        )

    partition_counts = partition_counts.sort_values("n_genes", ascending=False).reset_index(drop=True)
    return partition, partition_counts


def rebalance(partition, partition_counts):
    """Stub: a greedy pass that would move genes from oversized groups to an
    undersized candidate group (from n_candidate_groups > 1) until the max
    group size stops shrinking. Not implemented - rarest-wins already tends
    toward even groups; revisit only if that proves insufficient."""
    raise NotImplementedError

"""'Split' support: break one slim group down into the more specific raw GO
terms actually carried by its member genes, so an over-broad group (e.g.
1500+ genes under 'membrane') can be pulled apart by promoting a few of its
descendants into the working slim."""

import pandas as pd


def compute_term_breakdown(go_id, per_pair, ann, anc, names, depth):
    """For genes currently assigned to slim term `go_id` (multi-membership),
    tally the more specific raw GO terms they carry - any raw term for which
    `go_id` is an ancestor. Returns (breakdown_df, n_genes_in_group,
    n_genes_at_leaf), where n_genes_at_leaf counts genes whose most specific
    annotation *is* go_id itself (no more specific term to split into)."""
    if per_pair.empty:
        genes_in_group = set()
    else:
        genes_in_group = set(per_pair.loc[per_pair["slim_id"] == go_id, "gene"])

    genes_by_term = {}
    genes_with_more_specific_term = set()
    for gene in genes_in_group:
        for term in ann.get(gene, ()):
            if term == go_id:
                continue
            if go_id in anc.get(term, ()):
                genes_by_term.setdefault(term, set()).add(gene)
                genes_with_more_specific_term.add(gene)

    rows = [
        {"go_id": t, "go_name": names.get(t), "n_genes": len(genes), "depth": depth.get(t, 0)}
        for t, genes in genes_by_term.items()
    ]
    df = pd.DataFrame(rows, columns=["go_id", "go_name", "n_genes", "depth"])
    df = df.sort_values("n_genes", ascending=False).reset_index(drop=True)
    n_genes_at_leaf = len(genes_in_group - genes_with_more_specific_term)
    return df, len(genes_in_group), n_genes_at_leaf

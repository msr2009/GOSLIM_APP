"""Build downloadable TSV/text buffers from grouping results."""

import io

import pandas as pd


def build_full_gene_table(genes, grouping, ann, aspect):
    """One row per input gene: raw GO terms, slim terms, assigned group,
    and a status flag. Used for the 'download gene -> group TSV' button."""
    per_pair = grouping["per_pair"]
    partition = grouping["partition"]
    unmatched = set(grouping["unmatched_symbols"])

    slim_terms_by_gene = per_pair.groupby("gene")["slim_id"].apply(list).to_dict()
    slim_names_by_gene = per_pair.groupby("gene")["slim_name"].apply(list).to_dict()
    group_by_gene = partition.set_index("gene")[["group_id", "group_name"]].to_dict("index")

    aspect_ann = ann[aspect]
    rows = []
    for g in genes:
        if g in unmatched:
            status = "symbol_not_found"
        elif g not in aspect_ann or not aspect_ann[g]:
            status = "no_annotation"
        elif g in group_by_gene:
            status = "covered"
        else:
            status = "unmapped"

        group_info = group_by_gene.get(g, {})
        rows.append({
            "gene": g,
            "raw_go_terms": ";".join(sorted(aspect_ann.get(g, ()))),
            "slim_terms": ";".join(slim_terms_by_gene.get(g, [])),
            "slim_names": ";".join(n for n in slim_names_by_gene.get(g, []) if n),
            "assigned_group": group_info.get("group_id"),
            "assigned_group_name": group_info.get("group_name"),
            "status": status,
        })

    return pd.DataFrame(rows)


def df_to_tsv_bytes(df):
    """Serialize a DataFrame to TSV bytes for a download handler."""
    buf = io.StringIO()
    df.to_csv(buf, sep="\t", index=False)
    return buf.getvalue().encode("utf-8")


def df_to_xlsx_bytes(df, sheet_name="Sheet1"):
    """Serialize a DataFrame to XLSX bytes for a download handler."""
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    return buf.getvalue()


def gene_list_to_text_bytes(genes):
    """One gene symbol per line - round-trips as a new gene-list input."""
    return ("\n".join(sorted(genes)) + "\n").encode("utf-8")


def genes_for_group(per_pair, slim_id):
    """All genes mapped to one slim group (from the multi-membership table)."""
    return sorted(per_pair.loc[per_pair["slim_id"] == slim_id, "gene"].unique())


def genes_for_partition_group(partition, group_id):
    """All genes assigned to one group in the forced-partition table."""
    return sorted(partition.loc[partition["group_id"] == group_id, "gene"].unique())

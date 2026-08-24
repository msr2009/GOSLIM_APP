"""
Map a list of genes to their Cellular Component GO annotations and GO Slim terms.

Inputs:
  - GAF file (e.g. gene_association.WS298.wb)
  - GO ontology OBO file (e.g. gene_ontology.WS298.obo) - matched release to the GAF
  - GO Slim OBO file(s) to build a merged/custom slim (e.g. goslim_generic.obo, goslim_agr.obo)
  - A gene list (one gene symbol per line, or a Python list)

Outputs (written as <out-prefix>_per_term.tsv and <out-prefix>_per_gene.tsv):
  - per_term: one row per (gene, go_slim_term) - columns gene, go_term, go_slim_term, cellular_component
  - per_gene: one row per gene - columns gene, go_term, go_slim_term, cellular_component,
    with multiple slim terms/names joined by ";"

Usage:
  python map_genes_to_go_cc.py \
      --gaf gene_association.WS298.wb \
      --obo gene_ontology.WS298.obo \
      --slim goslim_generic.obo goslim_agr.obo \
      --genes gene_list.txt \
      --out-prefix worm_disease_go
"""

import argparse
import sys
from collections import defaultdict

import pandas as pd
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim


# Extra terms to inject into the merged slim because neither goslim_generic
# nor goslim_agr contained them, despite being very common CC annotations
# in our worm GAF (see conversation history for how these were identified).
EXTRA_SLIM_TERMS = [
    "GO:0016020",  # membrane
    "GO:0005737",  # cytoplasm
    "GO:0043005",  # neuron projection (covers axon/dendrite via is_a)
    "GO:0045202",  # synapse
    "GO:0030054",  # cell junction (covers gap junction via is_a)
]

# The bare CC root term - carries no real localization information.
# Genes whose only annotation is this term are treated as "no informative
# CC data", same bucket as genes with zero CC annotations.
GO_CC_ROOT = "GO:0005575"

GAF_COLS = [
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
    "DB_Reference", "Evidence_Code", "With_From", "Aspect",
    "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
    "Taxon", "Date", "Assigned_By", "Annotation_Extension",
    "Gene_Product_Form_ID",
]

QUALIFIERS_TO_KEEP = {"located_in", "part_of"}


def build_replaced_by_map(obo_path):
    """Parse the obo file directly to build {obsolete_id: replacement_id}."""
    replaced_by = {}
    current_id = None
    with open(obo_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split("id: ")[1]
            elif line.startswith("replaced_by:"):
                replaced_by[current_id] = line.split("replaced_by: ")[1].strip()
    return replaced_by


def build_custom_slim_dag(slim_paths, godag, extra_terms):
    """Merge one or more slim OBO files into a single GODag, then patch in
    any extra terms (pulled from the full godag) that are missing."""
    if not slim_paths:
        raise ValueError("At least one --slim file is required")

    # Start from the first slim file as the base GODag object
    custom_slim_dag = GODag(slim_paths[0], optional_attrs=["relationship"])

    # Merge in terms from any additional slim files
    for path in slim_paths[1:]:
        extra_dag = GODag(path, optional_attrs=["relationship"])
        for go_id, term in extra_dag.items():
            if go_id not in custom_slim_dag:
                custom_slim_dag[go_id] = term

    # Patch in manually-identified missing terms
    for go_id in extra_terms:
        if go_id in godag and go_id not in custom_slim_dag:
            custom_slim_dag[go_id] = godag[go_id]

    return custom_slim_dag


def load_gaf(gaf_path):
    return pd.read_csv(
        gaf_path,
        sep="\t",
        comment="!",
        header=None,
        names=GAF_COLS,
        dtype=str,
    )


def load_genes(genes_arg):
    """genes_arg can be a path to a file (genes separated by newlines, tabs,
    commas, or any mix of whitespace), or a delimited string passed directly
    on the command line."""
    import os
    import re

    if os.path.isfile(genes_arg):
        with open(genes_arg) as f:
            text = f.read()
    else:
        text = genes_arg

    # Split on any run of whitespace, tabs, newlines, and/or commas
    genes = [g.strip() for g in re.split(r"[,\s]+", text) if g.strip()]
    return genes


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaf", required=True, help="Path to GAF file")
    parser.add_argument("--obo", required=True, help="Path to GO ontology OBO file (matched release to GAF)")
    parser.add_argument("--slim", required=True, nargs="+", help="One or more GO Slim OBO files to merge")
    parser.add_argument("--genes", required=True, help="Path to gene list file (one per line) or comma-separated list")
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Prefix for output TSVs; writes <prefix>_per_term.tsv and <prefix>_per_gene.tsv",
    )
    parser.add_argument(
        "--keep-root",
        action="store_true",
        help="Keep the bare GO cellular_component root term as its own slim category "
             "instead of treating it as missing/uninformative data (default: fold into None)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("Loading gene ontology DAG...", file=sys.stderr)
    godag = GODag(args.obo, optional_attrs=["relationship"])

    print("Building replaced_by map for obsolete terms...", file=sys.stderr)
    replaced_by = build_replaced_by_map(args.obo)

    print("Building merged custom slim DAG...", file=sys.stderr)
    custom_slim_dag = build_custom_slim_dag(args.slim, godag, EXTRA_SLIM_TERMS)

    print("Loading GAF...", file=sys.stderr)
    worm_go = load_gaf(args.gaf)

    print("Loading gene list...", file=sys.stderr)
    genes = load_genes(args.genes)
    print(f"{len(genes)} genes loaded", file=sys.stderr)

    def map_to_slim(go_id):
        go_id = replaced_by.get(go_id, go_id)  # substitute obsolete-with-replacement terms
        try:
            direct_ancestors, _all_ancestors = mapslim(go_id, godag, custom_slim_dag)
            return direct_ancestors
        except Exception:
            return set()

    def retrieve_go_name(go_id):
        try:
            return godag[go_id].name
        except KeyError:
            return None

    rows = []
    for g in genes:
        go_term_cc = worm_go[
            (worm_go["DB_Object_Symbol"] == g)
            & (worm_go["Qualifier"].isin(QUALIFIERS_TO_KEEP))
        ]["GO_ID"]

        if len(go_term_cc) == 0:
            rows.append({
                "gene": g,
                "go_term": None,
                "go_slim_term": None,
                "cellular_component": None,
            })
            continue

        for t in go_term_cc:
            slim_term = map_to_slim(t)

            if len(slim_term) == 0:
                rows.append({
                    "gene": g,
                    "go_term": t,
                    "go_slim_term": None,
                    "cellular_component": None,
                })
            else:
                st = slim_term.pop()

                # Fold the bare CC root term into "no informative data"
                # unless the user explicitly wants to keep it
                if st == GO_CC_ROOT and not args.keep_root:
                    rows.append({
                        "gene": g,
                        "go_term": t,
                        "go_slim_term": None,
                        "cellular_component": None,
                    })
                else:
                    rows.append({
                        "gene": g,
                        "go_term": t,
                        "go_slim_term": st,
                        "cellular_component": retrieve_go_name(st),
                    })

    result = pd.DataFrame(
        rows, columns=["gene", "go_term", "go_slim_term", "cellular_component"]
    )

    # Deduplicate: one row per unique (gene, go_slim_term) rather than one row
    # per raw go_term, so multiple specific terms collapsing to the same slim
    # bucket don't produce redundant rows for the same gene.
    per_term = (
        result
        .groupby(["gene", "go_slim_term", "cellular_component"], as_index=False, dropna=False)
        .agg(go_term=("go_term", lambda x: ";".join(sorted(set(v for v in x if pd.notna(v))))))
    )
    per_term = per_term[["gene", "go_term", "go_slim_term", "cellular_component"]]

    # Collapse further to one row per gene, joining multiple slim terms/names
    # (and their source go_terms) with ";" so every gene appears exactly once.
    per_gene = (
        per_term
        .groupby("gene", as_index=False)
        .agg(
            go_term=("go_term", lambda x: ";".join(sorted(set(v for v in x if pd.notna(v) and v)))),
            go_slim_term=("go_slim_term", lambda x: ";".join(sorted(set(v for v in x if pd.notna(v))))),
            cellular_component=("cellular_component", lambda x: ";".join(sorted(set(v for v in x if pd.notna(v))))),
        )
    )
    per_gene = per_gene[["gene", "go_term", "go_slim_term", "cellular_component"]].replace("", None)

    per_term_out = f"{args.out_prefix}_per_term.tsv"
    per_gene_out = f"{args.out_prefix}_per_gene.tsv"

    per_term.to_csv(per_term_out, sep="\t", index=False)
    print(f"Wrote {len(per_term)} rows to {per_term_out}", file=sys.stderr)

    per_gene.to_csv(per_gene_out, sep="\t", index=False)
    print(f"Wrote {len(per_gene)} rows to {per_gene_out}", file=sys.stderr)

    n_genes_no_data = per_term[per_term["go_slim_term"].isna()]["gene"].nunique()
    print(f"{n_genes_no_data} / {len(genes)} genes have no informative CC slim annotation", file=sys.stderr)


if __name__ == "__main__":
    main()

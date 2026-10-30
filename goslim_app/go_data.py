"""Raw file loaders: GAF, obsolete-term map, slim OBO term ids, gene lists."""

import os
import re

import pandas as pd
from goatools.obo_parser import GODag

GAF_COLS = [
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
    "DB_Reference", "Evidence_Code", "With_From", "Aspect",
    "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
    "Taxon", "Date", "Assigned_By", "Annotation_Extension",
    "Gene_Product_Form_ID",
]

# only these columns are needed downstream - skip the rest to save memory
GAF_USECOLS = ["DB_Object_Symbol", "Qualifier", "GO_ID", "Aspect", "Taxon"]


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


def load_godag(obo_path):
    """Load the full release-matched GO DAG, including part_of relationships."""
    return GODag(str(obo_path), optional_attrs=["relationship"])


def load_slim_ids(slim_path, aspect_namespace):
    """Return the set of GO ids in a slim OBO file that belong to one aspect
    namespace. Parses the OBO directly rather than building a second GODag,
    since only ids + namespace are needed."""
    ids = set()
    current_id = None
    current_namespace = None
    with open(slim_path) as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                if current_id and current_namespace == aspect_namespace:
                    ids.add(current_id)
                current_id = None
                current_namespace = None
            elif line.startswith("id: GO:"):
                current_id = line.split("id: ")[1]
            elif line.startswith("namespace:"):
                current_namespace = line.split("namespace: ")[1].strip()
    # flush the last stanza in the file
    if current_id and current_namespace == aspect_namespace:
        ids.add(current_id)
    return ids


def load_gaf(gaf_path):
    """Load the GAF, keeping only the columns the app needs."""
    return pd.read_csv(
        gaf_path,
        sep="\t",
        comment="!",
        header=None,
        names=GAF_COLS,
        usecols=GAF_USECOLS,
        dtype=str,
    )


def load_genes(genes_arg):
    """genes_arg can be a path to a file (genes separated by newlines, tabs,
    commas, or any mix of whitespace), or a delimited string passed directly."""
    if os.path.isfile(genes_arg):
        with open(genes_arg) as f:
            text = f.read()
    else:
        text = genes_arg

    genes = [g.strip() for g in re.split(r"[,\s]+", text) if g.strip()]
    return genes


def load_gene_classes(path, symbol_lookup=None):
    """Long-form gene<sep>class file -> {class: frozenset(genes)}, in first-seen class order.
    Tab-separated unless the name ends in .csv; header row optional; extra columns ignored.
    When symbol_lookup (casefolded -> canonical symbol) is given, genes are canonicalized."""
    if not path or not os.path.isfile(path):
        return {}
    sep = "," if str(path).lower().endswith(".csv") else "\t"
    df = pd.read_csv(
        path, sep=sep, header=None, usecols=[0, 1], names=["gene", "class"],
        dtype=str, keep_default_na=False, skip_blank_lines=True,
    )
    df = df.apply(lambda col: col.str.strip())
    # drop an optional header row and rows missing either field
    if len(df) and (df.iloc[0].str.casefold() == ["gene", "class"]).all():
        df = df.iloc[1:]
    df = df[(df["gene"] != "") & (df["class"] != "")]
    if symbol_lookup:
        df["gene"] = [symbol_lookup.get(g.casefold(), g) for g in df["gene"]]
    return {name: frozenset(grp["gene"]) for name, grp in df.groupby("class", sort=False)}

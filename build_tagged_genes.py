"""Convert a WormTagDB 'filtered strains' CSV export into wormtagdb_tagged_genes.tsv."""

import argparse
import re
from datetime import date
from pathlib import Path

import pandas as pd

REPO_DIR = Path(__file__).resolve().parent
FETCH_DATE_PATH = REPO_DIR / "goslim_app" / "data_fetched.txt"
OUT_COLS = ["gene", "wbid", "n_strains", "fluor_tags", "other_tags", "cgc_strains"]
FLUOR_CLASS = "WormTagDB fluorescent tags"
CGC_CLASS = "CGC-available fluorescent tags"


def strip_html(value):
    """Remove HTML tags (the export wraps IDs in <a> links) and nbsp padding."""
    return re.sub(r"<[^>]+>", "", value).replace("\xa0", " ").strip()


def split_tags(value):
    """Split a comma-joined tag cell into a list, dropping NA/empty entries."""
    tags = [t.strip() for t in value.split(",")]
    return [t for t in tags if t and t != "NA"]


def join_unique(values):
    """Sorted unique non-empty values joined with ';'."""
    return ";".join(sorted({v for v in values if v}))


def collapse_by_gene(df):
    """One row per gene from the per-strain export (tags unioned, CGC strains collected)."""
    df = df.copy()
    df["Gene"] = df["Gene"].map(strip_html)
    df["wbid"] = df["WormBase ID"].map(strip_html)
    # Only strains linked to cgc.umn.edu are CGC-orderable; 'NA' means no deposited strain
    df["cgc"] = [strip_html(s) if "cgc.umn.edu" in s else "" for s in df["Strain"].astype(str)]
    rows = []
    for gene, grp in df.groupby("Gene", sort=True):
        fluor = [t for cell in grp["Fluor"] for t in split_tags(cell)]
        other = [t for cell in grp["Other tags"] for t in split_tags(cell)]
        rows.append(
            {
                "gene": gene,
                "wbid": join_unique(grp["wbid"]),
                "n_strains": len(grp),
                "fluor_tags": join_unique(fluor),
                "other_tags": join_unique(other),
                "cgc_strains": join_unique(grp["cgc"]),
            }
        )
    return pd.DataFrame(rows, columns=OUT_COLS)


def to_long_form(wide):
    """Per-gene table -> long-form (gene, class) rows for the app's classification input."""
    fluor = wide[wide["fluor_tags"] != ""]
    cgc = fluor[fluor["cgc_strains"] != ""]
    parts = [
        pd.DataFrame({"gene": fluor["gene"], "class": FLUOR_CLASS}),
        pd.DataFrame({"gene": cgc["gene"], "class": CGC_CLASS}),
    ]
    return pd.concat(parts, ignore_index=True)


def report_changes(old, new):
    """Print genes added/removed and fluor/CGC status changes between two gene tables."""
    old_genes, new_genes = set(old["gene"]), set(new["gene"])
    print(f"genes: {len(old_genes)} -> {len(new_genes)}")
    print(f"  added:   {len(new_genes - old_genes)}")
    print(f"  removed: {len(old_genes - new_genes)}")
    merged = old.merge(new, on="gene", suffixes=("_old", "_new"))
    for col in ["fluor_tags", "other_tags", "cgc_strains"]:
        n = (merged[f"{col}_old"] != merged[f"{col}_new"]).sum()
        print(f"  {col} changed: {n}")


def fetch_date(csv_path):
    """Export date from a filtered_strains_YYYY-MM-DD*.csv name, else today."""
    m = re.search(r"(\d{4})-(\d{2})-(\d{2})", Path(csv_path).name)
    return date(*map(int, m.groups())) if m else date.today()


def main():
    """Parse args, collapse the export to per-gene rows, report changes, write the TSV."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", help="WormTagDB export (filtered_strains_*.csv)")
    ap.add_argument("--out", default=str(REPO_DIR / "wormtagdb_tagged_genes.tsv"))
    ap.add_argument("--classes-out", default=str(REPO_DIR / "gene_classes.tsv"))
    ap.add_argument("--dry-run", action="store_true", help="report changes without writing")
    args = ap.parse_args()

    df = pd.read_csv(args.csv, keep_default_na=False, dtype=str)
    new = collapse_by_gene(df)

    out = Path(args.out)
    if out.is_file():
        old = pd.read_csv(out, sep="\t", keep_default_na=False, dtype=str)
        report_changes(old, new)
    if args.dry_run:
        return
    new.to_csv(out, sep="\t", index=False)
    print(f"wrote {len(new)} genes -> {out}")
    long_form = to_long_form(new)
    long_form.to_csv(args.classes_out, sep="\t", index=False)
    print(f"wrote {len(long_form)} gene/class rows -> {args.classes_out}")
    # the app shows this date in its title bar
    FETCH_DATE_PATH.write_text(fetch_date(args.csv).isoformat() + "\n")
    print(f"wrote fetch date -> {FETCH_DATE_PATH}")


if __name__ == "__main__":
    main()

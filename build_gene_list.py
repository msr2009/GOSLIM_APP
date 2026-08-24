"""
One-off script: download the WormBase geneIDs file for C. elegans and derive
data/celegans_proteome.tsv, a plain list of protein-coding gene symbols.

The WormBase downloads site (downloads.wormbase.org) sits behind a Cloudflare
challenge that blocks plain curl/requests, but the EBI FTP mirror serves the
same release tree without it - use that instead.

Usage:
  mamba run -n go_terms python build_gene_list.py --release WS298
"""

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

EBI_MIRROR = (
    "https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/{release}/species/"
    "c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.{release}.geneIDs.txt.gz"
)


def download_gene_ids(release, dest_gz):
    """Fetch the gzipped WormBase geneIDs file from the EBI mirror."""
    url = EBI_MIRROR.format(release=release)
    print(f"Downloading {url}", file=sys.stderr)
    urllib.request.urlretrieve(url, dest_gz)


def extract_protein_coding_symbols(gz_path):
    """Parse the geneIDs CSV and return sorted public-name symbols for Live
    protein-coding genes, falling back to the sequence name when no public
    name is assigned."""
    symbols = set()
    with gzip.open(gz_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split(",")
            if len(fields) < 6:
                continue
            _taxon, _wbgene, public_name, seq_name, status, biotype = fields[:6]
            if status != "Live" or biotype != "protein_coding_gene":
                continue
            symbols.add(public_name if public_name else seq_name)
    return sorted(symbols)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", default="WS298", help="WormBase release tag")
    parser.add_argument(
        "--out", default="data/celegans_proteome.tsv", help="Output gene list path"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)
    gz_path = data_dir / f"celegans_geneIDs.{args.release}.txt.gz"

    download_gene_ids(args.release, gz_path)
    symbols = extract_protein_coding_symbols(gz_path)

    with open(args.out, "w") as f:
        f.write("\n".join(symbols) + "\n")
    print(f"Wrote {len(symbols)} protein-coding gene symbols to {args.out}", file=sys.stderr)

    gz_path.unlink()  # don't keep the raw download around


if __name__ == "__main__":
    main()

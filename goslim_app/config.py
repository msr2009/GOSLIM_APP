"""Paths and tunable settings for the GO slim explorer app."""

import json
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
REPO_DIR = APP_DIR.parent

# locally, data files live one level up in the repo; a cloud deploy (shinyapps.io)
# only bundles this app directory, so `prepare_deploy.sh` stages copies under
# APP_DIR/vendor with the same layout - use that instead when it's present
DATA_ROOT = APP_DIR / "vendor" if (APP_DIR / "vendor").is_dir() else REPO_DIR

GAF_PATH = DATA_ROOT / "gene_association.WS298.wb"
OBO_PATH = DATA_ROOT / "gene_ontology.WS298.obo"
SLIM_PATHS = {
    "agr": DATA_ROOT / "goslim_agr.obo",
    "generic": DATA_ROOT / "goslim_generic.obo",
}

GENE_LIST_DIR = DATA_ROOT / "data"

# every *.tsv/*.txt file under GENE_LIST_DIR is auto-loaded as a selectable gene
# list; key is the filename stem, label prettifies underscores for the dropdown
GENE_LIST_PATHS = {
    path.stem: path
    for path in sorted(GENE_LIST_DIR.glob("*"))
    if path.suffix in (".tsv", ".txt") and path.is_file()
}
GENE_LIST_LABELS = {stem: stem.replace("_", " ") for stem in GENE_LIST_PATHS}
DEFAULT_GENE_LIST = "celegans_proteome" if "celegans_proteome" in GENE_LIST_PATHS else next(
    iter(GENE_LIST_PATHS), None
)

SAVED_SLIMS_DIR = APP_DIR / "saved_slims"
EXCLUDE_FILE = DATA_ROOT / "excluded_terms.json"

# gene symbols with an endogenously tagged allele in WormTagDB (wormtagdb.rc.duke.edu);
# kept outside GENE_LIST_DIR since it has extra columns (wbid/n_strains/tags), not a
# plain one-symbol-per-line gene list
TAGGED_GENES_PATH = DATA_ROOT / "wormtagdb_tagged_genes.tsv"

# C. elegans taxon id as it appears in the GAF Taxon column
TAXON = "taxon:6239"

# whether the ancestor closure follows part_of edges in addition to is_a
USE_PART_OF_DEFAULT = True


def load_overrides():
    """Read optional goslim_app.json in the repo root to override any of the
    above paths/settings without editing this file."""
    override_path = REPO_DIR / "goslim_app.json"
    if not override_path.exists():
        return {}
    with open(override_path) as f:
        return json.load(f)

"""Turn the raw GAF into a per-aspect {gene: frozenset(go_ids)} index."""

from aspects import ASPECT_QUALIFIERS, ASPECT_ROOTS
from config import TAXON


def build_annotation_index(gaf_df, replaced_by, progress=print):
    """Reduce the GAF DataFrame to {aspect: {gene: frozenset(go_ids)}},
    filtered to C. elegans, positive qualifiers only, obsolete ids replaced,
    and the bare aspect root stripped out (uninformative)."""
    # only rows tagged for our taxon - the GAF spans many nematode species
    # and gene symbols collide across species without this filter
    taxon_rows = gaf_df["Taxon"].str.split("|").str[0] == TAXON
    df = gaf_df[taxon_rows]

    _warn_on_unexpected_qualifiers(df, progress)

    index = {}
    for aspect, qualifiers in ASPECT_QUALIFIERS.items():
        aspect_df = df[(df["Aspect"] == aspect) & (df["Qualifier"].isin(qualifiers))]
        root = ASPECT_ROOTS[aspect]

        gene_terms = {}
        for gene, go_id in zip(aspect_df["DB_Object_Symbol"], aspect_df["GO_ID"]):
            go_id = replaced_by.get(go_id, go_id)  # substitute obsolete ids
            if go_id == root:
                continue  # the bare root carries no informative signal
            gene_terms.setdefault(gene, set()).add(go_id)

        index[aspect] = {g: frozenset(terms) for g, terms in gene_terms.items()}

    return index


def _warn_on_unexpected_qualifiers(df, progress):
    """Flag any qualifier value not covered by ASPECT_QUALIFIERS or a NOT|
    negation, so a future GAF release doesn't silently drop annotations."""
    known = set()
    for qualifiers in ASPECT_QUALIFIERS.values():
        known |= qualifiers
    seen = set(df["Qualifier"].dropna().unique())
    unexpected = {q for q in seen if q not in known and not q.startswith("NOT")}
    if unexpected:
        progress(f"WARNING: unrecognized GAF qualifiers seen: {sorted(unexpected)}")

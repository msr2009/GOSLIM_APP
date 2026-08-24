"""GO aspect tables: qualifiers, roots, namespaces, and curated extra terms."""

ASPECT_LABELS = {
    "C": "Cellular Component",
    "F": "Molecular Function",
    "P": "Biological Process",
}

ASPECT_NAMESPACE = {
    "C": "cellular_component",
    "F": "molecular_function",
    "P": "biological_process",
}

# GAF Qualifier values that count as a positive, physically-informative
# annotation for each aspect. NOT|* qualifiers are excluded separately.
ASPECT_QUALIFIERS = {
    "C": {"located_in", "part_of"},
    "F": {"enables", "contributes_to"},
    "P": {
        "involved_in",
        "acts_upstream_of",
        "acts_upstream_of_positive_effect",
        "acts_upstream_of_negative_effect",
        "acts_upstream_of_or_within",
        "acts_upstream_of_or_within_positive_effect",
        "acts_upstream_of_or_within_negative_effect",
    },
}

# the bare aspect root term - carries no localization/function/process info
ASPECT_ROOTS = {
    "C": "GO:0005575",
    "F": "GO:0003674",
    "P": "GO:0008150",
}

# terms observed to be common in the worm GAF but missing from both stock
# slims; curated by hand while building the original CC-only script
EXTRA_SLIM_TERMS = {
    "C": [
        "GO:0016020",  # membrane
        "GO:0005737",  # cytoplasm
        "GO:0043005",  # neuron projection
        "GO:0045202",  # synapse
        "GO:0030054",  # cell junction
    ],
    "F": [],
    "P": [],
}

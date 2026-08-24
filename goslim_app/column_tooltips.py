"""Hover tooltips for data grid column headers.

Shiny's DataGrid renders plain HTML <th> elements with no built-in per-column
tooltip support, so this injects a small client-side script that sets the
`title` attribute on any header whose text matches one of the labels below.
One shared mapping covers every grid in the app, since the same column name
means the same thing wherever it appears.
"""

import json

from shiny import ui

COLUMN_TOOLTIPS = {
    "gene": "Gene symbol, matched against the GAF for the selected aspect.",
    "slim_id": "GO id of the slim (group) term.",
    "slim_name": "Human-readable name of the slim term.",
    "slim_terms": "Slim term ids this gene maps to (semicolon-separated).",
    "n_genes": "Genes mapped to this term. In the multi-membership table a gene "
               "may count toward more than one term, so this column can sum to "
               "more than the total gene count.",
    "n_genes_exclusive": "Genes for which this is the ONLY slim term they map to.",
    "depth": "Depth of this term in the GO DAG - larger means more specific.",
    "group_id": "GO id of the term this gene was assigned to under the forced "
                "single-group partition.",
    "group_name": "Name of the assigned partition group.",
    "n_candidate_groups": "How many slim terms this gene could have mapped to "
                          "before the rarest-wins tie-break picked one.",
    "go_id": "GO id of the raw (unslimmed) annotation term.",
    "go_name": "Human-readable name of the raw GO term.",
    "raw_go_terms": "All raw GO ids annotated to this gene in the current aspect.",
    "raw_go_names": "Names of the raw GO terms.",
    "n_uncovered_genes": "Genes not covered by the current slim that carry this "
                         "raw term - promote high-count rows to grow coverage.",
    "n_total_genes": "Total genes (covered or not) that carry this raw term.",
    "assigned_group": "The gene's single assigned group under the forced partition.",
    "assigned_group_name": "Name of the gene's assigned partition group.",
    "flags": "Why this term was flagged: obsolete, replaced, missing from the "
             "release OBO, uninformative (near-root or covers most genes), or "
             "empty (maps zero genes in the current gene list).",
    "replacement": "If flagged 'replaced', the GO id this term should be "
                   "substituted with.",
    "name": "Human-readable GO term name.",
    "reason": "Why this term is on the exclude list.",
    "status": "covered / unmapped (annotated but no slim term matches) / "
              "no_annotation (never accessible in this aspect) / "
              "symbol_not_found (gene symbol not in the GAF).",
}


def tooltip_script():
    """A script tag that keeps column-header title attributes in sync as
    Shiny's data grids re-render (on sort, filter, or data change)."""
    mapping_json = json.dumps(COLUMN_TOOLTIPS)
    return ui.tags.script(f"""
        (function() {{
            var tooltips = {mapping_json};
            function applyTooltips() {{
                document.querySelectorAll('table thead th').forEach(function(th) {{
                    var label = th.textContent.trim();
                    var tip = tooltips[label];
                    if (tip) {{
                        th.title = tip;
                        var div = th.querySelector('div');
                        if (div) {{ div.title = tip; }}
                    }}
                }});
            }}
            var observer = new MutationObserver(applyTooltips);
            observer.observe(document.body, {{childList: true, subtree: true}});
            document.addEventListener('DOMContentLoaded', applyTooltips);
            applyTooltips();
        }})();
    """)

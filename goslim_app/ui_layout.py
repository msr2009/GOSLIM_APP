"""App UI: sidebar controls + tabbed main panel."""

import config
from aspects import ASPECT_LABELS
from column_tooltips import tooltip_script
from shiny import ui
from shinywidgets import output_widget

ASPECT_CHOICES = {k: v for k, v in ASPECT_LABELS.items()}
# static choices always present in the "Base slim" dropdown; saved slims for
# the current aspect are appended dynamically (see server_main._refresh_base_slim_choices)
STATIC_BASE_SLIM_CHOICES = {"agr": "GO slim AGR", "generic": "GO slim generic", "none": "Empty (build from scratch)"}
GENE_LIST_CHOICES = config.GENE_LIST_LABELS

sidebar = ui.sidebar(
    ui.input_select("aspect", "GO aspect", ASPECT_CHOICES, selected="C"),
    ui.input_select("base_slim", "Base slim", STATIC_BASE_SLIM_CHOICES, selected="agr"),
    ui.input_select("gene_list", "Gene list", GENE_LIST_CHOICES, selected=config.DEFAULT_GENE_LIST),
    ui.input_file("gene_upload", "...or upload a gene list", accept=[".txt", ".tsv", ".csv"]),
    ui.input_checkbox("use_extra_terms", "Include curated extra CC terms", True),
    ui.input_checkbox("use_part_of", "Follow part_of edges (not just is_a)", config.USE_PART_OF_DEFAULT),
    ui.hr(),
    ui.output_text("slim_summary"),
    ui.download_button("dl_working_slim_xlsx", "Download working slim (.xlsx)", class_="btn-sm"),
    ui.input_action_button("btn_reset", "Reset to base slim", class_="btn-sm"),
    ui.hr(),
    ui.input_select("load_slim_choice", "Load saved slim", {}),
    ui.input_action_button("btn_load", "Load", class_="btn-sm"),
    ui.hr(),
    ui.input_text("save_name", "Save current slim as"),
    ui.input_text_area("save_notes", "Notes", rows=2),
    ui.input_action_button("btn_save", "Save slim", class_="btn-sm btn-primary"),
    ui.hr(),
    ui.download_button("dl_full_table", "Download gene -> group TSV"),
    width=340,
)

coverage_tab = ui.nav_panel(
    "Coverage",
    ui.layout_column_wrap(
        ui.value_box("Input genes", ui.output_text("vb_input")),
        ui.value_box("Covered", ui.output_text("vb_covered")),
        ui.value_box("Genes without annotation in current slim mapping", ui.output_text("vb_unmapped")),
        ui.value_box("No annotation (never accessible)", ui.output_text("vb_no_annotation")),
        ui.value_box("Symbol not found", ui.output_text("vb_unmatched")),
        width=1 / 5,
    ),
    output_widget("group_size_plot"),
)

slim_groups_tab = ui.nav_panel(
    "Slim groups",
    ui.layout_columns(
        ui.card(
            ui.card_header("Multi-membership counts"),
            ui.output_data_frame("term_counts_grid"),
        ),
        ui.card(
            ui.card_header("Forced single-group partition"),
            ui.output_data_frame("partition_counts_grid"),
        ),
    ),
    ui.card(
        ui.card_header("Genes in selected group"),
        ui.output_text("selected_group_label"),
        ui.download_button("dl_group_genes", "Download this group's genes"),
        ui.output_data_frame("selected_group_genes_grid"),
    ),
)

add_terms_tab = ui.nav_panel(
    "Add terms",
    ui.card(
        ui.card_header("GO terms among uncovered genes, by frequency"),
        ui.input_action_button("btn_promote", "Add selected terms to slim", class_="btn-sm btn-primary"),
        ui.output_data_frame("uncovered_terms_grid"),
    ),
    ui.layout_columns(
        ui.card(
            ui.card_header("Current working slim (select a row to inspect its split, or remove)"),
            ui.input_action_button("btn_remove", "Remove selected terms", class_="btn-sm btn-danger"),
            ui.output_data_frame("working_slim_grid"),
        ),
        ui.card(
            ui.card_header("Split selected term"),
            ui.output_text("drill_label"),
            ui.input_action_button("btn_add_subterms", "Add selected subterms to slim", class_="btn-sm btn-primary"),
            ui.output_data_frame("subterms_grid"),
        ),
    ),
)

gene_detail_tab = ui.nav_panel(
    "Gene detail",
    ui.input_text("gene_search", "Filter genes"),
    ui.output_data_frame("gene_detail_grid"),
)

slim_hygiene_tab = ui.nav_panel(
    "Slim hygiene",
    ui.card(
        ui.card_header("Base slim audit (obsolete / replaced / missing / uninformative / empty)"),
        ui.input_action_button("btn_exclude", "Exclude selected terms", class_="btn-sm btn-danger"),
        ui.output_data_frame("audit_grid"),
    ),
    ui.card(
        ui.card_header("Current exclude list"),
        ui.input_action_button("btn_unexclude", "Remove selected from exclude list", class_="btn-sm"),
        ui.output_data_frame("exclude_grid"),
    ),
)

app_ui = ui.page_sidebar(
    sidebar,
    ui.busy_indicators.use(spinners=True, pulse=True),
    tooltip_script(),
    ui.navset_card_tab(
        coverage_tab,
        slim_groups_tab,
        add_terms_tab,
        gene_detail_tab,
        slim_hygiene_tab,
    ),
    title="GO Slim Explorer",
)

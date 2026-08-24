"""Reactive server: wires sidebar inputs -> working slim -> grouping -> outputs."""

import datetime

import config
import excluded as excluded_mod
import pandas as pd
import plotly.express as px
import saved_slims
from aspects import EXTRA_SLIM_TERMS
from drilldown import compute_term_breakdown
from exports import (
    build_full_gene_table,
    df_to_tsv_bytes,
    df_to_xlsx_bytes,
    gene_list_to_text_bytes,
    genes_for_group,
    genes_for_partition_group,
)
from go_data import load_genes
from shiny import reactive, render, ui
from shinywidgets import render_widget
from slimming import compute_grouping
from startup import DATA
from ui_layout import STATIC_BASE_SLIM_CHOICES

SAVED_BASE_PREFIX = "saved:"


def _go_id_cell(go_id, slim_ids):
    """go_id text, with a small blue 'in goslim' tag appended when go_id is
    already part of the working slim."""
    if go_id not in slim_ids:
        return go_id
    return ui.tags.span(
        go_id,
        ui.tags.span(" in goslim", style="color: #2563eb; font-size: 0.75em;"),
    )


def _go_id_value(cell):
    """Undo _go_id_cell: recover the plain go_id string from a grid cell that
    may be either the raw string or the wrapped '<go_id> in goslim' tag."""
    return cell.children[0] if hasattr(cell, "children") else cell


def server(input, output, session):
    # {"added": set(go_id), "removed": set(go_id)} - the delta on top of the base slim
    slim_state = reactive.value({"added": set(), "removed": set()})
    excluded_state = reactive.value(excluded_mod.load_excluded(config.EXCLUDE_FILE))
    selected_group = reactive.value(None)  # (source, go_id) where source is "multi" or "partition"
    drill_term = reactive.value(None)  # go_id currently selected for the "split" breakdown
    slim_files_version = reactive.value(0)  # bumped after every save, to refresh dropdowns

    # --- reset the working slim whenever aspect/base slim/extra-terms change ---
    @reactive.effect
    @reactive.event(input.aspect, input.base_slim, input.use_extra_terms)
    def _reset_slim_on_context_change():
        state = slim_state.get()
        if state["added"] or state["removed"]:
            ui.notification_show("Aspect/base slim changed - working slim reset.", type="warning")
        slim_state.set({"added": set(), "removed": set()})
        drill_term.set(None)

    @reactive.calc
    def base_slim_ids():
        # base slim ids for the current aspect, plus curated extras when requested
        aspect = input.aspect()
        base_choice = input.base_slim()
        if base_choice == "none":
            ids = set()
        elif base_choice.startswith(SAVED_BASE_PREFIX):
            ids = _saved_slim_term_ids(base_choice[len(SAVED_BASE_PREFIX):])
        else:
            ids = set(DATA["base_slims"].get((base_choice, aspect), frozenset()))
        if input.use_extra_terms() and aspect == "C":
            ids |= set(EXTRA_SLIM_TERMS["C"])
        return frozenset(ids)

    @reactive.calc
    def working_slim():
        # base slim, plus user-added terms, minus user-removed and globally excluded terms
        state = slim_state.get()
        ids = (set(base_slim_ids()) | state["added"]) - state["removed"]
        ids -= excluded_mod.excluded_ids_for_aspect(excluded_state.get(), input.aspect())
        return frozenset(ids)

    @reactive.calc
    def gene_symbols():
        # uploaded file wins over the dropdown selection
        upload = input.gene_upload()
        if upload:
            genes = load_genes(upload[0]["datapath"])
        else:
            genes = DATA["gene_lists"].get(input.gene_list(), [])
        return genes

    @reactive.calc
    def grouping():
        return compute_grouping(gene_symbols(), input.aspect(), working_slim(), DATA)

    # --- promote / remove / reset / load buttons mutate slim_state ---

    @reactive.effect
    @reactive.event(input.btn_promote)
    def _promote_terms():
        selected = uncovered_terms_grid.data_view(selected=True)
        if selected is None or selected.empty:
            ui.notification_show("Select one or more rows first.", type="warning")
            return
        state = slim_state.get()
        state = {"added": state["added"] | set(selected["go_id"]), "removed": state["removed"]}
        slim_state.set(state)

    @reactive.effect
    @reactive.event(input.btn_remove)
    def _remove_terms():
        selected = working_slim_grid.data_view(selected=True)
        if selected is None or selected.empty:
            ui.notification_show("Select one or more rows first.", type="warning")
            return
        state = slim_state.get()
        state = {"added": state["added"] - set(selected["go_id"]), "removed": state["removed"] | set(selected["go_id"])}
        slim_state.set(state)

    @reactive.effect
    @reactive.event(input.btn_reset)
    def _reset_button():
        slim_state.set({"added": set(), "removed": set()})
        drill_term.set(None)

    @reactive.effect
    def _track_drill_selection():
        selected = working_slim_grid.data_view(selected=True)
        if selected is not None and not selected.empty:
            drill_term.set(selected.iloc[0]["go_id"])

    @reactive.effect
    @reactive.event(input.btn_add_subterms)
    def _add_subterms():
        selected = subterms_grid.data_view(selected=True)
        if selected is None or selected.empty:
            ui.notification_show("Select one or more subterm rows first.", type="warning")
            return
        state = slim_state.get()
        new_ids = {_go_id_value(v) for v in selected["go_id"]}
        state = {"added": state["added"] | new_ids, "removed": state["removed"]}
        slim_state.set(state)

    @reactive.effect
    def _refresh_saved_slim_choices():
        choices = saved_slims.list_saved(config.SAVED_SLIMS_DIR)
        ui.update_select("load_slim_choice", choices=choices)

    @reactive.effect
    def _refresh_base_slim_choices():
        # re-run whenever the aspect changes or a new slim is saved, so saved
        # slims for the current aspect show up as full "base slim" choices
        slim_files_version.get()
        aspect = input.aspect()
        choices = dict(STATIC_BASE_SLIM_CHOICES)
        choices.update(saved_slims.list_saved_for_aspect(config.SAVED_SLIMS_DIR, aspect, key_prefix=SAVED_BASE_PREFIX))
        current = input.base_slim()
        selected = current if current in choices else None
        ui.update_select("base_slim", choices=choices, selected=selected)

    @reactive.effect
    @reactive.event(input.btn_load)
    def _load_saved_slim():
        slug = input.load_slim_choice()
        if not slug:
            return
        payload = saved_slims.read_slim(config.SAVED_SLIMS_DIR / f"{slug}.json")
        added = {t["id"] for t in payload.get("added_terms", [])}
        removed = {t["id"] for t in payload.get("removed_terms", [])}
        slim_state.set({"added": added, "removed": removed})
        ui.notification_show(f"Loaded '{payload.get('name', slug)}'.", type="message")

    @reactive.effect
    @reactive.event(input.btn_save)
    def _save_slim():
        name = input.save_name().strip()
        if not name:
            ui.notification_show("Enter a name before saving.", type="warning")
            return
        state = slim_state.get()
        g = grouping()
        names = DATA["names"]
        payload = {
            "schema_version": 1,
            "name": name,
            "aspect": input.aspect(),
            "base_slim": input.base_slim(),
            "include_extra_terms": input.use_extra_terms(),
            "use_part_of": input.use_part_of(),
            "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "notes": input.save_notes(),
            "source": {
                "obo": config.OBO_PATH.name,
                "gaf": config.GAF_PATH.name,
                "taxon": config.TAXON,
                "gene_list": input.gene_list(),
            },
            "added_terms": [{"id": t, "name": names.get(t)} for t in sorted(state["added"])],
            "removed_terms": [{"id": t, "name": names.get(t)} for t in sorted(state["removed"])],
            "terms": [{"id": t, "name": names.get(t)} for t in sorted(working_slim())],
            "stats": {
                "n_terms": len(working_slim()),
                "n_genes_covered": g["coverage"]["n_covered"],
                "n_unmapped": g["coverage"]["n_annotated_but_unmapped"],
                "n_no_annotation": g["coverage"]["n_no_annotation"],
                "largest_group": int(g["term_counts"]["n_genes"].max()) if not g["term_counts"].empty else 0,
            },
        }
        path = saved_slims.write_slim(config.SAVED_SLIMS_DIR, payload)
        slim_files_version.set(slim_files_version.get() + 1)
        ui.notification_show(f"Saved to {path.name}", type="message")

    # --- slim hygiene: exclude / un-exclude buttons ---

    @reactive.effect
    @reactive.event(input.btn_exclude)
    def _exclude_terms():
        selected = audit_grid.data_view(selected=True)
        if selected is None or selected.empty:
            ui.notification_show("Select one or more rows first.", type="warning")
            return
        state = excluded_state.get()
        aspect = input.aspect()
        for row in selected.itertuples():
            excluded_mod.add_exclude(state, aspect, row.go_id, row.name, "excluded via Slim hygiene tab")
        excluded_state.set(dict(state))
        excluded_mod.save_excluded(config.EXCLUDE_FILE, state)

    @reactive.effect
    @reactive.event(input.btn_unexclude)
    def _unexclude_terms():
        selected = exclude_grid.data_view(selected=True)
        if selected is None or selected.empty:
            ui.notification_show("Select one or more rows first.", type="warning")
            return
        state = excluded_state.get()
        aspect = input.aspect()
        for row in selected.itertuples():
            excluded_mod.remove_exclude(state, aspect, row.go_id)
        excluded_state.set(dict(state))
        excluded_mod.save_excluded(config.EXCLUDE_FILE, state)

    # --- coverage tab ---

    @render.text
    def slim_summary():
        state = slim_state.get()
        return f"{len(working_slim())} terms = {len(base_slim_ids())} base + {len(state['added'])} added - {len(state['removed'])} removed"

    @render.text
    def vb_input():
        return str(grouping()["coverage"]["n_input"])

    @render.text
    def vb_covered():
        return str(grouping()["coverage"]["n_covered"])

    @render.text
    def vb_unmapped():
        return str(grouping()["coverage"]["n_annotated_but_unmapped"])

    @render.text
    def vb_no_annotation():
        return str(grouping()["coverage"]["n_no_annotation"])

    @render.text
    def vb_unmatched():
        return str(grouping()["coverage"]["n_unmatched_symbol"])

    @render_widget
    def group_size_plot():
        counts = grouping()["partition_counts"]
        ordered = counts[counts["group_id"].notna()].sort_values("n_genes")
        melted = ordered.melt(
            id_vars=["group_name", "n_genes"],
            value_vars=["n_splittable", "n_leaf"],
            var_name="kind",
            value_name="n",
        )
        melted["n_splittable"] = melted["group_name"].map(ordered.set_index("group_name")["n_splittable"])
        kind_labels = {
            "n_splittable": "has a more specific term not in the slim",
            "n_leaf": "this term is their most specific annotation",
        }
        melted["kind"] = melted["kind"].map(kind_labels)
        fig = px.bar(
            melted,
            x="n",
            y="group_name",
            color="kind",
            orientation="h",
            color_discrete_map={
                kind_labels["n_splittable"]: "#93c5fd",
                kind_labels["n_leaf"]: "#1e3a8a",
            },
            labels={"n": "Genes", "group_name": "", "kind": ""},
            title="Partition group sizes (one gene per group)",
            custom_data=["group_name", "n_genes", "n_splittable"],
        )
        fig.update_traces(
            hovertemplate="%{customdata[0]}<br>%{customdata[1]} total genes<br>%{customdata[2]} genes can be further refined<extra></extra>"
        )
        fig.update_yaxes(categoryorder="array", categoryarray=ordered["group_name"].tolist())
        max_total = ordered["n_genes"].max()
        fig.update_xaxes(range=[0, max_total * 1.12])
        fig.update_layout(
            barmode="stack",
            height=min(600, max(300, 30 * len(ordered))),
            margin={"l": 0, "r": 10, "t": 40, "b": 10},
            legend=dict(orientation="v", yanchor="bottom", y=0.01, xanchor="right", x=0.99, bgcolor="rgba(255,255,255,0.7)"),
        )

        partition = grouping()["partition"]
        group_totals = ordered.set_index("group_name")["n_genes"]

        fig.add_scatter(
            x=group_totals.values,
            y=group_totals.index,
            mode="text",
            text=group_totals.values,
            textposition="middle right",
            showlegend=False,
            hoverinfo="skip",
        )

        def _add_tagged_marker(gene_set, name, color, hover_label):
            if partition.empty or not gene_set:
                return
            n = (
                partition[partition["gene"].isin(gene_set)]
                .groupby("group_name")["gene"]
                .nunique()
                .reindex(ordered["group_name"], fill_value=0)
            )
            pct = (n / group_totals.reindex(n.index) * 100).fillna(0)
            fig.add_scatter(
                x=n.values,
                y=n.index,
                mode="markers",
                marker=dict(symbol="diamond", size=10, color=color, line=dict(width=1, color="rgba(0,0,0,0.33)")),
                name=name,
                customdata=list(zip(n.values, pct.values)),
                hovertemplate=f"%{{y}}<br>%{{customdata[0]}} {hover_label}, %{{customdata[1]:.1f}}% of group<extra></extra>",
            )

        _add_tagged_marker(
            DATA["tagged_genes"], "WormTagDB fluorescent tags", "#f97316", "tagged genes (WormTagDB)"
        )
        _add_tagged_marker(
            DATA["tagged_genes_cgc"],
            "CGC-available fluorescent tags",
            "#22c55e",
            "tagged genes with a CGC strain",
        )
        return fig

    # --- slim groups tab ---

    @render.data_frame
    def term_counts_grid():
        return render.DataGrid(grouping()["term_counts"], selection_mode="row", filters=True, width="100%")

    @render.data_frame
    def partition_counts_grid():
        return render.DataGrid(grouping()["partition_counts"], selection_mode="row", filters=True, width="100%")

    @reactive.effect
    def _track_multi_selection():
        selected = term_counts_grid.data_view(selected=True)
        if selected is not None and not selected.empty:
            selected_group.set(("multi", selected.iloc[0]["slim_id"]))

    @reactive.effect
    def _track_partition_selection():
        selected = partition_counts_grid.data_view(selected=True)
        if selected is not None and not selected.empty:
            selected_group.set(("partition", selected.iloc[0]["group_id"]))

    @render.text
    def selected_group_label():
        sel = selected_group.get()
        return f"Selected: {sel[1]}" if sel else "Click a row in either table above."

    @render.data_frame
    def selected_group_genes_grid():
        sel = selected_group.get()
        if sel is None:
            return render.DataGrid(_empty_gene_frame(), width="100%")
        source, go_id = sel
        g = grouping()
        genes = genes_for_group(g["per_pair"], go_id) if source == "multi" else genes_for_partition_group(g["partition"], go_id)
        return render.DataGrid(pd.DataFrame({"gene": genes}), width="100%")

    @render.download_button(filename=lambda: f"group_{(selected_group.get() or (None, 'none'))[1]}.txt")
    def dl_group_genes():
        sel = selected_group.get()
        if sel is None:
            yield gene_list_to_text_bytes([])
            return
        source, go_id = sel
        g = grouping()
        genes = genes_for_group(g["per_pair"], go_id) if source == "multi" else genes_for_partition_group(g["partition"], go_id)
        yield gene_list_to_text_bytes(genes)

    # --- add terms tab ---

    @render.data_frame
    def uncovered_terms_grid():
        # uncovered terms are, by construction, never already in the working slim,
        # but excluded terms (e.g. aspect roots) shouldn't be offered as candidates either
        excluded_ids = excluded_mod.excluded_ids_for_aspect(excluded_state.get(), input.aspect())
        df = grouping()["uncovered_terms"]
        df = df[~df["go_id"].isin(excluded_ids)]
        return render.DataGrid(df, selection_mode="rows", filters=True, width="100%")

    @reactive.calc
    def working_slim_df():
        names = DATA["names"]
        depth = DATA["depth"]
        term_counts = grouping()["term_counts"]
        n_genes_by_term = term_counts.set_index("slim_id")["n_genes"].to_dict()
        rows = [
            {"go_id": t, "name": names.get(t), "n_genes": n_genes_by_term.get(t, 0), "depth": depth.get(t, 0)}
            for t in sorted(working_slim())
        ]
        cols = ["go_id", "name", "n_genes", "depth"]
        return pd.DataFrame(rows, columns=cols)

    @render.data_frame
    def working_slim_grid():
        return render.DataGrid(working_slim_df(), selection_mode="rows", filters=True, width="100%")

    @render.download_button(filename=lambda: f"working_slim_{input.aspect()}.xlsx")
    def dl_working_slim_xlsx():
        yield df_to_xlsx_bytes(working_slim_df(), sheet_name="working_slim")

    @render.text
    def drill_label():
        go_id = drill_term.get()
        if go_id is None:
            return "Select a row in 'Current working slim' above to see its split."
        aspect = input.aspect()
        breakdown, n_genes, n_leaf = compute_term_breakdown(
            go_id, grouping()["per_pair"], DATA["ann"][aspect], DATA["anc"], DATA["names"], DATA["depth"]
        )
        excluded_ids = excluded_mod.excluded_ids_for_aspect(excluded_state.get(), aspect)
        breakdown = breakdown[~breakdown["go_id"].isin(excluded_ids)]
        name = DATA["names"].get(go_id, go_id)
        if go_id not in working_slim():
            return f"{name} ({go_id}) is no longer in the working slim."
        return (
            f"Splitting {name} ({go_id}) - {n_genes} genes currently assigned here, "
            f"{len(breakdown)} more specific raw terms found among them below. "
            f"{n_leaf} of those genes have {go_id} as their most specific term and can't be split any further."
        )

    @render.data_frame
    def subterms_grid():
        go_id = drill_term.get()
        empty_cols = ["go_id", "go_name", "n_genes", "depth"]
        if go_id is None or go_id not in working_slim():
            return render.DataGrid(pd.DataFrame(columns=empty_cols), selection_mode="rows", filters=True, width="100%")
        aspect = input.aspect()
        breakdown, _n_genes, _n_leaf = compute_term_breakdown(
            go_id, grouping()["per_pair"], DATA["ann"][aspect], DATA["anc"], DATA["names"], DATA["depth"]
        )
        excluded_ids = excluded_mod.excluded_ids_for_aspect(excluded_state.get(), aspect)
        breakdown = breakdown[~breakdown["go_id"].isin(excluded_ids)].copy()
        breakdown["go_id"] = [_go_id_cell(t, working_slim()) for t in breakdown["go_id"]]
        return render.DataGrid(breakdown, selection_mode="rows", filters=True, width="100%")

    # --- gene detail tab ---

    @render.data_frame
    def gene_detail_grid():
        aspect = input.aspect()
        ann = DATA["ann"][aspect]
        names = DATA["names"]
        needle = (input.gene_search() or "").casefold()
        genes = [g for g in gene_symbols() if needle in g.casefold()] if needle else gene_symbols()

        per_pair = grouping()["per_pair"]
        partition = grouping()["partition"]
        slim_terms_by_gene = per_pair.groupby("gene")["slim_id"].apply(list).to_dict()
        group_by_gene = partition.set_index("gene")["group_name"].to_dict() if not partition.empty else {}

        rows = []
        for g in genes:
            raw = sorted(ann.get(g, ()))
            rows.append({
                "gene": g,
                "raw_go_terms": ";".join(raw),
                "raw_go_names": ";".join(names.get(t, "") for t in raw),
                "slim_terms": ";".join(slim_terms_by_gene.get(g, [])),
                "assigned_group": group_by_gene.get(g),
            })
        cols = ["gene", "raw_go_terms", "raw_go_names", "slim_terms", "assigned_group"]
        return render.DataGrid(pd.DataFrame(rows, columns=cols), filters=True, width="100%")

    # --- slim hygiene tab ---

    @render.data_frame
    def audit_grid():
        term_counts = grouping()["term_counts"]
        term_gene_counts = term_counts.set_index("slim_id")["n_genes"].to_dict()
        df = excluded_mod.audit_slim(
            base_slim_ids(), DATA["godag"], DATA["replaced_by"], DATA["names"],
            term_gene_counts, grouping()["coverage"]["n_input"], DATA["depth"],
        )
        return render.DataGrid(df, selection_mode="rows", filters=True, width="100%")

    @render.data_frame
    def exclude_grid():
        state = excluded_state.get()
        aspect = input.aspect()
        rows = list(state.get("global", [])) + list(state.get(aspect, []))
        df = pd.DataFrame(rows, columns=["id", "name", "reason"]).rename(columns={"id": "go_id"})
        return render.DataGrid(df, selection_mode="rows", filters=True, width="100%")

    # --- downloads ---

    @render.download_button(filename=lambda: f"gene_groups_{input.aspect()}.tsv")
    def dl_full_table():
        genes = gene_symbols()
        g = grouping()
        table = build_full_gene_table(genes, g, DATA["ann"], input.aspect())
        yield df_to_tsv_bytes(table)


def _empty_gene_frame():
    return pd.DataFrame({"gene": []})


def _saved_slim_term_ids(slug):
    """Resolve a 'saved:<slug>' base slim choice to its full term id set."""
    path = config.SAVED_SLIMS_DIR / f"{slug}.json"
    if not path.exists():
        return set()
    payload = saved_slims.read_slim(path)
    return {t["id"] for t in payload.get("terms", [])}

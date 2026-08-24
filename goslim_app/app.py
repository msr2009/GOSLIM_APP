"""Entry point: run with `shiny run --reload goslim_app/app.py`."""

from server_main import server
from shiny import App
from ui_layout import app_ui

app = App(app_ui, server)

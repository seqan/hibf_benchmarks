"""Creates a plot."""

import csv
import os

from bokeh.layouts import column, row
from bokeh.models import Div, TabPanel, Tabs
from bokeh.palettes import Set2_4, Set2_6
from bokeh.plotting import curdoc, figure, output_file, save
from bokeh.themes import Theme

from components.convert_data import convert_size_data, convert_time_data
from components.helpers import convert_list_to_string, get_max_result
from components.log_init import log_init
from components.plot_css_html import create_latex_text, create_vercel_div, get_global_style, get_tab_style
from components.plot_style import add_legend, add_second_y_axis, configure_size_plot, configure_time_plot

log_init(snakemake.log[0])

TIME_FORMAT = snakemake.config["PLOT_TIME_FORMAT"]
TIME_NAMES = snakemake.config["TIME_NAMES"]
SIZE_FORMAT = snakemake.config["PLOT_SIZE_FORMAT"]
SIZE_NAMES = snakemake.config["SIZE_NAMES"]
KEYS = snakemake.config["KEYS"]
KEYS_NAMES = snakemake.config["KEYS_NAMES"]


def create_time_plot(time_data, y_range, max_result_time, file_name, scale_in_minutes):
    """Creates the time plot."""
    plot = figure(
        y_range=y_range,
        x_range=(max_result_time, 0),
        tools="",
    )
    renderers = plot.hbar_stack(
        stackers=TIME_FORMAT[1:], y=TIME_FORMAT[0], height=0.4, source=(time_data), color=Set2_6
    )
    add_legend(plot, renderers, file_name, TIME_NAMES, SIZE_NAMES, "TIME_FORMAT", "left")
    configure_time_plot(plot, scale_in_minutes)
    add_second_y_axis(plot, y_range)
    return plot


def create_size_plot(size_data, y_range, max_result_size, file_name):
    """Creates the size plot."""
    plot = figure(
        y_range=y_range,
        x_range=(0, max_result_size),
        x_axis_label="size in GB",
        tools="",
    )
    renderers = plot.hbar_stack(SIZE_FORMAT[1:], y=SIZE_FORMAT[0], height=0.4, source=(size_data), color=Set2_4)
    configure_size_plot(plot)
    add_legend(plot, renderers, file_name, TIME_NAMES, SIZE_NAMES, "SIZE_FORMAT", "right")
    return plot


def create_plot():
    """Creates the final plot."""
    output_file(filename="results/html/index.html", title="HIBF Benchmarks")
    curdoc().theme = Theme(filename="workflow/scripts/bokeh_plot/plot_theme.yaml")
    tabs = []
    for file_name_index, file_name in enumerate(KEYS):
        with open(os.path.join("results", "prepared_time", file_name), "r", encoding="utf-8") as timing_file, open(
            os.path.join("results", "prepared_size", file_name), "r", encoding="utf-8"
        ) as size_file:
            time_reader = csv.reader(timing_file, delimiter="\t")
            size_reader = csv.reader(size_file, delimiter="\t")
            time_data_list, size_data_list = list(time_reader), list(size_reader)
            if len(time_data_list) == 0 or len(size_data_list) == 0:
                continue
            time_data = convert_time_data(time_data_list, file_name, TIME_FORMAT)
            size_data = convert_size_data(size_data_list, file_name, SIZE_FORMAT)
            max_result_time = get_max_result(time_data_list[1:], 1.01)
            max_result_size = get_max_result(size_data_list[1:5], 1.01)
            scale_in_minutes = max_result_time > 120

            size_y_range = convert_list_to_string(size_data["SUBKEY"])
            time_y_range = convert_list_to_string(time_data["SUBKEY"])
            y_range = size_y_range if len(size_y_range) > len(time_y_range) else time_y_range

            plot1 = create_time_plot(time_data, y_range, max_result_time, file_name, scale_in_minutes)
            plot2 = create_size_plot(size_data, plot1.y_range, max_result_size, file_name)
            both_plots = row(plot1, plot2, sizing_mode="scale_both")

            vercel_div = create_vercel_div()
            all_elements = column(both_plots, vercel_div, sizing_mode="scale_both")
            tabs.append(TabPanel(child=all_elements, title=KEYS_NAMES[file_name_index]))
    latex_text = create_latex_text()
    tabs.append(
        TabPanel(child=Div(text=latex_text, styles={"color": "white", "font-size": "14px"}), title="Description")
    )
    tab_style = get_tab_style()
    global_style = get_global_style()
    save(Tabs(tabs=tabs, sizing_mode="scale_both", stylesheets=[tab_style, global_style]))


create_plot()
print("Plot created successfully.")

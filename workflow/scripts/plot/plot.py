"""Creates a plot."""

import pandas as pd
from bokeh.layouts import row
from bokeh.models import TabPanel
from bokeh.palettes import Set2_4, Set2_6
from bokeh.plotting import curdoc, figure, output_file
from bokeh.themes import Theme

from components.convert_data import prepare_size_data, prepare_time_data
from components.log_init import log_init
from components.plot_style import add_legend, add_second_y_axis, configure_size_plot, configure_time_plot, save_tabs

log_init(snakemake.log[0]) # type: ignore

SIZE_INPUT = snakemake.input["SIZE_INPUT"] # type: ignore
TIME_INPUT = snakemake.input["TIME_INPUT"] # type: ignore

PLOT_FILE = snakemake.output["PLOT_FILE"] # type: ignore

THEME = snakemake.params["THEME"] # type: ignore
KEYS = snakemake.params["KEYS"] # type: ignore
TIME = snakemake.params["TIME"] # type: ignore
SIZE = snakemake.params["SIZE"] # type: ignore

TIME_NAMES = [TIME["NAMES"].get(key, key) for key in TIME["FORMAT"]]
SIZE_NAMES = [SIZE["NAMES"].get(key, key) for key in SIZE["FORMAT"]]


def create_time_plot(data, y_range, x_range, file_name):
    """Creates the time plot."""
    plot = figure(
        y_range=y_range,
        x_range=(x_range, 0),
        toolbar_location="left",
        tools="",
    )
    renderers = plot.hbar_stack(stackers=TIME["FORMAT"], y=("SUBKEY"), height=0.4, source=(data), color=Set2_6)
    add_legend(plot, renderers, file_name, TIME_NAMES, SIZE_NAMES, "TIME_FORMAT", "left")
    configure_time_plot(plot, x_range > 120)
    add_second_y_axis(plot, y_range)
    return plot


def create_size_plot(size_data, y_range, max_result_size, file_name):
    """Creates the size plot."""
    plot = figure(
        y_range=y_range,
        x_range=(0, max_result_size),
        x_axis_label="size in GB",
        toolbar_location="right",
        tools="",
    )
    renderers = plot.hbar_stack(stackers=SIZE["FORMAT"], y=("SUBKEY"), height=0.4, source=(size_data), color=Set2_4)
    configure_size_plot(plot)
    add_legend(plot, renderers, file_name, TIME_NAMES, SIZE_NAMES, "SIZE_FORMAT", "right")
    return plot


def create_plot():
    """Creates the final plot."""
    output_file(filename=PLOT_FILE, title="HIBF Benchmarks")
    curdoc().theme = Theme(filename=THEME)
    tabs = []
    with open(SIZE_INPUT, "r", encoding="utf-8") as size_file, open(TIME_INPUT, "r", encoding="utf-8") as timing_file:
        time_data = pd.read_csv(timing_file, delimiter="\t")
        size_data = pd.read_csv(size_file, delimiter="\t")

    for key in KEYS.keys():
        time_dic = prepare_time_data(time_data, (key, KEYS[key]), TIME)
        size_dic = prepare_size_data(size_data, (key, KEYS[key]), SIZE)

        time_x_range = round(max(time_dic["TOTAL_TIME"]) * 1.03, 3)
        size_x_range = round(max(size_dic["GB_TOTAL_SIZE"]) * 1.03, 3)

        size_y_range = [str(i) for i in size_dic["SUBKEY"]]
        time_y_range = [str(i) for i in time_dic["SUBKEY"]]
        y_range = size_y_range if len(size_y_range) > len(time_y_range) else time_y_range

        plot1 = create_time_plot(time_dic, y_range, time_x_range, key)
        # plot1.y_range causes panning/zooming to be synchronized between the two plots
        plot2 = create_size_plot(size_dic, plot1.y_range, size_x_range, key)
        both_plots = row(plot1, plot2, sizing_mode="scale_both")

        tabs.append(TabPanel(child=both_plots, title=KEYS[key]))

    save_tabs(tabs)


create_plot()
print("Plot created successfully.")

#!/usr/bin/env python
"""Creates a plot."""

import csv
import os

from bokeh.io import output_notebook
from bokeh.layouts import row
from bokeh.models import AdaptiveTicker, CustomJSTickFormatter, FactorRange, HoverTool, LinearAxis, TabPanel, Tabs
from bokeh.palettes import Set2_4, Set2_7
from bokeh.plotting import figure, output_file, save, show

from shared import html_file, output_sizes, output_timings

time_structure = [
    "subkeys",
    "query_file_io_in_seconds",
    "determine_query_length_in_seconds",
    "compute_minimiser_avg_per_thread_in_seconds",
    "generate_results_avg_per_thread_in_seconds",
    "load_index_in_seconds",
    "query_ibf_avg_per_thread_in_seconds",
    "wall_clock_time_in_seconds",
]
size_structure = [
    "subkeys",
    "size_level_0",
    "size_level_1",
    "size_level_2",
    "size_level_3",
]


def convert_list_to_float(data):
    """Returns a list containing the values of the given list as floats."""
    return [float(i) for i in data]


def convert_list_to_string(data):
    """Returns a list containing the values of the given list as strings."""
    return [str(i) for i in data]


def convert_time_data(data, key):
    """Returns a dictionary containing the given data in the format for the time plot."""
    export = {}
    export[time_structure[0]] = [f"{key} = {value}" for value in data[0]]
    export["value"] = data[0]
    export["all_times"] = [sum(float(i) for i in sublist) for sublist in zip(*data[1:])]
    for i, element in enumerate(data[1:]):
        export[time_structure[i + 1]] = convert_list_to_float(element)
        export[f"{time_structure[i+1]}_percentage"] = [
            round((float(i) / export["all_times"][0]) * 100, 2) for i in export[time_structure[i + 1]]
        ]
    return export


def convert_size_data(data, key):
    """Returns a dictionary containing the given data in the format for the size plot."""
    export = {}
    export[size_structure[0]] = [f"{key} = {value}" for value in data[0]]
    export["value"] = data[0]
    export["sizes"] = [sum(float(i) for i in sublist) for sublist in zip(*data[1:])]
    for i, element in enumerate(data[1:]):
        export[size_structure[i + 1]] = convert_list_to_float(element)
        export[f"{size_structure[i+1]}_percentage"] = [
            round((float(i) / export["sizes"][0]) * 100, 2) for i in export[size_structure[i + 1]]
        ]
    return export


def add_arrays(data):
    """Returns a list containing the sum of the given lists."""
    return [sum(float(i) for i in sublist) for sublist in zip(*data)]


def get_max_result(data, factor):
    """Returns the maximum value of the given list multiplied by the given factor."""
    return round(max(add_arrays(data)) * factor)


def hex_to_rgb(hex_code):
    """Returns the rgb values of the given hex code."""
    hex_code = hex_code.lstrip("#")
    return tuple(int(hex_code[i : i + 2], 16) for i in (0, 2, 4))


def mix_with_white(color, alpha):
    """Simulates applying the given alpha value to the given color."""
    white = (255, 255, 255)
    r = (1 - alpha) * color[0] + alpha * white[0]
    g = (1 - alpha) * color[1] + alpha * white[1]
    b = (1 - alpha) * color[2] + alpha * white[2]
    return (int(r), int(g), int(b))


def create_plot(interactive=False):
    """Creates the plot."""
    output_file(filename=html_file, title="Static HTML file")
    if interactive:
        output_notebook()
    tabs = []
    files_names = ["alpha", "hash", "k", "r-fpr", "none", "U", "U+R"]
    special_files_names = ["none", "U", "U+R"]
    files_names_titles = ["alpha", "hash", "k-mer", "r-fpr", "no U+R", "U", "U+R"]
    for i, file_name in enumerate(files_names):
        with open(os.path.join(output_timings, file_name + ".csv"), "r", encoding="utf-8") as timing_file, open(
            os.path.join(output_sizes, file_name + ".csv"), "r", encoding="utf-8"
        ) as size_file:
            if file_name in special_files_names:
                file_name = "t_max"
            time_reader, size_reader = csv.reader(timing_file), csv.reader(size_file)
            time_data_list, size_data_list = list(time_reader), list(size_reader)
            time_data, size_data = convert_time_data(time_data_list, file_name), convert_size_data(
                size_data_list[:5], file_name
            )
            max_result_time, max_result_size = get_max_result(time_data_list[1:], 1.8), get_max_result(
                size_data_list[1:], 1.8
            )
            p1 = figure(
                y_range=convert_list_to_string(size_data[time_structure[0]]),
                height=500,
                width=800,
                x_range=(max_result_time, 0),
                x_axis_label="time in minutes",
                toolbar_location="left",
                tools="",
            )
            renderers1 = p1.hbar_stack(
                stackers=time_structure[1:],
                y=time_structure[0],
                height=0.4,
                source=(time_data),
                color=Set2_7,
                legend_label=time_structure[1:],
            )
            for r in renderers1:
                key, tag = r.name, r.name
                if key == "wall_clock_time_in_seconds":
                    tag = "time_left"
                p1.add_tools(
                    HoverTool(
                        tooltips=[
                            (f"{file_name}", "@value"),
                            ("wall_clock_time_in_seconds", "@all_times{0.00} sek"),
                            (tag, "@$name{0.00} sek"),
                            ("Percentage", f"@{key}_percentage{{0.00}}%"),
                        ],
                        renderers=[r],
                    )
                )
            p1.toolbar_location = None
            p1.toolbar.logo = None
            p1.yaxis.visible = False
            p1.extra_y_ranges = {
                "zusätzliche_achse": FactorRange(factors=convert_list_to_string(size_data[time_structure[0]]))
            }
            zweite_y_achse = LinearAxis(y_range_name="zusätzliche_achse")
            zweite_y_achse.major_label_text_font_size = "1pt"
            zweite_y_achse.major_label_text_color = "#ffffff"
            p1.add_layout(zweite_y_achse, "right")
            p1.xaxis.ticker = AdaptiveTicker(base=60, mantissas=[1, 2, 5], min_interval=60, max_interval=600)
            p1.xaxis.formatter = CustomJSTickFormatter(
                code="""
                return (tick / 60);
            """
            )
            p1.y_range.range_padding = 0.1
            p1.ygrid.grid_line_color = None
            p1.axis.minor_tick_line_color = None
            p1.legend.location = "top_left"
            p1.legend.label_text_font_size = "10pt"
            p1.legend.glyph_height = 12
            p1.legend.glyph_width = 12
            p1.legend.border_line_width = 2
            p1.legend.title = "Standard parameters:\n\t\t\tt_max = 192\n\t\tunion estimation (U) = yes\n\t\t\trearrangement (R) = yes\n\t\t\tk-mer size (k) = 32\n\t\t\tnumber of hash function (hash) = 2\n\t\t\talpha = 1.2\n\t\t\trelaxed false positive rate (r-fpr) = 0.5\n\t\t\tmaximum false positive rate (fpr) = 0.05"
            p1.legend.title_text_font_size = "10pt"
            p1.yaxis.major_tick_line_color = None
            p1.yaxis.minor_tick_line_color = None
            p1.outline_line_color = None
            p2 = figure(
                y_range=convert_list_to_string(size_data[size_structure[0]]),
                height=500,
                width=400,
                x_range=(0, max_result_size),
                x_axis_label="size in GB",
                toolbar_location="right",
                tools="",
            )
            renderers2 = p2.hbar_stack(
                size_structure[1:],
                y=size_structure[0],
                height=0.4,
                source=(size_data),
                color=Set2_4,
                legend_label=size_structure[1:],
            )
            for r in renderers2:
                key = r.name
                p2.add_tools(
                    HoverTool(
                        tooltips=[
                            (f"{file_name}", "@value"),
                            ("all_level_size", "@sizes{0.00}GB"),
                            (key, "@$name{0.00}GB"),
                            ("Percentage", f"@{key}_percentage{{0.00}}%"),
                        ],
                        renderers=[r],
                    )
                )
            p2.toolbar.logo = None
            p2.toolbar_location = None
            p2.y_range.range_padding = 0.1
            p2.ygrid.grid_line_color = None
            p2.axis.minor_tick_line_color = None
            p2.outline_line_color = None
            p2.legend.location = "top_right"
            p2.yaxis.major_tick_line_color = None
            p2.yaxis.minor_tick_line_color = None
            p2.yaxis.major_label_standoff = 15
            both_plots = row(p1, p2)
            tabs.append(TabPanel(child=both_plots, title=files_names_titles[i]))
    if interactive:
        show(Tabs(tabs=tabs))
    else:
        save(Tabs(tabs=tabs))


if __name__ == "__main__":
    create_plot()

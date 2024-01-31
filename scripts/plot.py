#!/usr/bin/env python
"""Creates a plot."""

import csv
import os

from bokeh import events
from bokeh.io import output_notebook
from bokeh.layouts import column, row
from bokeh.models import (
    AdaptiveTicker,
    CustomJS,
    CustomJSTickFormatter,
    Div,
    FactorRange,
    HoverTool,
    Legend,
    LinearAxis,
    TabPanel,
    Tabs,
)
from bokeh.palettes import Set2_4, Set2_6
from bokeh.plotting import curdoc, figure, output_file, save, show

from shared import html_dir, html_file, output_sizes, output_timings

time_structure = [
    "subkeys",
    "query_file_io_in_seconds",
    "determine_query_length_in_seconds",
    "compute_minimiser_avg_per_thread_in_seconds",
    "generate_results_avg_per_thread_in_seconds",
    "load_index_in_seconds",
    "query_ibf_avg_per_thread_in_seconds",
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
        export[f"{time_structure[i+1]}_percentage"] = devide_arrays_in_percentage(
            export[time_structure[i + 1]], export["all_times"]
        )
    return export


def convert_size_data(data, key):
    """Returns a dictionary containing the given data in the format for the size plot."""
    export = {}
    export[size_structure[0]] = [f"{key} = {value}" for value in data[0]]
    export["value"] = data[0]
    export["sizes"] = [sum(float(i) for i in sublist) for sublist in zip(*data[1:])]
    for i, element in enumerate(data[1:5]):
        export[size_structure[i + 1]] = convert_list_to_float(element)
        export[f"{size_structure[i+1]}_percentage"] = devide_arrays_in_percentage(
            export[size_structure[i + 1]], export["sizes"]
        )
    for i, element in enumerate(data[5:]):
        export[f"{size_structure[i+1]}_avg_load_factor"] = element
    return export


def add_arrays(data):
    """Returns a list containing the sum of the given lists."""
    return [sum(float(i) for i in sublist) for sublist in zip(*data)]


def get_max_result(data, factor):
    """Returns the maximum value of the given list multiplied by the given factor."""
    return round(max(add_arrays(data)) * factor)


def devide_arrays_in_percentage(list1, list2):
    """Returns the percentage of each element from list1 to list2."""
    return [round(float(i) / float(j) * 100, 2) for i, j in zip(list1, list2)]


def create_plot(interactive=False):
    """Creates the plot."""
    try:
        os.makedirs(html_dir, exist_ok=True)
        print(f"Directory {html_dir} created")
    except OSError:
        print("files replaced in directory {html_dir}")
    output_file(filename=html_file, title="HIBF Benchmarks")
    curdoc().theme = "dark_minimal"
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
                size_data_list, file_name
            )
            max_result_time, max_result_size = get_max_result(time_data_list[1:], 1.01), get_max_result(
                size_data_list[1:], 1.01
            )
            plot1 = figure(
                y_range=convert_list_to_string(size_data[time_structure[0]]),
                x_range=(max_result_time, 0),
                x_axis_label="time in minutes",
                toolbar_location="left",
                tools="",
            )
            renderers1 = plot1.hbar_stack(
                stackers=time_structure[1:], y=time_structure[0], height=0.4, source=(time_data), color=Set2_6
            )
            legend_items = []
            for renderer in renderers1:
                key, tag = renderer.name, renderer.name
                legend_items.append((key, [renderer]))
                if key == "wall_clock_time_in_seconds":
                    tag = "time_left"
                plot1.add_tools(
                    HoverTool(
                        tooltips=[
                            (f"{file_name}", "@value"),
                            ("wall_clock_time_in_seconds", "@all_times{0.00} sek"),
                            (tag, "@$name{0.00} sek"),
                            ("Percentage", f"@{key}_percentage{{0.00}}%"),
                        ],
                        renderers=[renderer],
                    )
                )
            plot1.toolbar_location = None
            plot1.toolbar.logo = None
            plot1.yaxis.visible = False
            plot1.extra_y_ranges = {
                "zusätzliche_achse": FactorRange(factors=convert_list_to_string(size_data[time_structure[0]]))
            }
            zweite_y_achse = LinearAxis(y_range_name="zusätzliche_achse")
            zweite_y_achse.major_label_text_font_size = "1pt"
            zweite_y_achse.major_label_text_color = "#15191c"
            plot1.add_layout(zweite_y_achse, "right")
            plot1.xaxis.ticker = AdaptiveTicker(base=60, mantissas=[1, 2, 5], min_interval=60, max_interval=600)
            plot1.xaxis.formatter = CustomJSTickFormatter(
                code="""
                return (tick / 60);
            """
            )
            plot1.y_range.range_padding = 0.1
            plot1.ygrid.grid_line_color = None
            plot1.axis.minor_tick_line_color = None

            plot1.yaxis.major_tick_line_color = None
            plot1.yaxis.minor_tick_line_color = None
            plot1.outline_line_color = None
            plot1.sizing_mode = "scale_both"
            plot1.title.text = "Double click on legend/plot to hide/show the legend"
            plot1.title.align = "right"

            legend = Legend(items=legend_items)
            legend.location = "center"
            legend.orientation = "vertical"
            legend.glyph_height = 12
            legend.glyph_width = 12
            legend.click_policy = "mute"
            # legend.nrows = 3
            legend.title = "Default parameters:\n\t\t\tt_max = 192\n\t\tunion estimation (U) = yes\n\t\t\trearrangement (R) = yes\n\t\t\tk-mer size (k) = 32\n\t\t\tnumber of hash function (hash) = 2\n\t\t\talpha = 1.2\n\t\t\trelaxed false positive rate (r-fpr) = 0.5\n\t\t\tmaximum false positive rate (fpr) = 0.05"
            legend.title_text_color = "#e0e0e0"
            plot1.add_layout(legend, "left")
            toggle_legend_js = CustomJS(
                args={"legend": legend},
                code="""
                    legend.visible = !legend.visible
            """,
            )
            plot1.js_on_event(events.DoubleTap, toggle_legend_js)

            plot2 = figure(
                y_range=convert_list_to_string(size_data[size_structure[0]]),
                x_range=(0, max_result_size),
                x_axis_label="size in GB",
                toolbar_location="right",
                tools="",
            )
            legend_items = []
            renderers2 = plot2.hbar_stack(
                size_structure[1:], y=size_structure[0], height=0.4, source=(size_data), color=Set2_4
            )
            for renderer in renderers2:
                key = renderer.name
                legend_items.append((key, [renderer]))
                plot2.add_tools(
                    HoverTool(
                        tooltips=[
                            (f"{file_name}", "@value"),
                            ("all_level_size", "@sizes{0.00}GB"),
                            (key, "@$name{0.00}GB"),
                            ("Percentage", f"@{key}_percentage{{0.00}}%"),
                            ("avg_load_factor", f"@{key}_avg_load_factor{{0.00}}"),
                        ],
                        renderers=[renderer],
                    )
                )
            plot2.toolbar.logo = None
            plot2.toolbar_location = None
            plot2.y_range.range_padding = 0.1
            plot2.ygrid.grid_line_color = None
            plot2.axis.minor_tick_line_color = None
            plot2.outline_line_color = None
            plot2.title.text = "Click on legend entries to mute the corresponding bars"

            legend = Legend(items=legend_items)
            legend.location = "center"
            legend.glyph_height = 12
            legend.glyph_width = 12
            legend.click_policy = "mute"
            plot2.add_layout(legend, "right")
            toggle_legend_js = CustomJS(
                args={"legend": legend},
                code="""
                    legend.visible = !legend.visible
            """,
            )
            plot2.js_on_event(events.DoubleTap, toggle_legend_js)

            plot2.yaxis.major_tick_line_color = None
            plot2.yaxis.minor_tick_line_color = None
            plot2.yaxis.major_label_standoff = 15
            plot2.sizing_mode = "scale_both"
            both_plots = row(plot1, plot2)
            both_plots.sizing_mode = "scale_both"
            vercel_logo = """
                <svg width="209" height="40" viewBox="0 0 209 40" fill="none" xmlns="http://www.w3.org/2000/svg">
                <path d="M0 5C0 2.23858 2.23858 0 5 0H204C206.761 0 209 2.23858 209 5V35C209 37.7614 206.761 40 204 40H5C2.23858 40 0 37.7614 0 35V5Z" fill="black"/>
                <path fill-rule="evenodd" clip-rule="evenodd" d="M20 13L28 27H12L20 13Z" fill="white"/>
                <line x1="40.5" y1="2.18556e-08" x2="40.5" y2="40" stroke="#333333"/>
                <path d="M53.2784 26H55.0341V21.9091H57.4205C60.1193 21.9091 61.4545 20.2784 61.4545 18.1307C61.4545 15.9886 60.1307 14.3636 57.4261 14.3636H53.2784V26ZM55.0341 20.4205V15.8693H57.2386C58.9773 15.8693 59.6875 16.8125 59.6875 18.1307C59.6875 19.4489 58.9773 20.4205 57.2614 20.4205H55.0341ZM66.9432 26.1761C69.4034 26.1761 71.0114 24.375 71.0114 21.6761C71.0114 18.9602 69.4034 17.1591 66.9432 17.1591C64.483 17.1591 62.875 18.9602 62.875 21.6761C62.875 24.375 64.483 26.1761 66.9432 26.1761ZM66.9489 24.75C65.3409 24.75 64.5909 23.3466 64.5909 21.6705C64.5909 20 65.3409 18.5795 66.9489 18.5795C68.5455 18.5795 69.2955 20 69.2955 21.6705C69.2955 23.3466 68.5455 24.75 66.9489 24.75ZM74.5341 26H76.2614L78.0341 19.6989H78.1648L79.9375 26H81.6705L84.233 17.2727H82.4773L80.7784 23.6534H80.6932L78.9886 17.2727H77.233L75.517 23.6818H75.4318L73.7216 17.2727H71.9659L74.5341 26ZM89.3409 26.1761C91.2443 26.1761 92.5909 25.2386 92.9773 23.8182L91.3693 23.5284C91.0625 24.3523 90.3239 24.7727 89.358 24.7727C87.9034 24.7727 86.9261 23.8295 86.8807 22.1477H93.0852V21.5455C93.0852 18.392 91.1989 17.1591 89.2216 17.1591C86.7898 17.1591 85.1875 19.0114 85.1875 21.6932C85.1875 24.4034 86.767 26.1761 89.3409 26.1761ZM86.8864 20.875C86.9545 19.6364 87.8523 18.5625 89.233 18.5625C90.5511 18.5625 91.4148 19.5398 91.4205 20.875H86.8864ZM94.9702 26H96.669V20.6705C96.669 19.5284 97.5497 18.7045 98.7543 18.7045C99.1065 18.7045 99.5043 18.767 99.6406 18.8068V17.1818C99.4702 17.1591 99.1349 17.142 98.919 17.142C97.8963 17.142 97.0213 17.7216 96.7031 18.6591H96.6122V17.2727H94.9702V26ZM104.56 26.1761C106.463 26.1761 107.81 25.2386 108.196 23.8182L106.588 23.5284C106.281 24.3523 105.543 24.7727 104.577 24.7727C103.122 24.7727 102.145 23.8295 102.099 22.1477H108.304V21.5455C108.304 18.392 106.418 17.1591 104.44 17.1591C102.009 17.1591 100.406 19.0114 100.406 21.6932C100.406 24.4034 101.986 26.1761 104.56 26.1761ZM102.105 20.875C102.173 19.6364 103.071 18.5625 104.452 18.5625C105.77 18.5625 106.634 19.5398 106.639 20.875H102.105ZM113.456 26.1705C115.047 26.1705 115.672 25.1989 115.979 24.642H116.121V26H117.78V14.3636H116.081V18.6875H115.979C115.672 18.1477 115.092 17.1591 113.467 17.1591C111.359 17.1591 109.808 18.8239 109.808 21.6534C109.808 24.4773 111.337 26.1705 113.456 26.1705ZM113.831 24.7216C112.314 24.7216 111.524 23.3864 111.524 21.6364C111.524 19.9034 112.297 18.6023 113.831 18.6023C115.314 18.6023 116.109 19.8125 116.109 21.6364C116.109 23.4716 115.297 24.7216 113.831 24.7216ZM124.575 26H126.234V24.642H126.376C126.683 25.1989 127.308 26.1705 128.899 26.1705C131.013 26.1705 132.547 24.4773 132.547 21.6534C132.547 18.8239 130.99 17.1591 128.882 17.1591C127.263 17.1591 126.678 18.1477 126.376 18.6875H126.274V14.3636H124.575V26ZM126.24 21.6364C126.24 19.8125 127.036 18.6023 128.518 18.6023C130.058 18.6023 130.831 19.9034 130.831 21.6364C130.831 23.3864 130.036 24.7216 128.518 24.7216C127.058 24.7216 126.24 23.4716 126.24 21.6364ZM135.216 29.25C136.619 29.25 137.511 28.517 138.011 27.1648L141.619 17.2898L139.784 17.2727L137.574 24.0455H137.483L135.273 17.2727H133.455L136.648 26.1136L136.438 26.6932C136.006 27.8239 135.398 27.9261 134.466 27.6705L134.057 29.0625C134.261 29.1591 134.705 29.25 135.216 29.25ZM149.426 14.3636H146.693L150.71 26H153.881L157.892 14.3636H155.165L152.347 23.2045H152.239L149.426 14.3636ZM162.224 26.1705C164.384 26.1705 165.838 25.1193 166.179 23.5L163.94 23.3523C163.696 24.017 163.071 24.3636 162.264 24.3636C161.054 24.3636 160.287 23.5625 160.287 22.2614V22.2557H166.23V21.5909C166.23 18.625 164.435 17.1591 162.128 17.1591C159.56 17.1591 157.895 18.983 157.895 21.6761C157.895 24.4432 159.537 26.1705 162.224 26.1705ZM160.287 20.7557C160.338 19.7614 161.094 18.9659 162.168 18.9659C163.219 18.9659 163.946 19.7159 163.952 20.7557H160.287ZM167.81 26H170.23V21.0625C170.23 19.9886 171.014 19.25 172.082 19.25C172.418 19.25 172.878 19.3068 173.105 19.3807V17.233C172.889 17.1818 172.588 17.1477 172.344 17.1477C171.366 17.1477 170.565 17.7159 170.247 18.7955H170.156V17.2727H167.81V26ZM177.893 26.1705C180.217 26.1705 181.678 24.8068 181.791 22.8011H179.507C179.365 23.733 178.751 24.2557 177.922 24.2557C176.791 24.2557 176.058 23.3068 176.058 21.6364C176.058 19.9886 176.797 19.0455 177.922 19.0455C178.808 19.0455 179.376 19.6307 179.507 20.5H181.791C181.689 18.483 180.161 17.1591 177.882 17.1591C175.234 17.1591 173.598 18.9943 173.598 21.6705C173.598 24.3239 175.206 26.1705 177.893 26.1705ZM187.318 26.1705C189.477 26.1705 190.932 25.1193 191.273 23.5L189.034 23.3523C188.79 24.017 188.165 24.3636 187.358 24.3636C186.148 24.3636 185.381 23.5625 185.381 22.2614V22.2557H191.324V21.5909C191.324 18.625 189.528 17.1591 187.222 17.1591C184.653 17.1591 182.989 18.983 182.989 21.6761C182.989 24.4432 184.631 26.1705 187.318 26.1705ZM185.381 20.7557C185.432 19.7614 186.188 18.9659 187.261 18.9659C188.312 18.9659 189.04 19.7159 189.045 20.7557H185.381ZM195.324 14.3636H192.903V26H195.324V14.3636Z" fill="white"/>
                </svg>
            """
            vercel_div = Div(
                text=f"<a href='https://vercel.com/?utm_source=seqan&utm_campaign=oss'><svg height=40px>{vercel_logo}</svg></a>",
                styles={
                    "margin": "0px",
                    "height": "40px",
                    "display": "block" if "VERCEL" in os.environ else "none",
                    "width": "100vw",
                    "text-align": "center",
                },
            )
            all_elements = column(both_plots, vercel_div, sizing_mode="scale_both")
            tabs.append(TabPanel(child=all_elements, title=files_names_titles[i]))
    if interactive:
        show(Tabs(tabs=tabs, sizing_mode="scale_both"))
    else:
        tab_style = """
            :host(.bk-Tabs) {
            background-color: #15191c;
            }

            :host(.bk-Tabs) .bk-header {
            color: #e0e0e0;
            border-bottom: 1px solid #9c9c9c;
            }

            .bk-tab:hover {
                background-color: #7B7D7E;
            }
            """
        global_style = """
            body { background: #15191c; }
            div { max-height: 96vh; max-width: 100vw; }
        """
        save(Tabs(tabs=tabs, sizing_mode="scale_both", stylesheets=[tab_style, global_style]))


if __name__ == "__main__":
    create_plot()

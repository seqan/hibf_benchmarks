"""Functions to configure the style of the plots."""

from bokeh import events
from bokeh.models import (
    AdaptiveTicker,
    BoxZoomTool,
    CustomJS,
    CustomJSTickFormatter,
    FactorRange,
    HoverTool,
    Legend,
    LinearAxis,
    PanTool,
    ResetTool,
    WheelZoomTool,
)


def add_hover_tool(plot, renderer, key, display_key, file_name, format_kind):
    """Adds a hover tool to the plot."""
    if file_name in ("none", "U", "U+R"):
        file_name = "tmax"
    if format_kind == "TIME_FORMAT":
        plot.add_tools(
            HoverTool(
                tooltips=[
                    (file_name, "@value"),
                    ("Wall clock time", "@all_times{0.00} sek"),
                    (display_key, "@$name{0.00} sek"),
                    ("Percentage", f"@{key}_percentage{{0.00}}%"),
                ],
                renderers=[renderer],
                visible=False,
            )
        )
    else:
        plot.add_tools(
            HoverTool(
                tooltips=[
                    (file_name, "@value"),
                    ("Size", "@sizes{0.00}GiB"),
                    (display_key, "@$name{0.00}GiB"),
                    ("Percentage", f"@{key}_percentage{{0.00}}%"),
                    ("Load Factor (avg)", f"@{key}_avg_load_factor{{0.00}}"),
                ],
                renderers=[renderer],
                visible=False,
            )
        )


def add_legend(plot, renderers, file_name, time_names, size_names, format_kind, location):
    """Adds a legend to the plot."""
    key_format = time_names if format_kind == "TIME_FORMAT" else size_names
    legend_items = []
    for i, renderer in enumerate(renderers):
        display_key = key_format[i]
        key = renderer.name
        legend_items.append((display_key, [renderer]))
        add_hover_tool(plot, renderer, key, display_key, file_name, format_kind)
    legend = Legend(items=legend_items)
    legend.location = "left"
    legend.orientation = "vertical"
    legend.glyph_height = 4
    legend.glyph_width = 12
    legend.spacing = 0
    legend.click_policy = "mute"
    plot.add_layout(legend, location)
    toggle_legend_js = CustomJS(
        args={"legend": legend},
        code="""
            legend.visible = !legend.visible
            """,
    )
    plot.js_on_event(events.DoubleTap, toggle_legend_js)


def add_second_y_axis(plot, y_range):
    """Adds a second y-axis on the right side to the plot."""
    plot.extra_y_ranges = {"additional_axis": FactorRange(factors=y_range)}
    second_y_axis = LinearAxis(y_range_name="additional_axis")
    second_y_axis.major_label_text_font_size = "1pt"
    second_y_axis.major_label_text_color = "#15191c"
    plot.add_layout(second_y_axis, "right")


def configure_time_plot(plot, scale_in_minutes):
    """Configures the time plot."""
    if scale_in_minutes:
        plot.xaxis.ticker = AdaptiveTicker(base=60)
        plot.xaxis.axis_label = "time in minutes"
        plot.xaxis.formatter = CustomJSTickFormatter(code="return (tick / 60);")
    else:
        plot.xaxis.ticker = AdaptiveTicker(base=10)
        plot.xaxis.axis_label = "time in seconds"
    plot.toolbar.logo = None
    plot.toolbar_location = "below"
    plot.toolbar.autohide = True
    zoom_tool = WheelZoomTool(maintain_focus=False)
    plot.add_tools(PanTool(), zoom_tool, BoxZoomTool(), ResetTool())
    plot.toolbar.active_scroll = zoom_tool
    plot.x_range.bounds = (0, float("inf"))

    plot.yaxis.visible = False
    plot.y_range.range_padding = 0.1
    plot.ygrid.grid_line_color = None
    plot.sizing_mode = "scale_both"


def configure_size_plot(plot):
    """Configures the size plot."""
    plot.toolbar.logo = None
    plot.toolbar.autohide = True
    plot.toolbar_location = "below"
    plot.toolbar.autohide = True
    zoom_tool = WheelZoomTool(maintain_focus=False)
    plot.add_tools(PanTool(), zoom_tool, BoxZoomTool(), ResetTool())
    plot.toolbar.active_scroll = zoom_tool
    plot.x_range.bounds = (0, float("inf"))

    plot.y_range.range_padding = 0.1
    plot.ygrid.grid_line_color = None
    plot.axis.minor_tick_line_color = None
    plot.outline_line_color = None
    plot.axis.minor_tick_line_color = None
    plot.yaxis.major_tick_line_color = None
    plot.yaxis.minor_tick_line_color = None
    plot.yaxis.major_label_standoff = 15
    plot.sizing_mode = "scale_both"

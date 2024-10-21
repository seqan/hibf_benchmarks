"""Functions for creating CSS and HTML components for the Bokeh plot."""

# Todo: Config
descriptions = {
    "Alpha (&alpha;)": "Influences the ratio of merged bins and split bins",
    "Hash": "The number of hash functions for Bloom Filters",
    "k-mer": "Choosing window and k-mer size",
    "fpr": "Sets an upper bound for Bloom Filter false positives",
    "Estimate Union (U)": "Algorithm estimates sequence similarity between input data",
    "Rearrangement (U+R)": "Change order of sequences based on their estimated similarity",
}

default_params = {
    "&alpha;": 1.2,
    "t<sub>max</sub>": 192,
    "hash": 2,
    "k-mer": "32",
    "r-fpr": 0.5,
    "fpr": 0.05,
    "U": True,
    "R": True,
}

legend_desc = {
    "Determine query length": "Time to determine the length of the query",
    "Queryfile IO": "Time to read the query file",
    "Load index": "Time to load the index",
    "Compute minimizer (max)": "Average time to compute minimizer per thread",
    "Query IBF (max)": "Average time to query the IBF per thread",
    "Generate results (max)": "Average time to generate results per thread",
    "Level $$n$$": "Size of the $$n$$-th level of the index",
    "Avg load factor": "Average load factor of the index",
}

dataset = {
    "Description": "Simulated dataset with 1 haplotype and 1M reads",
    "Parameters": {
        "Type": "Simulated",
        "Sequence size": "512 MiB",
        "Number of bins": 1024,
        "Number of haplotypes": 1,
        "Number of reads": 1048576,
        "Read length": 250,
        "Read errors": 2,
    },
}


def create_latex_text():
    """Creates a div containing the description of the plot."""
    desc = "<div><h2>Description:</h2><ul>"
    for key, value in descriptions.items():
        desc += f"<li><span><strong>{key}:</strong> {value}</span></li>"
    desc += "</ul>"
    desc += "<h2>Default parameters:</h2><ul>"
    for key, value in default_params.items():
        desc += f"<li><span><strong>{key}</strong> = {value}</span></li>"
    desc += "</ul>"
    desc += "<h2>Legend:</h2><ul>"
    for key, value in legend_desc.items():
        desc += f"<li><span><strong>{key}:</strong> {value}</span></li>"
    desc += "</ul></div>"
    return desc


def create_dataset_text():
    """Creates a div containing the description of the plot."""
    desc = f'<div id="dataset" style="margin: 5px;"><h2 style="margin-bottom: 0;">Dataset:</h2><h4 style="margin-top: 0;"><strong>{dataset["Description"]}</strong></h4><div style="display: table; border-spacing: 2px;">'
    for key, value in dataset["Parameters"].items():
        desc += f'<div style="display: table-row;"><div style="display: table-cell; text-align: right; font-weight: bold">{key}:</div><div style="display: table-cell; text-align: left;"> {value}</div></div>'
    desc += "</div></div>"
    return desc


def get_tab_style():
    """Returns the CSS style for the tabs."""
    return """
        :host(.bk-Tabs) {
            background-color: #15191c;
        }

        .bk-tab {
            border-right: 1px solid #404040;
            border-width: 3px 1px 0px 1px;
            border-color: #15191c #404040 #404040 #404040;
            border-radius: var(--border-radius) var(--border-radius) 0 0;
            transition: background-color 0.1s ease-in-out, border-color 0.1s ease-in-out, color 0.1s ease-in-out;
        }

        :host(.bk-Tabs) .bk-header {
            color: #d0d0d0;
            border-bottom: 1px solid #404040;
            font-size: 1.1em;
        }

        .bk-tab.bk-active {
            background-color: #d0d0d0;
            border-color: #d0d0d0;
            color: black;
        }

        .bk-tab:not(.bk-active):hover {
            background-color: #404040;
        }

        .bk-tab:focus {
            outline: none;
        }

        .bk-tab:active {
            outline: none;
        }
        """


def get_button_style():
    """Returns the CSS style for the toggle button."""
    return [
        """
        .bk-btn, .bk-btn-success, .bk-btn:hover, .bk-btn:active, .bk-btn:focus {
            background: #500000;
            border: solid 1px #777777;
            color: #d0d0d0;
            cursor: pointer;
            outline: none;
            box-shadow: none;
        }
        .bk-btn.bk-btn-success.bk-active {
            background: #003509;
            color: #d0d0d0;
            border: solid 1px #777777;
        }
        .bk-btn:hover, .bk-btn.bk-btn-success.bk-active:hover {
            border: solid 1px #888888;
        }
    """
    ]


def get_global_style():
    """Returns the global CSS style."""
    return """
        body { background-color: #15191c; }
        div { max-height: 96vh; max-width: 100vw; }
        """


def get_hover_code():
    """Returns the JavaScript code for changing the hover tooltips."""
    return """
        var current_state = button.active;
        var time_description, size_description;

        if (current_state === false) {
            time_description = normal_time_description_list;
            size_description = normal_size_description_list;
        } else {
            time_description = advanced_time_description_list;
            size_description = advanced_size_description_list;
        }

        for (var i = 0; i < time_plot_hovers.length; i++) {
            time_plot_hovers[i].tooltips = time_description[i];
        }
        for (var i = 0; i < size_plot_hovers.length; i++) {
            size_plot_hovers[i].tooltips = size_description[i];
        }
        """


def landing_page_css():
    """Returns the CSS style for the landing page."""
    return """
    body {
        font-family: Arial, sans-serif;
        background-color: #15191c;
        margin: 0;
        padding: 20px;
    }

    .header {
        text-align: center;
        padding: 20px;
    }

    .header h1 {
        margin: 0;
        font-size: 2em;
        color: #ffffff;
    }

    .gallery {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
        gap: 20px;
        padding: 20px;
    }

    .gallery-item {
        height: 200px;
        background-color: #000000;
        border-radius: 8px;
        box-shadow: 0 0 10px rgba(255, 255, 255, 0.1);
        overflow: hidden;
        position: relative;
        text-align: center;
        display: flex;
        transform: scale(0.97);
        transition: 300ms ease-in-out;
    }

    .gallery-item:hover {
        box-shadow: 0 0 20px rgba(255, 255, 255, 0.134);
        transform: scale(1);
    }

    .gallery-item img {
        max-width: 100%;
        height: auto;
        display: block;
    }

    .gallery-item h4 {
        margin: 0;
        padding: 20px;
        text-align: center;
        display: flex;
        align-items: center;
        justify-content: center;
        color: #ffffff;
    }


    .description-box {
        position: absolute;
        height: calc(100% - 20px);
        bottom: 0;
        left: 0;
        right: 0;
        background: rgba(0, 0, 0, 0.8);
        color: white;
        padding: 10px;
        transform: translateY(100%);
        transition: transform 0.3s ease-in-out;
        display: flex;
        flex-direction: column;
        justify-content: center;
    }

    .gallery-item:hover .description-box {
        transform: translateY(0);
    }

    .description {
        max-height: 100%;
        overflow: auto;
        scrollbar-width: none;
        -ms-overflow-style: none;
        -webkit-overflow-scrolling: touch;
    }

    .description::-webkit-scrollbar {
        display: none;
    }

    .gallery-item a {
        text-decoration: none;
        color: #ffffff;
        display: block;
        padding: 10px;
        font-size: 1em;
        transition: 400ms ease-in-out;
    }

    .gallery-item a:hover {
        color: #c2c2c2;
    }

    .gallery-item h4 {
        font-size: 1.3em;
        margin-bottom: 5px;
        width: 100%;
    }
    """

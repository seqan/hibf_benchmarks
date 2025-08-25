"""
Creates a landing page gallery with screenshots of Bokeh plot HTML files.
"""

import os
import re

from bokeh_plot.components.log_init import log_init
from html2image import Html2Image

html_files = snakemake.input["PLOT_FILE"]
output_file = snakemake.output["OUTPUT_FILE"]
extra_file_plotting = snakemake.params["EXTRA_FILE_PLOTTING"]
html_dir = snakemake.params["HTML_DIR"]

log_init(snakemake.log[0])

# get html files
if extra_file_plotting:
    html_files = [os.path.join(html_dir, f) for f in os.listdir(html_dir) if f.endswith(".html")]

# get html names
html_names = [re.sub(".html", "", html_file) for html_file in html_files]

# create png for each html
hti = Html2Image(output_path=html_dir, custom_flags=["--headless", "--disable-gpu"])
if os.getenv("CI", "false") == "true":
    hti.browser.flags += ["--no-sandbox"]

for html_file in html_files:
    hti.screenshot(html_file=html_file, save_as=re.sub(".html", ".png", os.path.basename(html_file)))

# all parts of the landing page
LIST_OF_PARTS = "\n".join(
    [
        f"""
        <div class="gallery-item">
            <a href="{os.path.basename(html_name)}.html">
                <img src="{os.path.basename(html_name)}.png" alt="Bokeh Plot {ihtml_name+1}">
                Plot {ihtml_name+1}
            </a>
        </div>
        """
        for (ihtml_name, html_name) in enumerate(html_names)
    ]
)

# create landing page
CSS_TEXT = """
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
        grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
        gap: 20px;
        padding: 20px;
    }

    .gallery-item {
        background-color: #000000;
        border-radius: 8px;
        box-shadow: 0 0 10px rgba(255, 255, 255, 0.1);
        overflow: hidden;
        text-align: center;
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
        margin-bottom: 10px;
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
    """

HTML_TEXT = (
    """
    <!DOCTYPE html>
    <html lang="en">

    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Landing Page for Bokeh Plots</title>
        <link rel="stylesheet" href="style.css">
    </head>

    <body>
        <div class="header">
            <h1>Plot Gallery</h1>
        </div>
        <div class="gallery">
    """
    + LIST_OF_PARTS
    + """
        </div>
    </body>
    </html>
    """
)

# save landing page
with open(os.path.join(os.path.dirname(output_file), "style.css"), "w", encoding="utf-8") as f:
    f.write(CSS_TEXT)

with open(output_file, "w", encoding="utf-8") as f:
    f.write(HTML_TEXT)

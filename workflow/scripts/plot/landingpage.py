"""
Creates a landing page gallery with screenshots of Bokeh plot HTML files.
"""

import json
import os

from bs4 import BeautifulSoup

from components.log_init import log_init
from components.plot_css_html import landing_page_css

html_files = snakemake.input["PLOT_FILE"]  # type: ignore
output_file = snakemake.output["OUTPUT_FILE"]  # type: ignore
extra_file_plotting = snakemake.params["EXTRA_FILE_PLOTTING"]  # type: ignore
html_dir = snakemake.params["HTML_DIR"]  # type: ignore

log_init(snakemake.log[0])  # type: ignore


def find_texts(obj, texts):
    """recursiv text extraction from json"""
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key == "text":
                texts.append(value)
            elif isinstance(value, (dict, list)):
                find_texts(value, texts)
    elif isinstance(obj, list):
        for item in obj:
            if isinstance(item, (dict, list)):
                find_texts(item, texts)


def get_html_name(html_file):
    """get the html name from the html file"""
    return html_file.split("/")[-1].replace(".html", "")


def clean_html(html_file):
    """prepare html for extraction"""
    with open(html_file, "r", encoding="utf-8") as f:
        html = str(f.read())
    soup = BeautifulSoup(html, "html.parser")
    extracted_script = soup.find("script", {"type": "application/json"})
    converted_html = extracted_script.get_text().replace("&lt;", "<").replace("&gt;", ">").replace('"', '"')
    json_data = json.loads(converted_html)
    texts = []
    find_texts(json_data, texts)
    text_string = "\n".join(texts)
    converted_soup = BeautifulSoup(text_string, "html.parser")
    dataset = converted_soup.find("div", {"id": "dataset"})
    return dataset


def extract_description(html_file):
    """extract dataset details from html file"""
    cleaned_html = clean_html(html_file)
    if cleaned_html:
        for headline in cleaned_html.find_all("h2") + cleaned_html.find_all("h4"):
            headline.decompose()
    filename = get_html_name(html_file)
    return cleaned_html if cleaned_html else filename


def extract_headline(html_file):
    """extract headline from html file"""
    cleaned_html = clean_html(html_file)
    if cleaned_html:
        for headline in cleaned_html.find_all("h4"):
            return headline.get_text() if headline else html_file
    filename = get_html_name(html_file)
    return filename


if extra_file_plotting:
    """get html files"""
    html_files = [os.path.join(html_dir, f) for f in os.listdir(html_dir) if f.endswith(".html")]


# get html names
html_names = [get_html_name(html_file) for html_file in html_files]


# all gallery items for the landing page, headline of the data-description will become the title of the gallery-item,
# the description of the data-set will become the hover box
LIST_OF_PARTS = "\n".join(
    [
        f"""
        <div class="gallery-item">
            <h4>{extract_headline(f"results/html/{html_name}.html")}</h4>
            <a href="{html_name}.html">
                <div class="description-box">
                    <div class="description">
                        {extract_description(f"results/html/{html_name}.html")}
                    </div>
                </div>
            </a>
        </div>
        """
        for html_name in html_names
    ]
)


# html template
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


# save landing page and css
with open(os.path.join(os.path.dirname(output_file), "style.css"), "w", encoding="utf-8") as f:
    f.write(landing_page_css())

with open(output_file, "w", encoding="utf-8") as f:
    f.write(HTML_TEXT)

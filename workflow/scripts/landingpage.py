"""
Creates a landing page gallery with screenshots of Bokeh plot HTML files.
"""

import os
import json
from bs4 import BeautifulSoup
from bokeh_plot.components.log_init import log_init
from html2image import Html2Image

html_files = snakemake.input["PLOT_FILE"] # type: ignore
output_file = snakemake.output["OUTPUT_FILE"] # type: ignore
extra_file_plotting = snakemake.params["EXTRA_FILE_PLOTTING"] # type: ignore
html_dir = snakemake.params["HTML_DIR"] # type: ignore

log_init(snakemake.log[0]) # type: ignore

# recursiv text extraction from json
def find_texts(obj, texts):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key == 'text':
                texts.append(value)
            elif isinstance(value, (dict, list)):
                find_texts(value, texts)
    elif isinstance(obj, list):
        for item in obj:
            if isinstance(item, (dict, list)):
                find_texts(item, texts)

# extract dataset details from html file
def extract_description(html_file):
    with open(html_file, "r", encoding="utf-8") as f:
        html = str(f.read())
    soup = BeautifulSoup(html, 'html.parser')
    extracted_script = soup.find("script", {"type": "application/json"})
    converted_html = extracted_script.get_text().replace("&lt;", "<").replace("&gt;", ">").replace('\"', '"')
    json_data = json.loads(converted_html)
    texts = []
    find_texts(json_data, texts)
    text_string = "\n".join(texts)
    converted_soup = BeautifulSoup(text_string, 'html.parser')
    dataset = converted_soup.find("div", {"id": "dataset"})
    if dataset:
        for headline in dataset.find_all("h2"):
            headline.decompose()
    print(dataset)
    return dataset if dataset else html_file.split("/")[-1].replace(".html", "")

# get html files
if extra_file_plotting:
    html_files = [os.path.join(html_dir, f) for f in os.listdir(html_dir) if f.endswith(".html")]

# get html names
html_names = [
    html_file.split("/")[-1].replace(".html", "")
    for html_file in html_files
]

# create png for each html
hti = Html2Image(output_path=html_dir, custom_flags=["--headless", "--disable-gpu"])

for html_file in html_files:
    hti.screenshot(html_file=html_file, save_as=os.path.basename(html_file).replace(".html", ".png"))

# all parts of the landing page
LIST_OF_PARTS = "\n".join(
    [
        f"""
        <div class="gallery-item">
            <a href="{html_name}.html">
                <img src="{html_name}.png" alt="Bokeh Plot {ihtml_name+1}">
                <div class="description-box">
                    <div class="description">
                        {extract_description(f"results/html/{html_name}.html")}
                    </div>
                </div>
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
        grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
        gap: 20px;
        padding: 20px;
    }

    .gallery-item {
        background-color: #000000;
        border-radius: 8px;
        box-shadow: 0 0 10px rgba(255, 255, 255, 0.1);
        overflow: hidden;
        position: relative;
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

# HIBF Benchmarks

🚧 Work in progress 🚧

## Quickstart

### Edit configuration

* General configuration: [config/config.yaml](config/config.yaml)
* Plot configuration: [config/plot_config.yaml](config/plot_config.yaml)
* Meta information: [workflow/scripts/bokeh_plot/components/plot_css_html.py](workflow/scripts/bokeh_plot/components/plot_css_html.py)

### Install

```bash
conda env create --channel conda-forge --name snakemake bioconda::snakemake r-base r-here r-yaml r-tidyr bokeh
```

### Run

```bash
conda activate snakemake
snakemake --cores 32
```

The results can be found in `results/html/index.html`.

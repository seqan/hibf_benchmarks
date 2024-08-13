# HIBF Benchmarks

🚧 Work in progress 🚧

## Quickstart

### Edit configuration

* General configuration: [config/config.yaml](config/config.yaml)
* Plot configuration: [config/plot_config.yaml](config/plot_config.yaml)
* Meta information: [workflow/scripts/bokeh_plot/components/plot_css_html.py](workflow/scripts/bokeh_plot/components/plot_css_html.py)

### Setup

* [Conda](https://docs.anaconda.com/miniconda/miniconda-install) or [Mamba](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

```bash
# Create an environment called `snakemake` with snakemake installed
conda env create --channel conda-forge --name snakemake bioconda::snakemake
```

### Run

```bash
conda activate snakemake
snakemake --use-conda --cores 2
```

The results can be found in `results/html/index.html`.

### Testing/Running without Raptor

To test the workflow without running Raptor, data is provided in `.github/data`:

```bash
# run from repository root directory
mkdir --parents results/raw_data
tar xf .github/data/data.tar.zst --directory=results/raw_data
export CI=true # Skips running raptor
snakemake --use-conda --cores 2 --forceall
```

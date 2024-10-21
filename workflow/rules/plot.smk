rule plot_data:
    input:
        SIZE_INPUT="results/summarized/size.tsv",
        TIME_INPUT="results/summarized/time.tsv",
    output:
        PLOT_FILE=f"results/html/{config['PLOT_NAME']}.html",
    params:
        THEME="workflow/scripts/plot/plot_theme.yaml",
        KEYS=config["KEYS"],
        TIME=config["TIME"],
        SIZE=config["SIZE"],
    log:
        "logs/plot_data/plot_data.log",
    conda:
        "../envs/bokeh.yml"
    script:
        "../scripts/plot/plot.py"


rule plot_landingpage:
    input:
        PLOT_FILE=f"results/html/{config['PLOT_NAME']}.html",
    output:
        OUTPUT_FILE="results/html/index.html",
    params:
        EXTRA_FILE_PLOTTING=config["EXTRA_FILE_PLOTTING"],
        HTML_DIR="results/html",
    log:
        "logs/plot_landingpage/plot_landingpage.log",
    conda:
        "../envs/landingpage.yml"
    script:
        "../scripts/plot/landingpage.py"

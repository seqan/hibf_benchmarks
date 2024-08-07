rule plot_data:
    input:
        SIZE_INPUT="results/size",
        TIME_INPUT="results/time",
    output:
        PLOT_FILE=f"results/html/{config['PLOT_NAME']}.html",
    params:
        THEME="workflow/scripts/bokeh_plot/plot_theme.yaml",
        KEYS=config["KEYS"],
        TIME=config["TIME"],
        SIZE=config["SIZE"],
    log:
        "logs/plot_data/plot_data.log",
    conda:
        "../envs/bokeh.yml"
    script:
        "../scripts/bokeh_plot/plot.py"


rule plot_landingpage:
    input:
        PLOT_FILE=f"results/html/{config['PLOT_NAME']}.html",
    output:
        OUTPUT_FILE="results/html/index.html",
        PNG_FILE=f"results/html/{config['PLOT_NAME']}.png",
    params:
        EXTRA_FILE_PLOTTING=config["EXTRA_FILE_PLOTTING"],
        HTML_DIR="results/html",
    log:
        "logs/plot_landingpage/plot_landingpage.log",
    conda:
        "../envs/landingpage.yml"
    script:
        "../scripts/landingpage.py"

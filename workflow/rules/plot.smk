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

rule extract_data_for_plot:
    input:
        SIZE_INPUT="results/size",
        TIME_INPUT="results/time",
    output:
        expand(
            "results/prepared_{format}/{key}",
            format=["time", "size"],
            key=config["KEYS"],
        ),
    log:
        "logs/extract_data_for_plot/extract_data_for_plot.log",
    conda:
        "../envs/r.yml"
    script:
        "../scripts/extract_results.R"


rule plot_data:
    input:
        expand(
            "results/prepared_{format}/{key}",
            format=["time", "size"],
            key=config["KEYS"],
        ),
    output:
        "results/html/index.html",
    log:
        "logs/plot_data/plot_data.log",
    conda:
        "../envs/bokeh.yml"
    script:
        "../scripts/bokeh_plot/plot.py"

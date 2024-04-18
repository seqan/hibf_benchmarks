rule extract_data_for_plot:
    input:
        SIZE_INPUT=f"{config['BUILD_DIR']}/size",
        TIME_INPUT=f"{config['BUILD_DIR']}/time",
    output:
        expand(
            f"{config['BUILD_DIR']}/prepared_{{format}}/{{key}}",
            format=["time", "size"],
            key=config["KEYS"],
        ),
    log:
        f"{config['LOG_DIR']}/extract_data_for_plot.log",
    conda:
        "../../envs/r.yml"
    script:
        "../extract_results.R"


rule plot_data:
    input:
        expand(
            f"{config['BUILD_DIR']}/prepared_{{format}}/{{key}}",
            format=["time", "size"],
            key=config["KEYS"],
        ),
    output:
        f"{config['PLOT_FILE']}",
    log:
        f"{config['LOG_DIR']}/plot_data.log",
    conda:
        "../../envs/bokeh.yml"
    script:
        "../bokeh_plot/plot.py"

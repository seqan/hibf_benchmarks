rule extract_data_for_plot:
    input:
        SIZE_INPUT=f"{config['BUILD_DIR']}/size",
        TIME_INPUT=f"{config['BUILD_DIR']}/time",
    output:
        expand(
            f"{config['BUILD_DIR']}/prepared_{{format}}/{{key}}",
            format=["time", "size"],
            key=plot_config["KEYS_FORMAT"],
        ),
    params:
        BUILD_DIR=config["BUILD_DIR"],
        TIME_FORMAT=config["TIME_FORMAT"],
        SIZE_FORMAT=config["SIZE_FORMAT"],
        KEYS_FORMAT=plot_config["KEYS_FORMAT"],
    log:
        f"{config['LOG_DIR']}/extract_data_for_plot.log",
    conda:
        "../../envs/r_basic_env.yaml"
    script:
        "../extract_results.R"


rule plot_data:
    input:
        expand(
            f"{config['BUILD_DIR']}/prepared_{{format}}/{{key}}",
            format=["time", "size"],
            key=plot_config["KEYS_FORMAT"],
        ),
    output:
        f"{plot_config['PLOT_FILE']}",
    params:
        BUILD_DIR=config["BUILD_DIR"],
        PLOT_FILE=plot_config["PLOT_FILE"],
        THEME=plot_config["THEME"],
        KEYS_FORMAT=plot_config["KEYS_FORMAT"],
        KEYS_NAMES=plot_config["KEYS_NAMES"],
        TIME_FORMAT=plot_config["TIME_FORMAT"],
        TIME_NAMES=plot_config["TIME_NAMES"],
        SIZE_FORMAT=plot_config["SIZE_FORMAT"],
        SIZE_NAMES=plot_config["SIZE_NAMES"],
    log:
        f"{config['LOG_DIR']}/plot_data.log",
    conda:
        "../../envs/bokeh_env.yaml"
    script:
        "../bokeh_plot/plot.py"

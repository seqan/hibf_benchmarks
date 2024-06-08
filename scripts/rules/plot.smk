# rule extract_data_for_plot:
#     input:
#         SIZE_INPUT=f"{config['BUILD_DIR']}/size",
#         TIME_INPUT=f"{config['BUILD_DIR']}/time",
#     output:
#         expand(
#             f"{config['BUILD_DIR']}/prepared_{{format}}/{{key}}",
#             format=["time", "size"],
#             key=plot_config["KEYS_FORMAT"],
#         ),
#     params:
#         BUILD_DIR=config["BUILD_DIR"],
#         TIME_FORMAT=config["TIME_FORMAT"],
#         SIZE_FORMAT=config["SIZE_FORMAT"],
#         KEYS_FORMAT=plot_config["KEYS_FORMAT"],
#     log:
#         f"{config['LOG_DIR']}/extract_data_for_plot.log",
#     conda:
#         "../../envs/r_basic_env.yaml"
#     script:
#         "../extract_results.R"


rule plot_data:
    input:
        SIZE_INPUT=f"{config['BUILD_DIR']}/size",
        TIME_INPUT=f"{config['BUILD_DIR']}/time",
    output:
        PLOT_FILE = f"{plot_config['PLOT_FILE']}",
    params:
        THEME=plot_config["THEME"],
        KEYS=plot_config["KEYS"],
        TIME=plot_config["TIME"],
        SIZE=plot_config["SIZE"],
    log:
        f"{config['LOG_DIR']}/plot_data.log",
    conda:
        "../../envs/bokeh_env.yaml"
    script:
        "../bokeh_plot/plot.py"

include: "common.py"


rule summarize_timings:
    input:
        INPUT_FILES=expand(
            f"{config['BUILD_DIR']}/{{param}}/out.time",
            param=get_params(True),
        ),
    output:
        OUTPUT_FILE=f"{config['BUILD_DIR']}/time",
    params:
        FORMAT=config["TIME_FORMAT"],
    log:
        f"{config['LOG_DIR']}/summarize_timings.log",
    conda:
        "../../envs/r_basic_env.yaml"
    script:
        "../summarize_results.R"


rule summarize_sizes:
    input:
        INPUT_FILES=expand(f"{config['BUILD_DIR']}/{{param}}/out.sizes", param=get_params(False)),
    output:
        OUTPUT_FILE=f"{config['BUILD_DIR']}/size",
    params:
        FORMAT=config["SIZE_FORMAT"],
    log:
        f"{config['LOG_DIR']}/summarize_sizes.log",
    conda:
        "../../envs/r_basic_env.yaml"
    script:
        "../summarize_results.R"

rule raptor_layout:
    output:
        LAYOUT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
        LAYOUT_TIME=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout.time",
    threads: config["NUM_THREADS"]
    priority: 2
    message:
        "Running raptor build for {wildcards.key}={wildcards.param}."
    log:
        f"{config['LOG_DIR']}/raptor_layout/{{key}}_{{param}}.log",
    conda:
        "../../envs/raptor_env.yaml"
    script:
        "run_layout.py"


rule raptor_build:
    input:
        f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
    output:
        INDEX_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/index",
        INDEX_TIME=f"{config['BUILD_DIR']}/{{key}}={{param}}/index.time",
    threads: config["NUM_THREADS"]
    priority: 1
    log:
        f"{config['LOG_DIR']}/raptor_build/{{key}}_{{param}}.log",
    conda:
        "../../envs/raptor_env.yaml"
    params:
        RAPTOR_BINARY=config["RAPTOR_BINARY"],
        WINDOW_SIZE=config["DEFAULT_PARAMS"]["WINDOW_SIZE"],
    shell:
        """
    (
    {params.RAPTOR_BINARY} build \
      --input {input} \
      --output {output.INDEX_FILE} \
      --window {params.WINDOW_SIZE} \
      --quiet \
      --timing-output {output.INDEX_TIME} \
      --threads {threads}
    ) &>> {log}
    """


rule raptor_search:
    input:
        INDEX_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/index",
    output:
        RESULT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/out",
        RESULT_TIME=f"{config['BUILD_DIR']}/{{key}}={{param}}/out.time",
    threads: config["NUM_THREADS"]
    log:
        f"{config['LOG_DIR']}/raptor_search/{{key}}_{{param}}.log",
    conda:
        "../../envs/raptor_env.yaml"
    params:
        RAPTOR_BINARY=config["RAPTOR_BINARY"],
        QUERY_FILE=config["QUERY_FILE"],
        QUERY_ERRORS=config["DATA_PARAMETERS"]["QUERY_ERRORS"],
    shell:
        """
    (
    {params.RAPTOR_BINARY} search \
      --index {input.INDEX_FILE} \
      --query {params.QUERY_FILE} \
      --output {output.RESULT_FILE} \
      --error {params.QUERY_ERRORS} \
      --quiet \
      --timing-output {output.RESULT_TIME} \
      --threads {threads}
    ) &>> {log}
    """


rule display_layout:
    input:
        LAYOUT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
    output:
        SIZE_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/out.sizes",
    threads: config["NUM_THREADS"]
    log:
        f"{config['LOG_DIR']}/display_layout/{{key}}_{{param}}.log",
    conda:
        "../../envs/raptor_env.yaml"
    params:
        DISPLAY_LAYOUT_BINARY=config["DISPLAY_LAYOUT_BINARY"],
        WINDOW_SIZE=config["DEFAULT_PARAMS"]["WINDOW_SIZE"],
    shell:
        """
    (
    {params.DISPLAY_LAYOUT_BINARY} sizes \
      --input {input.LAYOUT_FILE} \
      --output {output.SIZE_FILE} \
      --threads {threads}
    ) &>> {log}
    """

rule raptor_layout:
    output:
        LAYOUT_FILE="results/raw_data/{key}={param}/layout",
        LAYOUT_TIME="results/raw_data/{key}={param}/layout.time",
    threads: config["NUM_THREADS"]
    priority: 2
    message:
        "Running raptor build for {wildcards.key}={wildcards.param}."
    log:
        "logs/raptor_layout/{key}_{param}.log",
    conda:
        "../envs/raptor.yml"
    script:
        "../scripts/run_layout.py"


rule raptor_build:
    input:
        "results/raw_data/{key}={param}/layout",
    output:
        INDEX_FILE="results/raw_data/{key}={param}/index",
        INDEX_TIME="results/raw_data/{key}={param}/index.time",
    threads: config["NUM_THREADS"]
    priority: 1
    log:
        "logs/raptor_build/{key}_{param}.log",
    conda:
        "../envs/raptor.yml"
    params:
        RAPTOR_BINARY=config["RAPTOR_BINARY"],
    shell:
        """
    (
    {params.RAPTOR_BINARY} build \
      --input {input} \
      --output {output.INDEX_FILE} \
      --timing-output {output.INDEX_TIME} \
      --threads {threads}
    ) &>> {log}
    """


rule raptor_search:
    input:
        INDEX_FILE="results/raw_data/{key}={param}/index",
    output:
        RESULT_FILE="results/raw_data/{key}={param}/out",
        RESULT_TIME="results/raw_data/{key}={param}/out.time",
    threads: config["NUM_THREADS"]
    log:
        "logs/raptor_search/{key}_{param}.log",
    conda:
        "../envs/raptor.yml"
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
      --timing-output {output.RESULT_TIME} \
      --threads {threads}
    ) &>> {log}
    """


rule display_layout:
    input:
        LAYOUT_FILE="results/raw_data/{key}={param}/layout",
    output:
        SIZE_FILE="results/raw_data/{key}={param}/out.sizes",
    threads: config["NUM_THREADS"]
    log:
        "logs/display_layout/{key}_{param}.log",
    conda:
        "../envs/raptor.yml"
    params:
        DISPLAY_LAYOUT_BINARY=config["DISPLAY_LAYOUT_BINARY"],
    shell:
        """
    (
    {params.DISPLAY_LAYOUT_BINARY} sizes \
      --input {input.LAYOUT_FILE} \
      --output {output.SIZE_FILE} \
      --threads {threads}
    ) &>> {log}
    """

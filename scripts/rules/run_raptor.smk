rule raptor_layout:
    output:
        LAYOUT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
        LAYOUT_TIME=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout.time",
    threads: config["NUM_THREADS"]
    priority: 2
    log:
        f"{config['LOG_DIR']}/raptor_layout/{{key}}_{{param}}.log",
    conda:
        "../../envs/raptor_env.yaml"
    params:
        RAPTOR_BINARY=config["RAPTOR_BINARY"],
        FILENAMES_FILE=config["FILENAMES_FILE"],
        KMER=config["DEFAULT_PARAMS"]["KMER_SIZE"],
        WINDOW_SIZE=config["DEFAULT_PARAMS"]["WINDOW_SIZE"],
        T_MAX=config["DEFAULT_PARAMS"]["T_MAX"],
        R_FPR=config["DEFAULT_PARAMS"]["R_FPR"],
        NUM_HASHES=config["DEFAULT_PARAMS"]["NUM_HASHES"],
        MAXIMUM_FPR=config["DEFAULT_PARAMS"]["M_FPR"],
        ALPHA=config["DEFAULT_PARAMS"]["ALPHA"],
        MODE=config["DEFAULT_PARAMS"]["MODE"],
    shell:
        """(
        kmer=$([[ {wildcards.key} == "kmer" ]] && echo {wildcards.param} || echo {params.KMER})
        mode=$([[ {wildcards.key} == "U+R" || {wildcards.key} == "U" || {wildcards.key} == "none" ]] && echo {wildcards.key} || echo {params.MODE})
        tmax=$([[ {wildcards.key} == "U+R" || {wildcards.key} == "U" || {wildcards.key} == "none" ]] && echo {wildcards.param} || echo {params.T_MAX})
        hash=$([[ {wildcards.key} == "hash" ]] && echo {wildcards.param} || echo {params.NUM_HASHES})
        relaxed_fpr=$([[ {wildcards.key} == "r-relaxed" ]] && echo {wildcards.param} | sed 's/_/./g' || echo {params.R_FPR})
        alpha=$([[ {wildcards.key} == "alpha" ]] && echo {wildcards.param} | sed 's/_/./g' || echo {params.ALPHA})
        disable_estimate_union=''
        disable_rearrangement=''
        if [[ $mode == 'U' ]]; then
        disable_rearrangement='--disable-rearrangement'
        elif [[ $mode == 'none' ]]; then
        disable_estimate_union='--disable-estimate-union'
        fi
        echo "[$(date +"%Y-%m-%d %T")] Running raptor layout for {wildcards.key}={wildcards.param}."
        {params.RAPTOR_BINARY} layout \
        --input {params.FILENAMES_FILE} \
        --output {output.LAYOUT_FILE} \
        --kmer $kmer \
        --window {params.WINDOW_SIZE} \
        --tmax $tmax \
        --fpr {params.MAXIMUM_FPR} \
        --relaxed-fpr $relaxed_fpr \
        --hash $hash \
        $disable_estimate_union \
        $disable_rearrangement \
        --alpha $alpha \
        --max-rearrangement-ratio 0.5 \
        --sketch-bits 12 \
        --timing-output {output.LAYOUT_TIME} \
        --threads {threads}
        ) &>> {log}
        """


rule raptor_build:
    input:
        LAYOUT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
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
        (echo "[$(date +"%Y-%m-%d %T")] Running raptor build for {wildcards.key}={wildcards.param}."
        {params.RAPTOR_BINARY} build \
        --input {input.LAYOUT_FILE} \
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
        (echo "[$(date +"%Y-%m-%d %T")] Running raptor search for {wildcards.key}={wildcards.param}."
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
        (echo "[$(date +"%Y-%m-%d %T")] Running display_layout for {wildcards.key}={wildcards.param}."
        {params.DISPLAY_LAYOUT_BINARY} sizes \
        --input {input.LAYOUT_FILE} \
        --output {output.SIZE_FILE} \
        --threads {threads}
        ) &>> {log}
        """

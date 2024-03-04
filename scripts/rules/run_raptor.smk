rule raptor_layout:
    input:
        "raptor_binary_checked",
        "input_file_checked",
    output:
        LAYOUT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
        LAYOUT_TIME=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout.time",
    threads: config["NUM_THREADS"]
    log: "log/raptor_layout_{key}_{param}.log",
    conda: "../envs/raptor_env.yaml"
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
    shell:
        """
        (kmer=$([[ {wildcards.key} == "kmer" ]] && echo {wildcards.param} || echo {params.KMER})
        tmax=$([[ {wildcards.key} == "tmax" ]] && echo {wildcards.param} || echo {params.T_MAX})
        hash=$([[ {wildcards.key} == "hash" ]] && echo {wildcards.param} || echo {params.NUM_HASHES})
        relaxed_fpr=$([[ {wildcards.key} == "r-relaxed" ]] && echo {wildcards.param} | sed 's/_/./g' || echo {params.R_FPR})
        alpha=$([[ {wildcards.key} == "alpha" ]] && echo {wildcards.param} | sed 's/_/./g' || echo {params.ALPHA})
        echo "[$(date +"%Y-%m-%d %T")] Running raptor layout for {wildcards.key}={wildcards.param}."
        disable_estimate_union='--disable-estimate-union'
        disable_rearrangement='--disable-rearrangement'
        if [[ '{wildcards.key}' == 'U' || '{wildcards.key}' == 'U+R' ]]; then
            disable_estimate_union=''
        fi
        if [[ '{wildcards.key}' == 'U+R' ]]; then
            disable_rearrangement='--disable-rearrangement'
        fi
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
        "raptor_binary_checked",
        LAYOUT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
    output:
        INDEX_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/index",
        INDEX_TIME=f"{config['BUILD_DIR']}/{{key}}={{param}}/index.time",
    threads: config["NUM_THREADS"]
    log: "log/raptor_build_{key}_{param}.log",
    conda: "../envs/raptor_env.yaml"
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
        "raptor_binary_checked",
        "input_file_checked",
        INDEX_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/index",
    output:
        RESULT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/out",
        RESULT_TIME=f"{config['BUILD_DIR']}/{{key}}={{param}}/out.time",
    threads: config["NUM_THREADS"]
    log: "log/raptor_search_{key}_{param}.log",
    conda: "../envs/raptor_env.yaml"
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
        "display_layout_binary_checked",
        LAYOUT_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/layout",
    output:
        SIZE_FILE=f"{config['BUILD_DIR']}/{{key}}={{param}}/out.sizes",
    threads: config["NUM_THREADS"]
    log: "log/display_layout_{key}_{param}.log",
    conda: "../envs/raptor_env.yaml"
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
rule store_timings:
    input:
        files=expand(
            f"{config['BUILD_DIR']}/{{param}}/out.time",
            param=[f"alpha={str(param).replace('.', '_')}" for param in config["PARAMS"]["ALPHA"]]
            + [
                f"{param1}={param2}"
                for param1 in config["PARAMS"]["MODE"] or config["DEFAULT_PARAMS"]["MODE"]
                for param2 in config["PARAMS"]["T_MAX"] or config["DEFAULT_PARAMS"]["T_MAX"]
            ]
            + [f"hash={param}" for param in config["PARAMS"]["NUM_HASHES"] if param < 6]
            + [f"kmer={param}" for param in config["PARAMS"]["KMER_SIZE"]]
            + [f"relaxed-fpr={str(param).replace('.', '_')}" for param in config["PARAMS"]["RELAXED_FPR"]],
        ),
    output:
        f"{config['BUILD_DIR']}/time",
    log:
        "log/store_timings.log",
    conda:
        "../envs/r_basic_env.yaml"
    shell:
        """
        (echo "[$(date +"%Y-%m-%d %T")] Storing timings."
        Rscript summarize_results.r "time_format" {input.files}
        ) &>> {log}
        """


rule store_sizes:
    input:
        files=expand(
            f"{config['BUILD_DIR']}/{{param}}/out.sizes",
            param=[f"alpha={str(param).replace('.', '_')}" for param in config["PARAMS"]["ALPHA"]]
            + [
                f"{param1}={param2}"
                for param1 in config["PARAMS"]["MODE"] or config["DEFAULT_PARAMS"]["MODE"]
                for param2 in config["PARAMS"]["T_MAX"] or config["DEFAULT_PARAMS"]["T_MAX"]
            ]
            + [f"hash={param}" for param in config["PARAMS"]["NUM_HASHES"]]
            + [f"kmer={param}" for param in config["PARAMS"]["KMER_SIZE"]]
            + [f"relaxed-fpr={str(param).replace('.', '_')}" for param in config["PARAMS"]["RELAXED_FPR"]],
        ),
    output:
        f"{config['BUILD_DIR']}/size",
    log:
        "log/store_sizes.log",
    conda:
        "../envs/r_basic_env.yaml"
    shell:
        """
        (echo "[$(date +"%Y-%m-%d %T")] Storing sizes."
        Rscript summarize_results.r "size_format" {input.files}
        ) &>> {log}
        """

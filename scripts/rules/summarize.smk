rule summerize_timings:
    input:
        INPUT_FILES=expand(
            f"{config['BUILD_DIR']}/{{param}}/out.time",
            param=[f"alpha={str(param).replace('.', '_')}" for param in config["PARAMS"]["ALPHA"]]
            + [f"{param1}={param2}"
                for param1 in config["PARAMS"]["MODE"] or config["DEFAULT_PARAMS"]["MODE"]
                for param2 in config["PARAMS"]["T_MAX"] or config["DEFAULT_PARAMS"]["T_MAX"]]
            + [f"hash={param}" for param in config["PARAMS"]["NUM_HASHES"] if param < 6]
            + [f"kmer={param}" for param in config["PARAMS"]["KMER_SIZE"]]
            + [f"relaxed-fpr={str(param).replace('.', '_')}" for param in config["PARAMS"]["RELAXED_FPR"]],
        ),
    output:
        OUTPUT_FILE = f"{config['BUILD_DIR']}/time",
    log:
        f"{config['LOG_DIR']}/summerize_timings.log",
    conda:
        "../../envs/r_basic_env.yaml"
    script:
        "../summarize_results.R"


rule summerize_sizes:
    input:
        INPUT_FILES=expand(
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
        OUTPUT_FILE=f"{config['BUILD_DIR']}/size",
    log:
        f"{config['LOG_DIR']}/summerize_sizes.log",
    conda:
        "../../envs/r_basic_env.yaml"
    script:
        "../summarize_results.R"

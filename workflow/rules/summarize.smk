include: "../scripts/common.py"


rule summarize_timings:
    input:
        INPUT_FILES=expand(
            "results/{param}/out.time",
            param=get_params(True),
        ),
    output:
        OUTPUT_FILE="results/time",
    log:
        "logs/summarize_timings/summarize_timings.log",
    conda:
        "../envs/r.yml"
    script:
        "../scripts/summarize_results.R"


rule summarize_sizes:
    input:
        INPUT_FILES=expand("results/{param}/out.sizes", param=get_params(False)),
    output:
        OUTPUT_FILE="results/size",
    log:
        "logs/summarize_sizes/summarize_sizes.log",
    conda:
        "../envs/r.yml"
    script:
        "../scripts/summarize_results.R"

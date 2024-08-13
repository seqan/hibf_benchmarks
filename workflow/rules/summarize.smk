include: "../scripts/common.py"


rule summarize_timings:
    input:
        INPUT_FILES=expand(
            "results/raw_data/{param}/out.time",
            param=get_params(True),
        ),
    output:
        OUTPUT_FILE="results/summarized/time.tsv",
    log:
        "logs/summarize_timings/summarize_timings.log",
    conda:
        "../envs/r.yml"
    script:
        "../scripts/summarize_results.R"


rule summarize_sizes:
    input:
        INPUT_FILES=expand("results/raw_data/{param}/out.sizes", param=get_params(False)),
    output:
        OUTPUT_FILE="results/summarized/size.tsv",
    log:
        "logs/summarize_sizes/summarize_sizes.log",
    conda:
        "../envs/r.yml"
    script:
        "../scripts/summarize_results.R"

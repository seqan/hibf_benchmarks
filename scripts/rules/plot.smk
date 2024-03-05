rule prepare_data_for_plot:
    input:
        "output_checked",
        SIZE_FILE=f"{config['BUILD_DIR']}/size",
        TIME_FILE=f"{config['BUILD_DIR']}/time",
    output:
        expand(
            f"{config['BUILD_DIR']}/prepared_{{format}}/{{key}}",
            format=["time", "size"],
            key=["alpha", "hash", "kmer", "relaxed-fpr", "none", "U", "U+R"],
        ),
    threads: config["NUM_THREADS"]
    log:
        "log/prepare_data_for_plot.log",
    conda:
        "../envs/r_basic_env.yaml"
    shell:
        """
        (echo "[$(date +"%Y-%m-%d %T")] Preparing size for plot."
        Rscript extract_results.r {input.TIME_FILE} {input.SIZE_FILE}
        ) &>> {log}
        """


rule plot_data:
    input:
        "output_checked",
        expand(
            f"{config['BUILD_DIR']}/prepared_{{format}}/{{key}}",
            format=["time", "size"],
            key=["alpha", "hash", "kmer", "relaxed-fpr", "none", "U", "U+R"],
        ),
    output:
        f"{config['BUILD_DIR']}/index.html",
    threads: config["NUM_THREADS"]
    log:
        "log/plot_data.log",
    conda:
        "../envs/bokeh_env.yaml"
    shell:
        """
        (echo "[$(date +"%Y-%m-%d %T")] Plotting data."
        python plot.py
        ) &>> {log}
        """

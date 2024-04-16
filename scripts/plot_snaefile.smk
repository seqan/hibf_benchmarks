configfile: "parameters.yaml"


include: "rules/summarize.smk"
include: "rules/plot.smk"


rule all:
    input:
        f"{config['PLOT_DIR']}/index.html",

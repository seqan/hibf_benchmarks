rule check_output_dir:
    output:
        (touch("output_checked")),
    log: "log/check_paths.log",
    params:
        OUTPUT_DIR=config['BUILD_DIR'],
    shell:
        """
        (if ! [[ -d {params.OUTPUT_DIR} ]]; then
            echo "Directory {params.OUTPUT_DIR} does not exist. Exiting."
            exit 1
        fi 
        echo "[$(date +"%Y-%m-%d %T")] Output directory {params.OUTPUT_DIR} exists."
        ) &>> {log}
        """


rule check_raptor_binary:
    output:
        (touch("raptor_binary_checked")),
    log: "log/check_paths.log",
    params:
        RAPTOR_BINARY=config["RAPTOR_BINARY"],
    shell:
        """
        (if ! [[ -s {params.RAPTOR_BINARY} ]]; then
            echo "Executable {params.RAPTOR_BINARY} does not exist. Exiting."
            exit 1
        fi
        if ! {params.RAPTOR_BINARY} --help &> /dev/null; then
            echo "Executable {params.RAPTOR_BINARY} cannot be run. Exiting."
            exit 1
        fi
        echo "[$(date +"%Y-%m-%d %T")] Raptor binary {params.RAPTOR_BINARY} looks good." 
        ) &>> {log}
        """


rule check_input_files:
    output:
        (touch("input_file_checked")),
    log: "log/check_paths.log",
    params:
        FILENAMES_FILE=config["FILENAMES_FILE"],
        QUERY_FILE=config["QUERY_FILE"],
    shell:
        """
        (if ! [[ -s {params.FILENAMES_FILE} ]]; then
            echo "Filenames file {params.FILENAMES_FILE} does not exist or is empty. Exiting."
            exit 1
        fi
        if ! [[ -s {params.QUERY_FILE} ]]; then
            echo "Query file {params.QUERY_FILE} does not exist or is empty. Exiting."
            exit 1
        fi
        echo "[$(date +"%Y-%m-%d %T")] Input files look good."
        ) &>> {log}
        """


rule display_layout_binary_check:
    output:
        (touch("display_layout_binary_checked")),
    log: "log/check_paths.log",
    params:
        DISPLAY_LAYOUT_BINARY=config["DISPLAY_LAYOUT_BINARY"],
    shell:
        """
        (if ! [[ -s {params.DISPLAY_LAYOUT_BINARY} ]]; then
            echo "Executable {params.DISPLAY_LAYOUT_BINARY} does not exist. Exiting."
            exit 1
        fi

        if ! {params.DISPLAY_LAYOUT_BINARY} --help &> /dev/null; then
            echo "Executable {params.DISPLAY_LAYOUT_BINARY} cannot be run. Exiting."
            exit 1
        fi
        echo "[$(date +"%Y-%m-%d %T")] Display layout binary {params.DISPLAY_LAYOUT_BINARY} looks good."
        ) &>> {log}
        """
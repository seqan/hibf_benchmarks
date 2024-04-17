#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euxo pipefail

set +x
source ~/.bashrc
conda activate
set -x

export WORK_DIR=`pwd`

wget --quiet --output-document yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64
chmod +x yq
./yq --inplace '.BUILD_DIR = strenv(WORK_DIR) + "/build"' ${WORK_DIR}/scripts/parameters.yaml
./yq --inplace '.LOG_DIR = strenv(WORK_DIR) + "/log"' ${WORK_DIR}/scripts/parameters.yaml
./yq --inplace '.PLOT_FILE = strenv(WORK_DIR) + "/build/html/index.html"' ${WORK_DIR}/scripts/plot_parameters.yaml
./yq --inplace '.THEME = strenv(WORK_DIR) + "/scripts/bokeh_plot/plot_theme.yaml"' ${WORK_DIR}/scripts/plot_parameters.yaml

mkdir -p ${WORK_DIR}/build
zstd --decompress ${WORK_DIR}/data/data.tar.zst -c | tar xf - -C ${WORK_DIR}/build

cd scripts
snakemake --use-conda --cores 2

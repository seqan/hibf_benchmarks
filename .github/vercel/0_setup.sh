#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euvo pipefail

MF_DIR="`pwd`/node_modules/miniforge"
MF_FILE="Miniforge3-$(uname)-$(uname -m).sh"
echo ${MF_DIR}
echo ${MF_FILE}

curl -LsO "https://github.com/conda-forge/miniforge/releases/latest/download/${MF_FILE}"
chmod +x ${MF_FILE}
./${MF_FILE} -b -u -p ${MF_DIR}
${MF_DIR}/bin/mamba init bash

set +v && source ~/.bashrc && set -v

conda config --set default_threads 4
conda config --set channel_priority strict
mamba env create --yes --file ./workflow/envs/snakemake.yml
mamba clean --all --yes

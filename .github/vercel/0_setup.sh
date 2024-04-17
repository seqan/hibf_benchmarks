#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euxo pipefail

WORK_DIR=`pwd`
CACHE_DIR="${WORK_DIR}/node_modules" # The node_modules directory is always cached.

mkdir -p ${CACHE_DIR}/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o ${CACHE_DIR}/miniconda3/miniconda.sh
chmod +x ${CACHE_DIR}/miniconda3/miniconda.sh
${CACHE_DIR}/miniconda3/miniconda.sh -b -u -p ${CACHE_DIR}/miniconda3
rm -rf ${CACHE_DIR}/miniconda3/miniconda.sh

${CACHE_DIR}/miniconda3/bin/conda init bash

set +x
source ~/.bashrc
set -x

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install --quiet --yes snakemake zstd wget mamba

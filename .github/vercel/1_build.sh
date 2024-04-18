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

mkdir -p ${WORK_DIR}/results
zstd --decompress ${WORK_DIR}/.github/data/data.tar.zst -c | tar xf - -C ${WORK_DIR}/results

snakemake --use-conda --cores 2

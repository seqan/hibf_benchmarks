#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -euvo pipefail

set +v && source ~/.bashrc && set -v
set +v && mamba activate snakemake && set -v

mkdir -p results
zstd --decompress .github/data/data.tar.zst -c | tar xf - -C results

snakemake --use-conda \
          --conda-frontend mamba \
          --conda-prefix `pwd`/node_modules/.snakemake \
          --conda-cleanup-pkgs cache \
          --cores 4

#!/usr/bin/env python
"""Generates CSV files with data timings from query.timings."""

## INPUT:
# 0         1       2                           3                           4                   5                                   6                           7                       8                                           9                                   10                                  11                          12                                          13                              14              15          16                  17                  18      19
# KEY	    SUBKEY	wall_clock_time_in_seconds	peak_memory_usage_in_KiB	index_size_in_KiB	determine_query_length_in_seconds	query_file_io_in_seconds	load_index_in_seconds	compute_minimiser_avg_per_thread_in_seconds	compute_minimiser_sum_in_seconds	query_ibf_avg_per_thread_in_seconds	query_ibf_sum_in_seconds	generate_results_avg_per_thread_in_seconds	generate_results_sum_in_seconds	time-v-seconds	time-v-time	time-v-memory-kbits	time-v-memory-human	command -
# k	        k=19	193.42	                    140326444	                136312729	        20.44	                            16.67	                    39.92	                2.70	                                    86.36	                            108.05	                            3457.46	                    3.05	                                    97.74		                    -               193.64	    3:13.64	            140326444	        133.83G "/srv/data/smehringer/raptor/build/bin/raptor search --index /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/k=19.index --query /srv/data/smehringer/RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq.only250.fastq --threshold 0.7 --output /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/k=19.result --timing-output /tmp/raptor.timings --threads 32"
# tmax=1024	none	769.75	                    212058284	                208038740           20.54	                            19.07	                    73.28	                2.85	                                    91.24	                            572.29	                            18313.29	                1.84	                                    58.99		                    -               770.01	    12:50.01	        212058284	        202.23G	"/srv/data/smehringer/raptor/build/bin/raptor search --index /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/tmax=1024-noU.index --query /srv/data/smehringer/RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq.only250.fastq --threshold 0.7 --output /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/tmax=1024-noU.result --timing-output /tmp/raptor.timings --threads 32"

## OUTPUT:
# [Subkeys]
# [Query File IO]
# [Determine Query Length]
# [Compute Minimiser Avg Per Thread]
# [Generate Results Avg Per Thread]
# [Load Index]
# [Query IBF Avg Per Thread (including time left))]

import csv
import os

from shared import input_timings, output_timings, write_csvs


def create_timing_csv():
    """Generates CSV files with data timings from query.timings."""
    keys = ["k", "alpha", "hash", "r-fpr", "none", "U+R", "U"]
    repeat = 3
    try:
        os.makedirs(output_timings, exist_ok=True)
        print(f"Directory {output_timings} created")
    except OSError:
        print(f"files replaced in directory {output_timings}")
    for key in keys:
        output_path = f"{output_timings}/{key}.csv"
        subkeys = []
        time_left = []
        query_ibf_avg_per_thread_in_seconds = []
        load_index_in_seconds = []
        generate_results_avg_per_thread_in_seconds = []
        compute_minimiser_avg_per_thread_in_seconds = []
        determine_query_length_in_seconds = []
        query_file_io_in_seconds = []
        with open(input_timings, "r", encoding="utf-8") as file:
            reader = csv.reader(file, delimiter="\t")
            i = 0
            wall_clock_time_in_seconds_i = 0
            query_ibf_avg_per_thread_in_seconds_i = 0
            load_index_in_seconds_i = 0
            generate_results_avg_per_thread_in_seconds_i = 0
            compute_minimiser_avg_per_thread_in_seconds_i = 0
            determine_query_length_in_seconds_i = 0
            query_file_io_in_seconds_i = 0
            for line in reader:
                if key in {line[0], line[1]}:
                    wall_clock_time_in_seconds_i += float(line[2])
                    query_ibf_avg_per_thread_in_seconds_i += float(line[10])
                    load_index_in_seconds_i += float(line[7])
                    generate_results_avg_per_thread_in_seconds_i += float(line[12])
                    compute_minimiser_avg_per_thread_in_seconds_i += float(line[8])
                    determine_query_length_in_seconds_i += float(line[5])
                    query_file_io_in_seconds_i += float(line[6])
                    i += 1
                    if i == repeat:
                        if line[1] == key:
                            subkeys.append(float(line[0].split("=")[1]))
                        else:
                            subkeys.append(float(line[1].split("=")[1]))
                        query_ibf_avg_per_thread_in_seconds.append(
                            round((query_ibf_avg_per_thread_in_seconds_i / repeat), 2)
                        )
                        load_index_in_seconds.append(round((load_index_in_seconds_i / repeat), 2))
                        generate_results_avg_per_thread_in_seconds.append(
                            round((generate_results_avg_per_thread_in_seconds_i / repeat), 2)
                        )
                        compute_minimiser_avg_per_thread_in_seconds.append(
                            round((compute_minimiser_avg_per_thread_in_seconds_i / repeat), 2)
                        )
                        determine_query_length_in_seconds.append(
                            round((determine_query_length_in_seconds_i / repeat), 2)
                        )
                        query_file_io_in_seconds.append(round((query_file_io_in_seconds_i / repeat), 2))
                        time_left.append(
                            round(
                                max(
                                    0,
                                    (
                                        wall_clock_time_in_seconds_i
                                        - query_ibf_avg_per_thread_in_seconds_i
                                        - load_index_in_seconds_i
                                        - generate_results_avg_per_thread_in_seconds_i
                                        - compute_minimiser_avg_per_thread_in_seconds_i
                                        - determine_query_length_in_seconds_i
                                        - query_file_io_in_seconds_i
                                    ),
                                )
                                / repeat,
                                2,
                            )
                        )
                        wall_clock_time_in_seconds_i = 0
                        query_ibf_avg_per_thread_in_seconds_i = 0
                        load_index_in_seconds_i = 0
                        generate_results_avg_per_thread_in_seconds_i = 0
                        compute_minimiser_avg_per_thread_in_seconds_i = 0
                        determine_query_length_in_seconds_i = 0
                        query_file_io_in_seconds_i = 0
                        i = 0
        (
            subkeys,
            time_left,
            query_ibf_avg_per_thread_in_seconds,
            load_index_in_seconds,
            generate_results_avg_per_thread_in_seconds,
            compute_minimiser_avg_per_thread_in_seconds,
            determine_query_length_in_seconds,
            query_file_io_in_seconds,
        ) = zip(
            *sorted(
                zip(
                    subkeys,
                    time_left,
                    query_ibf_avg_per_thread_in_seconds,
                    load_index_in_seconds,
                    generate_results_avg_per_thread_in_seconds,
                    compute_minimiser_avg_per_thread_in_seconds,
                    determine_query_length_in_seconds,
                    query_file_io_in_seconds,
                )
            )
        )
        query_ibf_avg_per_thread_in_seconds = [
            round(x + y, 2) for x, y in zip(query_ibf_avg_per_thread_in_seconds, time_left)
        ]
        all_lists = [
            subkeys,
            query_file_io_in_seconds,
            determine_query_length_in_seconds,
            compute_minimiser_avg_per_thread_in_seconds,
            generate_results_avg_per_thread_in_seconds,
            load_index_in_seconds,
            query_ibf_avg_per_thread_in_seconds,
        ]
        write_csvs(output_path, all_lists)


if __name__ == "__main__":
    create_timing_csv()

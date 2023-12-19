#!/usr/bin/env python
"""Shared variables."""

import csv
import os

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))
input_sizes = os.path.join(data_dir, "all.sizes")
input_timings = os.path.join(data_dir, "query.timings")
output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../build"))
output_sizes = os.path.join(output_dir, "sizes")
output_timings = os.path.join(output_dir, "timings")
html_dir = os.path.join(output_dir, "html")
html_file = os.path.join(html_dir, "index.html")


def write_csvs(output_path, all_lists):
    """Generates a CSV file for each list in all_lists."""
    with open(output_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        for element in all_lists:
            writer.writerow(element)

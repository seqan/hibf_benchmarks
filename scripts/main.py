#!/usr/bin/env python
"""Generates CSV files and creates plot."""

from get_data_sizes import create_size_csvs
from get_data_timings import create_timing_csv
from plot import create_plot


def main():
    """Main function."""
    create_size_csvs()
    create_timing_csv()
    create_plot()


if __name__ == "__main__":
    main()

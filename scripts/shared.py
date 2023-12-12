import os

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../data"))
input_sizes = os.path.join(data_dir, "all.sizes")
input_timings = os.path.join(data_dir, "query.timings")

output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../build"))
output_sizes = os.path.join(output_dir, "sizes")
output_timings = os.path.join(output_dir, "timings")

html_file = os.path.join(output_dir, "index.html")

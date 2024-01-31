# INPUTS
# peak_memory_usage_in_kibibytes	index_size_in_kibibytes	configured_threads	wall_clock_time_in_seconds	user_time_in_seconds	system_time_in_seconds	cpu_usage_in_percent	determine_query_length_in_seconds	complete_search_in_seconds	query_file_io_in_seconds	load_index_in_seconds	parallel_search_in_seconds	cpu_usage_parallel_search_in_percent	compute_minimiser_max_in_seconds	compute_minimiser_avg_in_seconds	query_ibf_max_in_seconds	query_ibf_avg_in_seconds	generate_results_max_in_seconds	generate_results_avg_in_seconds
# 1230196							884722					4					28.90						109.29					1.40					383.03					0.00								28.89						1.20						0.40					27.52						397.20									2.18								2.17								25.19						24.93						0.10							0.10

TIME_FORMAT <- c("peak_memory_usage_in_kibibytes",
	"index_size_in_kibibytes",
	"configured_threads",
	"wall_clock_time_in_seconds",
	"user_time_in_seconds",
	"system_time_in_seconds",
	"cpu_usage_in_percent",
	"determine_query_length_in_seconds",
	"complete_search_in_seconds",
	"query_file_io_in_seconds",
	"load_index_in_seconds",
	"parallel_search_in_seconds",
	"cpu_usage_parallel_search_in_percent",
	"compute_minimiser_max_in_seconds",
	"compute_minimiser_avg_in_seconds",
	"query_ibf_max_in_seconds",
	"query_ibf_avg_in_seconds",
	"generate_results_max_in_seconds",
	"generate_results_avg_in_seconds")

SIZE_FORMAT <- c("LEVEL",
	"BIT_SIZE",
	"IBFS",
	"AVG_LOAD_FACTOR",
	"TBS_TOO_BIG",
	"AVG_TBS_TOO_BIG_ELEMENTS",
	"AVG_MAX_ELEMENTS")

args <- commandArgs(trailingOnly = TRUE)

data_list <- list()

print(args[1])
SELECTED_FORMAT <- TIME_FORMAT
if (grepl("size", args[1])) {
	print("SIZE")
	SELECTED_FORMAT <- SIZE_FORMAT
}


for (file in args[2:length(args)]) {
    data <- read.delim(file, comment.char = "#", sep = "\t")
    data <- data[SELECTED_FORMAT]
    data_list[[file]] <- data
}

result_data <- do.call(rbind, data_list)
write.table(result_data, file = args[1], sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

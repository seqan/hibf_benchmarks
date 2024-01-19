# INPUTS
# peak_memory_usage_in_kibibytes	index_size_in_kibibytes	configured_threads	wall_clock_time_in_seconds	user_time_in_seconds	system_time_in_seconds	cpu_usage_in_percent	determine_query_length_in_seconds	complete_search_in_seconds	query_file_io_in_seconds	load_index_in_seconds	parallel_search_in_seconds	cpu_usage_parallel_search_in_percent	compute_minimiser_max_in_seconds	compute_minimiser_avg_in_seconds	query_ibf_max_in_seconds	query_ibf_avg_in_seconds	generate_results_max_in_seconds	generate_results_avg_in_seconds
# 1230196							884722					4					28.90						109.29					1.40					383.03					0.00								28.89						1.20						0.40					27.52						397.20									2.18								2.17								25.19						24.93						0.10							0.10

# OUTPUTS
# KEY		SUBKEY	wall_clock_time_in_seconds	peak_memory_usage_in_KiB	index_size_in_KiB	determine_query_length_in_seconds	query_file_io_in_seconds	load_index_in_seconds	compute_minimiser_avg_per_thread_in_seconds	compute_minimiser_sum_in_seconds	query_ibf_avg_per_thread_in_seconds	query_ibf_sum_in_seconds	generate_results_avg_per_thread_in_seconds	generate_results_sum_in_seconds	time-v-seconds	time-v-time	time-v-memory-kbits	time-v-memory-human	command
# k			k=19	193.42						140326444					136312729			20.44								16.67						39.92					2.70										86.36								108.05								3457.46						3.05										97.74							-				193.64		3:13.64				140326444	 		133.83G	"/srv/data/smehringer/raptor/build/bin/raptor search --index /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/k=19.index --query /srv/data/smehringer/RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq.only250.fastq --threshold 0.7 --output /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/k=19.result --timing-output /tmp/raptor.timings --threads 32"
# tmax=1024	none	769.75						212058284					208038740			20.54								19.07						73.28					2.85										91.24								572.29								18313.29					1.84										58.99							-				770.01		12:50.01			212058284	 		202.23G	"/srv/data/smehringer/raptor/build/bin/raptor search --index /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/tmax=1024-noU.index --query /srv/data/smehringer/RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq.only250.fastq --threshold 0.7 --output /srv/data/smehringer/refseq_parameter/relaxed_fpr_50/tmax=1024-noU.result --timing-output /tmp/raptor.timings --threads 32"

DIRECTORY <- "/srv/public/leonard/hibf_benchmarks/scripts/build"
FILE_NAME <- "out.time"
OUTPUT_FILE_NAME <- "out.time.csv"
FORMAT <- c("wall_clock_time_in_seconds", "peak_memory_usage_in_kibibytes", "index_size_in_kibibytes", "determine_query_length_in_seconds", "query_file_io_in_seconds", "load_index_in_seconds", "compute_minimiser_avg_in_seconds", "compute_minimiser_sum_in_seconds", "query_ibf_avg_in_seconds", "query_ibf_sum_in_seconds", "generate_results_avg_in_seconds", "generate_results_sum_in_seconds", "time-v-seconds", "time-v-time", "time-v-memory-kbits", "time-v-memory-human", "command")

result_data <- c(paste(FORMAT, collapse = "\t"))

sub_directories <- list.dirs(DIRECTORY, recursive = FALSE)
for (directory in sub_directories) {
	data <- read.delim(paste(directory, FILE_NAME, sep = "/"), header = TRUE, sep = "\t")

	wall_clock_time_in_seconds <- data$wall_clock_time_in_seconds
	peak_memory_usage_in_kibibytes <- data$peak_memory_usage_in_kibibytes
	index_size_in_kibibytes <- data$index_size_in_kibibytes
	determine_query_length_in_seconds <- data$determine_query_length_in_seconds
	query_file_io_in_seconds <- data$query_file_io_in_seconds
	load_index_in_seconds <- data$load_index_in_seconds
	compute_minimiser_avg_per_thread_in_seconds <- data$compute_minimiser_avg_in_seconds / data$configured_threads
	compute_minimiser_sum_in_seconds <- data$compute_minimiser_avg_in_seconds
	query_ibf_avg_per_thread_in_seconds <- data$query_ibf_avg_in_seconds / data$configured_threads
	query_ibf_sum_in_seconds <- data$query_ibf_avg_in_seconds
	generate_results_avg_per_thread_in_seconds <- data$generate_results_avg_in_seconds / data$configured_threads
	generate_results_sum_in_seconds <- data$generate_results_avg_in_seconds
	time_v_seconds <- data$wall_clock_time_in_seconds
	time_v_time <- data$wall_clock_time_in_seconds #???
	time_v_memory_kbits <- data$peak_memory_usage_in_kibibytes
	time_v_memory_human <- data$peak_memory_usage_in_kibibytes #???

	new_line <- paste(c(wall_clock_time_in_seconds,
		peak_memory_usage_in_kibibytes,
		index_size_in_kibibytes,
		determine_query_length_in_seconds,
		query_file_io_in_seconds,
		load_index_in_seconds,
		compute_minimiser_avg_per_thread_in_seconds,
		compute_minimiser_sum_in_seconds,
		query_ibf_avg_per_thread_in_seconds,
		query_ibf_sum_in_seconds,
		generate_results_avg_per_thread_in_seconds,
		generate_results_sum_in_seconds,
		time_v_seconds,
		time_v_time,
		time_v_memory_kbits,
		time_v_memory_human), collapse = "\t")
	result_data <- append(result_data, new_line)
}

writeLines(result_data, paste(DIRECTORY, OUTPUT_FILE_NAME, sep = "/"))

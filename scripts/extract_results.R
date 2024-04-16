log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

library(tidyr)

BUILD_DIR <- snakemake@params[["BUILD_DIR"]]
TIME_FORMAT <- snakemake@params[["TIME_FORMAT"]]
SIZE_FORMAT <- snakemake@params[["SIZE_FORMAT"]]
KEYS_FORMAT <- snakemake@params[["KEYS_FORMAT"]]
SIZE_INPUT <- snakemake@input[["SIZE_INPUT"]]
TIME_INPUT <- snakemake@input[["TIME_INPUT"]]


transform_size <- function(df) {
    df$BIT_SIZE <- df$BIT_SIZE / 8 / 1000^3
    df$AVG_LOAD_FACTOR <- df$AVG_LOAD_FACTOR / 100
    df_vollständig <- expand.grid(SUBKEY = unique(df$SUBKEY), LEVEL = 0:3)
    df_vollständig <- merge(df_vollständig, df, by = c("SUBKEY", "LEVEL"), all.x = TRUE)

    df_vollständig[is.na(df_vollständig)] <- 0
    df_wide <- pivot_wider(df_vollständig, names_from = LEVEL, values_from = c(BIT_SIZE, AVG_LOAD_FACTOR), names_sep = "_")
    df_wide <- data.frame(t(df_wide))
    return(df_wide)
}

transform_time <- function(df) {
    df <- t(df)
    return(df)
}

input_time <- read.delim(TIME_INPUT, comment.char = "#", sep = "\t", header = TRUE)
input_size <- read.delim(SIZE_INPUT, comment.char = "#", sep = "\t", header = TRUE)

new_subkeys <- strsplit(input_time$SUBKEY, "=")
input_time$SUBKEY <- sapply(new_subkeys, function(x) x[2])
new_subkeys <- strsplit(input_size$SUBKEY, "=")
input_size$SUBKEY <- sapply(new_subkeys, function(x) x[2])

for (key in KEYS_FORMAT) {
    time_data <- input_time[input_time$KEY == key, TIME_FORMAT[TIME_FORMAT != "KEY" & TIME_FORMAT != "wall_clock_time_in_seconds"]]
    time_data <- transform_time(time_data)
    output_time <- file.path(BUILD_DIR, "prepared_time", key)
    write.table(time_data, file = output_time, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    size_data <- input_size[input_size$KEY == key, SIZE_FORMAT[SIZE_FORMAT != "KEY"]]
    size_data <- transform_size(size_data)
    output_size <- file.path(BUILD_DIR, "prepared_size", key)
    write.table(size_data, file = output_size, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
print("Results have been extracted.")

sink(NULL)
sink(NULL, type = "message")

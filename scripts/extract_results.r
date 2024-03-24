# Output
# prepared_time/alpha:
# [Subkeys]
# [Query File IO]
# [Determine Query Length]
# [Compute Minimiser Avg Per Thread]
# [Generate Results Avg Per Thread]
# [Load Index]
# [Query IBF Avg Per Thread]

# prepared_size/alpha:
# [Subkeys]
# [Level 0]
# [Level 1]
# [Level 2]
# [Level 3]
# [Avg Load Factor]

# KEYS <- c("kmer", "alpha", "hash", "relaxed-fpr", "none", "U+R", "U")

library(tidyr)

config <- yaml::read_yaml("/srv/public/leonard/hibf_benchmarks/scripts/parameters.yaml")
TIME_FORMAT <- config$TIME_FORMAT
SIZE_FORMAT <- config$SIZE_FORMAT
KEYS <- config$KEYS

args <- commandArgs(trailingOnly = TRUE)

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

input_time <- read.delim(args[1], comment.char = "#", sep = "\t", header = TRUE)
input_size <- read.delim(args[2], comment.char = "#", sep = "\t", header = TRUE)

new_subkeys <- strsplit(input_time$SUBKEY, "=")
input_time$SUBKEY <- sapply(new_subkeys, function(x) x[2])
new_subkeys <- strsplit(input_size$SUBKEY, "=")
input_size$SUBKEY <- sapply(new_subkeys, function(x) x[2])

output_time <- file.path(config$BUILD_DIR, "prepared_time")
if (!dir.exists(output_time)) dir.create(output_time, recursive = TRUE)
output_size <- file.path(config$BUILD_DIR, "prepared_size")
if (!dir.exists(output_size)) dir.create(output_size, recursive = TRUE)

for (key in KEYS) {
    time_data <- input_time[input_time$KEY == key, TIME_FORMAT[TIME_FORMAT != "KEY" & TIME_FORMAT != "wall_clock_time_in_seconds"]]
    time_data <- transform_time(time_data)
    output_time <- file.path(config$BUILD_DIR, "prepared_time", key)
    write.table(time_data, file = output_time, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    size_data <- input_size[input_size$KEY == key, SIZE_FORMAT[SIZE_FORMAT != "KEY"]]
    size_data <- transform_size(size_data)
    output_size <- file.path(config$BUILD_DIR, "prepared_size", key)
    write.table(size_data, file = output_size, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

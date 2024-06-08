log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

INPUT_FILES <- snakemake@input[["INPUT_FILES"]]
OUTPUT_FILE <- snakemake@output[["OUTPUT_FILE"]]

extract_part <- function(file) {
    # Extract the key and subkey from the file name.
    parts <- strsplit(sub("^.*/(.+)/(.+)", "\\1", file), "=")[[1]]
    key <- parts[1]
    value <- as.numeric(gsub("_", ".", parts[2], fixed = TRUE))
    subkey <- paste(key, value, sep = "=")
    return(list(key, subkey))
}

get_key <- function(key) {
    # Get the key from the subkey.
    if (startsWith(key, "none=") || startsWith(key, "U=") || startsWith(key, "U+R=")) {
        params <- strsplit(key, "=")[[1]]
        key <- paste("tmax", params[2], sep = "=")
    }
    return(key)
}

read_and_modify <- function(FILE_PATH) {
    # Read the file and add the key and subkey.
    data <- read.delim(FILE_PATH, comment.char = "#", sep = "\t", header = TRUE)
    extracted_part <- extract_part(FILE_PATH)
    data <- data.frame(KEY = extracted_part[[1]], SUBKEY = get_key(extracted_part[[2]]), data)
    return(data)
}

list_of_dfs <- lapply(INPUT_FILES, read_and_modify)
combined_df <- do.call(rbind, list_of_dfs)

write.table(combined_df, file = OUTPUT_FILE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
print("Results have been summarized.")

sink(NULL)
sink(NULL, type = "message")

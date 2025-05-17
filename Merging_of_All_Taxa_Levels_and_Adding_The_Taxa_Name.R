# Load necessary library
library(dplyr)

# Set the directory containing your files
input_dir <- "your/path/to/miobiogen/exposure"
output_dir <- "your/path/to/miobiogen/exposure_final"  # You can use input_dir if you want to overwrite

# Ensure output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir)

# List all relevant files in the directory
file_list <- list.files(input_dir, full.names = TRUE, pattern = "^exposure_data_.*\\.csv$")

# Function to extract taxonomy level and name from filename
extract_taxonomy <- function(file_name) {
  # Identify taxonomy type and extract corresponding name
  taxonomy_levels <- c("genus", "order", "phyla", "class", "family")
  
  for (level in taxonomy_levels) {
    if (grepl(paste0("exposure_data_", level, "_"), file_name)) {
      tax_name <- sub(paste0("exposure_data_", level, "_([^\\.]+).*"), "\\1", file_name)
      return(list(level = level, name = tax_name))
    }
  }
  return(NULL)  # Return NULL if no match
}

# Function to process each file
process_file <- function(file_path) {
  # Extract filename without directory
  file_name <- basename(file_path)
  
  # Extract taxonomy info
  tax_info <- extract_taxonomy(file_name)
  
  if (is.null(tax_info)) {
    print(paste("Skipping file:", file_name, "- No taxonomy match found."))
    return(NULL)
  }
  
  # Read the data (adjust sep if TSV)
  data <- read.csv(file_path, sep = ",", stringsAsFactors = FALSE)
  
  # Add the "bac" column with extracted taxonomy name
  data$bac <- tax_info$name
  data$taxonomy_level <- tax_info$level  # Optional: store the level as well
  
  # Save the modified file
  output_path <- file.path(output_dir, file_name)
  write.csv(data, output_path, row.names = FALSE)
  
  return(output_path)
}

# Process all files
updated_files <- lapply(file_list, process_file)

print("Processing complete. Updated files saved in the output directory.")
#######################################################################################################

# Set the directory path containing CSV files
input_dir <- "your/path/to/miobiogen/exposure_final"   # Change this to your folder path

# List all CSV files in the directory
csv_files <- list.files(path = input_dir, pattern = "*.csv", full.names = TRUE)

# Read and merge all CSV files
merged_data <- do.call(rbind, lapply(csv_files, read.csv, stringsAsFactors = FALSE))

# Save the merged data to a new CSV file
write.csv(merged_data, "merged_output.csv", row.names = FALSE)

print("Merging complete! Check merged_output.csv")
################################################################################################

# Load your data (assuming it's a CSV file)
df <- read.csv("merged_output.csv", stringsAsFactors = FALSE)

# Create a new column with the desired format
df <- df %>%
  mutate(taxonomy_bac = paste0(taxonomy_level, ".", bac))

# View the updated data
head(df)

# Optionally, if you want to save it back
write.csv(df, "mibiogen_exposure.csv", row.names = FALSE)


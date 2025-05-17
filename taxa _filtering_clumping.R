#########################################LOADING PACKAGES################################################

library(data.table)
library(dplyr)
library(ieugwasr)
library(remotes)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(tibble)  # Contains as_tibble() function
library(tidyverse)
library(LDlinkR)
library(magrittr)
library(gt)
library(stringr)

######################################genus########################################################
output_dir <- "your/path/to/miobiogen"

# List all TSV files in your folder
file_list <- list.files(path = "your/path/to/miobiogen/genus", pattern = "*.txt", full.names = TRUE)
total_files <- length(file_list)

# Set your p-value threshold
p_threshold <- 1e-5

# Process each TSV file one by one
for (i in 1:131) {
  file <- file_list[i]
  
  # Extract only the relevant file name
  raw_name <- basename(file)  # Get the filename without path
  clean_name <- str_replace_all(raw_name, 
                                pattern = "genus\\.|\\.id\\.|\\.summary|\\.txt", 
                                replacement = "") %>% 
    str_replace_all(" ", "_")  # Replace spaces with underscores if any
  
  # Read one TSV file
  dt <- fread(file, sep = "\t")
  
  # Filter based on p-value threshold
  filtered_data_genus <- dt %>% filter(P.weightedSumZ < p_threshold)
  
  # Remove duplicate SNPs
  exposure_unique <- filtered_data_genus %>% distinct(rsID, .keep_all = TRUE)
  
  # Convert to data.frame
  exposure_unique <- data.frame(exposure_unique)
  
  # Format data for TwoSampleMR
  exposure_formatted <- format_data(
    exposure_unique,
    type = "exposure",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "eff.allele",
    other_allele_col = "ref.allele",
    pval_col = "P.weightedSumZ"
  )
  
  Sys.setenv(R_ENVIRON_USER = normalizePath("~/.Renviron"))
  Sys.getenv("R_ENVIRON_USER")
  
  
  #options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
  clumped_dat1 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat2 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "SAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat3 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat4 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AFR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat5 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AMR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat6 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "legacy",
    bfile = NULL,
    plink_bin = NULL
  )
  
  clumped_data <- bind_rows(
    clumped_dat1,
    clumped_dat2,
    clumped_dat3,
    clumped_dat4,
    clumped_dat5,
    clumped_dat6
  ) %>% 
    distinct(SNP, .keep_all = TRUE)
  
  # Remove palindromic SNPs
  is_palindromic <- function(a1, a2) {
    return((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
             (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  }
  
  exposure_data <- clumped_data %>%
    rowwise() %>%
    filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))
  
  # Save results with cleaned file name
  output_path <- file.path(output_dir, paste0("exposure_data_genus_", clean_name, ".csv"))
  
  # Save the file
  write.csv(exposure_data, output_path, row.names = FALSE)
  print(paste("Saved:", output_path))
  
  # Cleanup
  rm()
  rm()
  gc()
  gc()
}

print("Processing complete genus   HWAAAAHHHHHAHAHHAHAH!")



#################################order###################################################################


output_dir <- "E:/try/MS_MIBIOGEN/exposure"

# List all TSV files in your folder
file_list <- list.files(path = "your/path/to/miobiogen/order", pattern = "*.txt", full.names = TRUE)
total_files <- length(file_list)

# Set your p-value threshold
p_threshold <- 1e-5

# Process each TSV file one by one
for (i in 1:131) {
  file <- file_list[i]
  
  # Extract only the relevant file name
  raw_name <- basename(file)  # Get the filename without path
  clean_name <- str_replace_all(raw_name, 
                                pattern = "genus\\.|\\.id\\.|\\.summary|\\.txt", 
                                replacement = "") %>% 
    str_replace_all(" ", "_")  # Replace spaces with underscores if any
  
  # Read one TSV file
  dt <- fread(file, sep = "\t")
  
  # Filter based on p-value threshold
  filtered_data_genus <- dt %>% filter(P.weightedSumZ < p_threshold)
  
  # Remove duplicate SNPs
  exposure_unique <- filtered_data_genus %>% distinct(rsID, .keep_all = TRUE)
  
  # Convert to data.frame
  exposure_unique <- data.frame(exposure_unique)
  
  # Format data for TwoSampleMR
  exposure_formatted <- format_data(
    exposure_unique,
    type = "exposure",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "eff.allele",
    other_allele_col = "ref.allele",
    pval_col = "P.weightedSumZ"
  )
  
  Sys.setenv(R_ENVIRON_USER = normalizePath("~/.Renviron"))
  Sys.getenv("R_ENVIRON_USER")
  
  
  #options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
  clumped_dat1 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat2 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "SAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat3 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat4 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AFR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat5 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AMR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat6 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "legacy",
    bfile = NULL,
    plink_bin = NULL
  )
  
  clumped_data <- bind_rows(
    clumped_dat1,
    clumped_dat2,
    clumped_dat3,
    clumped_dat4,
    clumped_dat5,
    clumped_dat6
  ) %>% 
    distinct(SNP, .keep_all = TRUE)
  
  # Remove palindromic SNPs
  is_palindromic <- function(a1, a2) {
    return((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
             (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  }
  
  exposure_data <- clumped_data %>%
    rowwise() %>%
    filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))
  
  # Save results with cleaned file name
  output_path <- file.path(output_dir, paste0("exposure_data_order_", clean_name, ".csv"))
  
  # Save the file
  write.csv(exposure_data, output_path, row.names = FALSE)
  print(paste("Saved:", output_path))
  
  # Cleanup
  rm()
  rm()
  gc()
  gc()
}

print("Processing complete order   HWAAAAHHHHHAHAHHAHAH!")



#################################phyla###################################################################


output_dir <- "E:/try/MS_MIBIOGEN/exposure"

# List all TSV files in your folder
file_list <- list.files(path = "your/path/to/miobiogen/phyla", pattern = "*.txt", full.names = TRUE)
total_files <- length(file_list)

# Set your p-value threshold
p_threshold <- 1e-5

# Process each TSV file one by one
for (i in 1:131) {
  file <- file_list[i]
  
  # Extract only the relevant file name
  raw_name <- basename(file)  # Get the filename without path
  clean_name <- str_replace_all(raw_name, 
                                pattern = "phylum\\.|\\.id\\.|\\.summary|\\.txt", 
                                replacement = "") %>% 
    str_replace_all(" ", "_")  # Replace spaces with underscores if any
  
  # Read one TSV file
  dt <- fread(file, sep = "\t")
  
  # Filter based on p-value threshold
  filtered_data_genus <- dt %>% filter(P.weightedSumZ < p_threshold)
  
  # Remove duplicate SNPs
  exposure_unique <- filtered_data_genus %>% distinct(rsID, .keep_all = TRUE)
  
  # Convert to data.frame
  exposure_unique <- data.frame(exposure_unique)
  
  # Format data for TwoSampleMR
  exposure_formatted <- format_data(
    exposure_unique,
    type = "exposure",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "eff.allele",
    other_allele_col = "ref.allele",
    pval_col = "P.weightedSumZ"
  )
  
  Sys.setenv(R_ENVIRON_USER = normalizePath("~/.Renviron"))
  Sys.getenv("R_ENVIRON_USER")
  
  
  #options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
  clumped_dat1 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat2 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "SAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat3 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat4 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AFR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat5 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AMR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat6 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "legacy",
    bfile = NULL,
    plink_bin = NULL
  )
  
  clumped_data <- bind_rows(
    clumped_dat1,
    clumped_dat2,
    clumped_dat3,
    clumped_dat4,
    clumped_dat5,
    clumped_dat6
  ) %>% 
    distinct(SNP, .keep_all = TRUE)
  
  # Remove palindromic SNPs
  is_palindromic <- function(a1, a2) {
    return((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
             (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  }
  
  exposure_data <- clumped_data %>%
    rowwise() %>%
    filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))
  
  # Save results with cleaned file name
  output_path <- file.path(output_dir, paste0("exposure_data_phyla_", clean_name, ".csv"))
  
  # Save the file
  write.csv(exposure_data, output_path, row.names = FALSE)
  print(paste("Saved:", output_path))
  
  # Cleanup
  rm()
  rm()
  gc()
  gc()
}

print("Processing complete phyla HWAAAAHHHHHAHAHHAHAH!")

#################################class###################################################################


output_dir <- "E:/try/MS_MIBIOGEN/exposure"

# List all TSV files in your folder
file_list <- list.files(path = "your/path/to/miobiogen/class", pattern = "*.txt", full.names = TRUE)
total_files <- length(file_list)

# Set your p-value threshold
p_threshold <- 1e-5

# Process each TSV file one by one
for (i in 1:131) {
  file <- file_list[i]
  
  # Extract only the relevant file name
  raw_name <- basename(file)  # Get the filename without path
  clean_name <- str_replace_all(raw_name, 
                                pattern = "class\\.|\\.id\\.|\\.summary|\\.txt", 
                                replacement = "") %>% 
    str_replace_all(" ", "_")  # Replace spaces with underscores if any
  
  # Read one TSV file
  dt <- fread(file, sep = "\t")
  
  # Filter based on p-value threshold
  filtered_data_genus <- dt %>% filter(P.weightedSumZ < p_threshold)
  
  # Remove duplicate SNPs
  exposure_unique <- filtered_data_genus %>% distinct(rsID, .keep_all = TRUE)
  
  # Convert to data.frame
  exposure_unique <- data.frame(exposure_unique)
  
  # Format data for TwoSampleMR
  exposure_formatted <- format_data(
    exposure_unique,
    type = "exposure",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "eff.allele",
    other_allele_col = "ref.allele",
    pval_col = "P.weightedSumZ"
  )
  
  Sys.setenv(R_ENVIRON_USER = normalizePath("~/.Renviron"))
  Sys.getenv("R_ENVIRON_USER")
  
  
  #options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
  clumped_dat1 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat2 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "SAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat3 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat4 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AFR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat5 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AMR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat6 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "legacy",
    bfile = NULL,
    plink_bin = NULL
  )
  
  clumped_data <- bind_rows(
    clumped_dat1,
    clumped_dat2,
    clumped_dat3,
    clumped_dat4,
    clumped_dat5,
    clumped_dat6
  ) %>% 
    distinct(SNP, .keep_all = TRUE)
  
  # Remove palindromic SNPs
  is_palindromic <- function(a1, a2) {
    return((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
             (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  }
  
  exposure_data <- clumped_data %>%
    rowwise() %>%
    filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))
  
  # Save results with cleaned file name
  output_path <- file.path(output_dir, paste0("exposure_data_class_", clean_name, ".csv"))
  
  # Save the file
  write.csv(exposure_data, output_path, row.names = FALSE)
  print(paste("Saved:", output_path))
  
  # Cleanup
  rm()
  rm()
  gc()
  gc()
}

print("Processing complete class HWAAAAHHHHHAHAHHAHAH!")


#################################family###################################################################


output_dir <- "E:/try/MS_MIBIOGEN/exposure"

# List all TSV files in your folder
file_list <- list.files(path = "your/path/to/miobiogen/family", pattern = "*.txt", full.names = TRUE)
total_files <- length(file_list)

# Set your p-value threshold
p_threshold <- 1e-5

# Process each TSV file one by one
for (i in 1:131) {
  file <- file_list[i]
  
  # Extract only the relevant file name
  raw_name <- basename(file)  # Get the filename without path
  clean_name <- str_replace_all(raw_name, 
                                pattern = "family\\.|\\.id\\.|\\.summary|\\.txt", 
                                replacement = "") %>% 
    str_replace_all(" ", "_")  # Replace spaces with underscores if any
  
  # Read one TSV file
  dt <- fread(file, sep = "\t")
  
  # Filter based on p-value threshold
  filtered_data_genus <- dt %>% filter(P.weightedSumZ < p_threshold)
  
  # Remove duplicate SNPs
  exposure_unique <- filtered_data_genus %>% distinct(rsID, .keep_all = TRUE)
  
  # Convert to data.frame
  exposure_unique <- data.frame(exposure_unique)
  
  # Format data for TwoSampleMR
  exposure_formatted <- format_data(
    exposure_unique,
    type = "exposure",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "eff.allele",
    other_allele_col = "ref.allele",
    pval_col = "P.weightedSumZ"
  )
  
  Sys.setenv(R_ENVIRON_USER = normalizePath("~/.Renviron"))
  Sys.getenv("R_ENVIRON_USER")
  
  
  #options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
  clumped_dat1 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat2 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "SAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat3 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EAS",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat4 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AFR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat5 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "AMR",
    bfile = NULL,
    plink_bin = NULL
  )
  clumped_dat6 <- clump_data(
    exposure_formatted,
    clump_kb = 1000,
    clump_r2 = 0.01,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "legacy",
    bfile = NULL,
    plink_bin = NULL
  )
  
  clumped_data <- bind_rows(
    clumped_dat1,
    clumped_dat2,
    clumped_dat3,
    clumped_dat4,
    clumped_dat5,
    clumped_dat6
  ) %>% 
    distinct(SNP, .keep_all = TRUE)
  
  # Remove palindromic SNPs
  is_palindromic <- function(a1, a2) {
    return((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
             (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
  }
  
  exposure_data <- clumped_data %>%
    rowwise() %>%
    filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))
  
  # Save results with cleaned file name
  output_path <- file.path(output_dir, paste0("exposure_data_family_", clean_name, ".csv"))
  
  # Save the file
  write.csv(exposure_data, output_path, row.names = FALSE)
  print(paste("Saved:", output_path))
  
  # Cleanup
  rm()
  rm()
  gc()
  gc()
}

print("Processing complete family HWAAAAHHHHHAHAHHAHAH!")

########################################################################################################



print("ALL DONE YAYAYAYAYAYAY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

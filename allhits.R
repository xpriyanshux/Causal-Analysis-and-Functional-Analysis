##########################packages###################################################
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

########################pre-requisites:dataformat and munge proxies functiopn definition##########
# Define column types for summary statistics
coltypes = cols(
  ID = col_character(),
  CHROM = col_double(),
  POS = col_double(),
  REF = col_character(),
  ALT = col_character(),
  AF = col_double(),
  TRAIT = col_character(),
  BETA = col_double(),
  SE = col_double(),
  Z = col_double(),
  P = col_double(),
  N = col_double(),
  OR = col_double(),
  OR_L95 = col_double(),
  OR_U95 = col_double(),
  DIR = col_character(),
  G1000_ID = col_character(),
  G1000_VARIANT = col_character(),
  DBSNP_ID = col_character(),
  DBSNP_VARIANT = col_character(),
  OLD_ID = col_character(),
  OLD_VARIANT = col_character()
)


munge_proxies <- function(LDLink_file, outcome, outcome_clump){
  LDLink_file_path <- LDLink_file
  proxy_snps <- read_tsv(LDLink_file_path, skip = 1, col_names = F) %>%
    rename(id = X1, func = X2, proxy_snp = X3, coord = X4, alleles = X5, maf = X6, 
           distance = X7, dprime = X8, rsq = X9, correlated_alleles = X10, FORGEdb = X11, RegulomeDB = X12) %>%
    separate(coord, c('chr', 'pos'), sep = ":") %>%
    mutate(snp = ifelse(id == 1, proxy_snp, NA), 
           chr = str_replace(chr, 'chr', ""), 
           chr = as.numeric(chr), 
           pos = as.numeric(pos)) %>%
    fill(snp, .direction = 'down') %>%
    relocate(snp, .before = proxy_snp) %>%
    dplyr::select(-id, -func, -FORGEdb, -RegulomeDB) %>%
    filter(rsq >= 0.8)
  
  # Munge proxy SNP and outcome data
  proxy_outcome <- left_join(
    proxy_snps, outcome, by = c("proxy_snp" = "SNP")
  ) %>%
    separate(correlated_alleles, c("target_a1.outcome", "proxy_a1.outcome", 
                                   "target_a2.outcome", "proxy_a2.outcome"), sep = ",|=") %>%
    filter(!is.na(chr.outcome)) %>%
    arrange(snp, -rsq, abs(distance)) %>%
    group_by(snp) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      proxy.outcome = TRUE,
      target_snp.outcome = snp,
      proxy_snp.outcome = proxy_snp, 
    ) %>% 
    mutate(
      new_effect_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a1.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a2.outcome,
        TRUE ~ NA_character_
      ), 
      new_other_allele.outcome = case_when(
        proxy_a1.outcome == effect_allele.outcome & proxy_a2.outcome == other_allele.outcome ~ target_a2.outcome,
        proxy_a2.outcome == effect_allele.outcome & proxy_a1.outcome == other_allele.outcome ~ target_a1.outcome,
        TRUE ~ NA_character_
      ), 
      effect_allele.outcome = new_effect_allele.outcome, 
      other_allele.outcome = new_other_allele.outcome
    ) %>%
    dplyr::select(-proxy_snp, -chr, -pos, -alleles, -maf, -distance, -rsq, -dprime,  
                  -new_effect_allele.outcome, -new_other_allele.outcome) %>%
    relocate(target_a1.outcome, proxy_a1.outcome, target_a2.outcome, proxy_a2.outcome, .after = proxy_snp.outcome) %>%
    rename(SNP = snp) 
  
  # Ensure `samplesize.outcome` exists in both datasets
  outcome_clump <- outcome_clump %>%
    mutate(samplesize.outcome = ifelse("samplesize.outcome" %in% colnames(outcome_clump), samplesize.outcome, NA))
  
  proxy_outcome <- proxy_outcome %>%
    mutate(samplesize.outcome = NA) %>%  # Adding the missing column
    relocate(SNP, .after = samplesize.outcome)  # Now this line will not cause an error
  
  # Merge outcome and proxy outcomes
  outcome_dat <- bind_rows(
    outcome_clump, proxy_outcome
  ) %>% 
    arrange(chr.outcome, pos.outcome)
  
  return(outcome_dat)
}


##########exposure_formatting and ld clumping#############################################

exposure_path = "D:/final project/project work/exposure processing/mibiogen_exposure.csv"
exposure_ss <- fread(exposure_path, header = TRUE, sep = ",")

exposure_formatted <- exposure_ss %>%
  rename(
    SNP = SNP,
    beta = beta.exposure,
    se = se.exposure,
    effect_allele = effect_allele.exposure,
    other_allele = other_allele.exposure,
    pval = pval.exposure,
    eaf = eaf.exposure, # If available
    samplesize = id.exposure # This might need correction
  ) %>%
  select(SNP, beta, se, effect_allele, other_allele, pval, eaf, samplesize) %>%
  as_tibble()


# Save formatted exposure data
fwrite(exposure_formatted, "formatted_exposure.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Load formatted data into TwoSampleMR
exposure_data <- read_exposure_data(
  filename = "formatted_exposure.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "samplesize"
)

# Check the first few rows
head(exposure_data)

Sys.setenv(R_ENVIRON_USER = normalizePath("~/.Renviron"))
Sys.getenv("R_ENVIRON_USER")


#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
clumped_dat1 <- clump_data(
  exposure_data,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR",
  bfile = NULL,
  plink_bin = NULL
)
clumped_dat2 <- clump_data(
  exposure_data,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "SAS",
  bfile = NULL,
  plink_bin = NULL
)
clumped_dat3 <- clump_data(
  exposure_data,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EAS",
  bfile = NULL,
  plink_bin = NULL
)
clumped_dat4 <- clump_data(
  exposure_data,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "AFR",
  bfile = NULL,
  plink_bin = NULL
)
clumped_dat5 <- clump_data(
  exposure_data,
  clump_kb = 1000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "AMR",
  bfile = NULL,
  plink_bin = NULL
)
clumped_dat6 <- clump_data(
  exposure_data,
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
)

exposure_clump <- clumped_data
exposure_dat <- exposure_clump

##############################################outcome formatting#####################################################



#Define the coltypes object
coltypes <- cols(
  rsids = col_character(),
  `#chrom` = col_integer(),
  pos = col_integer(),
  ref = col_character(),
  alt = col_character(),
  af_alt = col_double(),
  beta = col_double(),
  sebeta = col_double(),
  pval = col_double(),
  nearest_genes = col_character()
)
outcome_path = "E:/try/summary_stats_release_finngen_R12_G6_ALS"
outcome_ss <- read_tsv(outcome_path, comment = "##", col_types = coltypes, 
                       col_select = c(rsids, `#chrom`, pos, ref, alt, af_alt, beta, sebeta, pval, nearest_genes))

# Format outcome
outcome <- outcome_ss %>%
  format_data(.,
              type = "outcome",
              snps = NULL,
              header = TRUE,
              phenotype_col = "nearest_genes",
              snp_col = "rsids",
              beta_col = "beta",
              se_col = "sebeta",
              eaf_col = "af_alt",
              effect_allele_col = "alt",
              other_allele_col = "ref",
              pval_col = "pval",
              chr_col = "#chrom",
              pos_col = "pos",
              log_pval = FALSE
  ) %>%
  as_tibble()

# extract exposure SNPs present in outcome
outcome_clump <- semi_join(
  outcome, exposure_dat, by = "SNP"
)

# Exposure SNPs not present in outomce
exp_snps_wo <- anti_join(
  exposure_dat, outcome, by = "SNP"
)
df1 <- exp_snps_wo[1:50, ]    # First 50 rows  
df2 <- exp_snps_wo[51:100, ]  # Next 50 rows  
df3 <- exp_snps_wo[101:150, ] # Next 50 rows  
df4 <- exp_snps_wo[151:200, ] # Next 50 rows  
df5 <- exp_snps_wo[201:250, ] # Next 50 rows  
df6 <- exp_snps_wo[251:252, ] # Next 50 rows 

# Remove rows where SNP column has "."

df3 <- df3[df3$SNP != ".", ]
#df5 <- df5[df5$SNP != ".", ]

#Sys.sleep(30)

# Use LDLinkR to identify proxy snps
LDproxy_batch(df1$SNP, 
              pop = "EUR",             
              r2d = "r2", 
              token = '182bf70ec476', 
              append = TRUE,           
              genome_build = "grch37") 
system("mv combined_query_snp_list_grch37.txt data/exposure_outcome_proxy_snps_1.txt")
Sys.sleep(10)

LDproxy_batch(df2$SNP, 
              pop = "ALL",             
              r2d = "r2", 
              token = '182bf70ec476', 
              append = TRUE,           
              genome_build = "grch37") 
system("mv combined_query_snp_list_grch37.txt data/exposure_outcome_proxy_snps_2.txt")
Sys.sleep(60)

LDproxy_batch(df3$SNP, 
              pop = "ALL",             
              r2d = "r2", 
              token = '182bf70ec476', 
              append = TRUE,           
              genome_build = "grch37") 
system("mv combined_query_snp_list_grch37.txt data/exposure_outcome_proxy_snps_3.txt")
Sys.sleep(10)

LDproxy_batch(df4$SNP, 
              pop = "ALL",             
              r2d = "r2", 
              token = '182bf70ec476', 
              append = TRUE,           
              genome_build = "grch37") 
system("mv combined_query_snp_list_grch37.txt data/exposure_outcome_proxy_snps_4.txt")
Sys.sleep(10)

LDproxy_batch(df5$SNP, 
              pop = "ALL",             
              r2d = "r2", 
              token = '182bf70ec476', 
              append = TRUE,           
              genome_build = "grch37") 
system("mv combined_query_snp_list_grch37.txt data/exposure_outcome_proxy_snps_5.txt")
Sys.sleep(10)

LDproxy_batch(df6$SNP, 
              pop = "ALL",             
              r2d = "r2", 
              token = '182bf70ec476', 
              append = TRUE,           
              genome_build = "grch37") 
system("mv combined_query_snp_list_grch37.txt data/exposure_outcome_proxy_snps_6.txt")
Sys.sleep(10)

########################merging the proxy snps#####################################

# Define the folder path
folder_path <- "D:/final project/project work/multiple sclerosis/data"

# List all matching files inside the "data" folder
files <- list.files(path = folder_path, pattern = "exposure_outcome_proxy_snps_.*\\.txt", full.names = TRUE)

# Read and merge all files
merged_data <- do.call(c, lapply(files, readLines))

# Write the merged content to a new file
writeLines(merged_data, "merged_snps.txt")

cat("Merging complete! File saved as merged_snps.txt\n")


# Munge proxy snp file
outcome_dat <- munge_proxies("merged_snps.txt", outcome, outcome_clump)

exposure_dat <- exposure_dat[exposure_dat$SNP != ".", ]

#length(intersect(exposure_dat$SNP, outcome_dat$SNP))
#sum(outcome_dat$SNP %in% exposure_dat$SNP)
#outcome_dat <- outcome_dat[!duplicated(outcome_dat$SNP), ]



is_palindromic <- function(a1, a2) {
  return((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
           (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
}

exposure_dat <- exposure_dat %>%
  rowwise() %>%
  filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))

outcome_dat <- outcome_dat %>%
  rowwise() %>%
  filter(!is_palindromic(effect_allele.outcome,other_allele.outcome))

#colnames(exposure_dat)[colnames(exposure_dat) == "beta"] <- "beta.exposure"
#colnames(exposure_dat)[colnames(exposure_dat) == "se"] <- "se.exposure"
#colnames(exposure_dat)[colnames(exposure_dat) == "effect_allele"] <- "effect_allele.exposure"
#colnames(exposure_dat)[colnames(exposure_dat) == "other_allele"] <- "other_allele.exposure"
#colnames(exposure_dat)[colnames(exposure_dat) == "eaf"] <- "eaf.exposure"

# Now check again
#colnames(exposure_dat)
exposure_dat$exposure <- "gut_microbiome"  # Change this label as needed



mr_dat <- harmonise_data(exposure_dat, outcome_dat)

# Run all desired MR methods
methods_to_run <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")

# Run MR
all_res <- mr(mr_dat, method_list = methods_to_run)

# Sensitivity analyses
all_het <- mr_heterogeneity(mr_dat)
all_pleio <- mr_pleiotropy_test(mr_dat)

# Extract unique SNPs for later
rsid_df <- mr_dat %>%
  select(SNP, id.exposure, id.outcome) %>%
  distinct()

# Loop through each method
for (method_name in unique(all_res$method)) {
  # Step 1: Filter results for this method
  method_res <- all_res %>% filter(method == method_name)
  
  # Step 2: Filter by nsnp <= 7
  filtered_res <- method_res %>% filter(nsnp <= 7)
  
  # Step 3: Get exposure-outcome pairs
  pairs <- filtered_res %>% select(id.exposure, id.outcome)
  
  # Step 4: Subset original harmonised data
  filtered_mr_dat <- mr_dat %>%
    semi_join(pairs, by = c("id.exposure", "id.outcome"))
  
  # Step 5: Get heterogeneity and pleiotropy for these pairs
  het <- all_het %>%
    semi_join(pairs, by = c("id.exposure", "id.outcome"))
  pleio <- all_pleio %>%
    semi_join(pairs, by = c("id.exposure", "id.outcome"))
  
  # Step 6: Combine all into final dataframe
  final_df <- filtered_res %>%
    select(id.exposure, id.outcome, pval) %>%
    left_join(het %>% select(id.exposure, id.outcome, qval = Q_pval), by = c("id.exposure", "id.outcome")) %>%
    left_join(pleio %>% select(id.exposure, id.outcome, pleiotropy_pval = pval), by = c("id.exposure", "id.outcome")) %>%
    left_join(rsid_df, by = c("id.exposure", "id.outcome")) %>%
    rename(rsid = SNP) %>%
    select(rsid, id.exposure, id.outcome, pval, qval, pleiotropy_pval) %>%
    filter(pval < 0.05) %>%
    distinct()
  
  # Step 7: Save only if not empty
  if (nrow(final_df) > 0) {
    file_name <- paste0("mr_sensitivity_filtered_", gsub(" ", "_", tolower(method_name)), ".csv")
    write.csv(final_df, file_name, row.names = FALSE)
    cat("✅ Saved:", file_name, "\n")
  } else {
    cat("⚠️ Skipped saving empty result for method:", method_name, "\n")
  }
}


####################################bacteria ma######################################

# Load the bac_map file once
bac_map <- read.csv("D:/final project/project work/exposure processing/mibiogen_exposure.csv")

# Create SNP-to-bac mapping
rsid_to_bac <- bac_map %>%
  group_by(SNP) %>%
  summarise(bac_all = paste(unique(bac), collapse = "; "), .groups = 'drop')

# List all method-specific MR result files
result_files <- list.files(pattern = "^mr_sensitivity_filtered_.*\\.csv$")

# Loop through and annotate each
for (file in result_files) {
  mr_results <- read.csv(file)
  
  if (nrow(mr_results) == 0) {
    cat("⏭️ Skipped empty file:", file, "\n")
    next
  }
  
  # Join with bac mapping
  mr_annotated <- mr_results %>%
    left_join(rsid_to_bac, by = c("rsid" = "SNP"))
  
  # Create new output filename
  output_file <- sub("\\.csv$", "_with_bac.csv", file)
  write.csv(mr_annotated, output_file, row.names = FALSE)
  
  cat("✅ Annotated file saved:", output_file, "\n")
}


#####################final list of bac post filteration and gene _bac count######################################



library(tidyr)
library(stringr)

# List all MR result files with bac annotations
result_files <- list.files(pattern = "^final_mr_results_with_bac_.*\\.csv$")

for (file in result_files) {
  # Extract method name for file naming
  method <- str_match(file, "final_mr_results_with_bac_(.*)\\.csv")[,2]
  
  # Step 1: Load and clean
  mr_results <- read_csv(file)
  mr_results_unique <- mr_results %>%
    distinct(rsid, id.exposure, id.outcome, .keep_all = TRUE)
  
  # Step 2: Filter significant (p < 0.05) and non-heterogeneous (qval > 0.05)
  filtered_mr <- mr_results_unique %>%
    filter(pval < 0.05, qval > 0.05)
  
  # Save filtered MR results
  write_csv(filtered_mr, paste0("filtered_mr_results_clean_final_", method, ".csv"))
  
  # Step 3: Extract and count bacterial taxa
  taxa_long <- filtered_mr %>%
    separate_rows(bac_all, sep = ";\\s*|,\\s*") %>%
    filter(!is.na(bac_all) & bac_all != "") %>%
    mutate(bac_all = str_trim(bac_all))
  
  top_taxa <- taxa_long %>%
    group_by(bac_all) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  write_csv(top_taxa, paste0("top_bacterial_taxa_clean_", method, ".csv"))
  
  # Step 4: Genus-only breakdown
  genus_only <- taxa_long %>%
    filter(str_detect(bac_all, "^genus\\.")) %>%
    mutate(genus = str_remove(bac_all, "^genus\\.")) %>%
    filter(genus != "")
  
  top_genus <- genus_only %>%
    group_by(genus) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  write_csv(top_genus, paste0("genus_level_counts_", method, ".csv"))
  
  cat("✅ Processed method:", method, "\n")
}



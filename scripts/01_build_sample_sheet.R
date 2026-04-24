# =============================================================================
# Building the sample sheet for PPMI Project 140
# Author: Melise Edwards
# Date: April 23, 2026
# Description: Merges link list with race/ethnicity metadata, constructs
#              Basename column for minfi, and filters to baseline samples
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)

# Load shared configuration
source("~/methylation/scripts/00_config.R")

# --- 1. Load metadata --------------------------------------------------------

cat("Loading link list...\n")
link_list <- read.csv(LINK_LIST, colClasses = c(SENTRIXID = "character", PATNO = "character"))
cat("Link list dimensions:", nrow(link_list), "rows x", ncol(link_list), "cols\n")
cat("Timepoints available:", paste(unique(link_list$EVENT_ID), collapse = ", "), "\n")

cat("\nLoading race/ethnicity data...\n")
race_data <- read.csv(RACE_FILE)
cat("Race/ethnicity dimensions:", nrow(race_data), "rows x", ncol(race_data), "cols\n")

# --- 2. Filter to baseline samples -------------------------------------------

cat("\nFiltering to baseline (BL) samples...\n")
baseline <- link_list %>%
  filter(EVENT_ID == PRIMARY_TIMEPOINT)
cat("Baseline samples:", nrow(baseline), "\n")

# --- 3. Merge with race/ethnicity metadata -----------------------------------

cat("\nMerging with race/ethnicity data...\n")
sample_sheet <- baseline %>%
  left_join(race_data %>% select(PATNO, HISPLAT_OL, RAASIAN_OL, 
                                  RABLACK_OL, RAHAWOPI_OL, RAWHITE_OL,
                                  RAINDALS_OL, AFICBERB_OL, ASHKJEW_OL),
            by = "PATNO")

cat("Sample sheet after merge:", nrow(sample_sheet), "rows\n")
cat("Samples missing race/ethnicity data:", 
    sum(is.na(sample_sheet$RAWHITE_OL)), "\n")

# --- 4. Construct Basename column --------------------------------------------

cat("\nConstructing Basename paths...\n")
sample_sheet <- sample_sheet %>%
  mutate(
    Basename = file.path(IDAT_DIR, SENTRIXID, 
                         paste0(SENTRIXID, "_", POSITION))
  )

# --- 5. Verify idat files exist ----------------------------------------------

cat("\nVerifying idat files exist on disk...\n")
sample_sheet <- sample_sheet %>%
  mutate(
    red_exists = file.exists(paste0(Basename, "_Red.idat")),
    grn_exists = file.exists(paste0(Basename, "_Grn.idat")),
    both_exist = red_exists & grn_exists
  )

cat("Samples with both Red and Green idat files:", 
    sum(sample_sheet$both_exist), "/", nrow(sample_sheet), "\n")

missing <- sample_sheet %>% filter(!both_exist)
if (nrow(missing) > 0) {
  cat("WARNING:", nrow(missing), "samples missing idat files:\n")
  print(missing %>% select(PATNO, SENTRIXID, POSITION))
} else {
  cat("All idat files found!\n")
}

# --- 6. Final clean sample sheet ---------------------------------------------

# Keep only samples with both idat files
sample_sheet_final <- sample_sheet %>%
  filter(both_exist) %>%
  select(-red_exists, -grn_exists, -both_exist)

cat("\nFinal sample sheet:", nrow(sample_sheet_final), "samples ready for QC\n")

# --- 7. Save sample sheet ----------------------------------------------------

output_file <- file.path(RESULTS_DIR, "sample_sheet_baseline.csv")
write.csv(sample_sheet_final, output_file, row.names = FALSE)
cat("Sample sheet saved to:", output_file, "\n")
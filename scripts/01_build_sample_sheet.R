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

# --- 2. Filter to baseline samples -------------------------------------------

cat("\nFiltering to baseline (BL) samples...\n")
baseline <- link_list %>%
  filter(EVENT_ID == "BL") # baseline only
cat("Baseline samples:", nrow(baseline), "\n")

# --- 3. Merge with participant status (diagnosis) ---------------------------

cat("\nLoading participant status...\n")
participant_status <- read.csv(PARTICIPANT_STATUS,
                               colClasses = c(PATNO = "character")) %>%
  select(PATNO, COHORT, COHORT_DEFINITION, ENROLL_AGE, ENROLL_STATUS,
         ENRLGBA, ENRLLRRK2, ENRLSNCA)

cat("Participant status dimensions:", nrow(participant_status), "rows\n")
cat("Cohort breakdown:\n")
print(table(participant_status$COHORT_DEFINITION))

# Merge baseline with participant status to create sample_sheet
sample_sheet <- baseline %>%
  left_join(participant_status, by = "PATNO")

# Filter to PD (1) and Healthy Control (2) only
sample_sheet <- sample_sheet %>%
  filter(COHORT %in% COHORTS_OF_INTEREST)

sample_sheet <- sample_sheet %>%
  mutate(DIAGNOSIS = ifelse(COHORT == 1, "PD", "HC"))

cat("Samples after filtering to PD + HC:", nrow(sample_sheet), "\n")
cat("PD:", sum(sample_sheet$COHORT == 1), "\n")
cat("HC:", sum(sample_sheet$COHORT == 2), "\n")

# --- 4. Merge with demographics (sex and race/ethnicity) -------------------

cat("\nLoading demographics (sex and race/ethnicity)...\n")
demographics <- read.csv(DEMOGRAPHICS,
                          colClasses = c(PATNO = "character")) %>%
  filter(EVENT_ID %in% c("SC", "TRANS")) %>%
  arrange(PATNO, EVENT_ID) %>% # SC before TRANS alphabetically, so SC kept by distinct()
  select(PATNO, SEX, BIRTHDT,
         HISPLAT, RAASIAN, RABLACK, RAHAWOPI,
         RAINDALS, RAWHITE, RAUNKNOWN, AFICBERB, ASHKJEW) %>%
  distinct(PATNO, .keep_all = TRUE)

sample_sheet <- sample_sheet %>%
  left_join(demographics, by = "PATNO")

sample_sheet <- sample_sheet %>%
  mutate(SEX_LABEL = case_when(
    SEX == 0 ~ "Female",
    SEX == 1 ~ "Male",
    TRUE ~ "Unknown"
  ))

cat("Sex breakdown:\n")
print(table(sample_sheet$SEX_LABEL, useNA = "always"))
cat("Samples missing race/ethnicity data:", 
    sum(is.na(sample_sheet$RAWHITE)), "\n")

# --- 5. Map SENTRIXID to actual directory paths ------------------------------

cat("\nMapping SENTRIXID to actual directory paths...\n")

# Find all directories recursively and extract their names
all_dirs <- list.dirs(IDAT_DIR, recursive = TRUE, full.names = TRUE)

# Build a lookup table: SENTRIXID -> full path
sentrix_map <- tibble(
  full_path = all_dirs,
  SENTRIXID = basename(all_dirs)) %>% 
  filter(SENTRIXID %in% sample_sheet$SENTRIXID) %>%
  distinct(SENTRIXID, .keep_all = TRUE)

cat("Found", nrow(sentrix_map), "matching SENTRIXID directories\n")

# Check for any SENTRIXIDs in our sample sheet not found on disk
missing_sentrix <- sample_sheet %>%
  filter(!SENTRIXID %in% sentrix_map$SENTRIXID) %>%
  pull(SENTRIXID) %>%
  unique()

if (length(missing_sentrix) > 0) {
  cat("WARNING:", length(missing_sentrix), 
      "SENTRIXIDs from link list not found on disk:\n")
  print(missing_sentrix)
} else {
  cat("All SENTRIXIDs found on disk!\n")
}

# Join sentrix_map to sample_sheet and construct Basename
sample_sheet <- sample_sheet %>%
  left_join(sentrix_map, by = "SENTRIXID") %>%
  mutate(
    Basename = file.path(full_path, paste0(SENTRIXID, "_", POSITION))
  ) %>%
  select(-full_path)

cat("Basename example:", sample_sheet$Basename[1], "\n")

# --- 6. Verify idat files exist ----------------------------------------------

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

# --- 7. Final clean sample sheet ---------------------------------------------

# Keep only samples with both idat files
sample_sheet_final <- sample_sheet %>%
  filter(both_exist) %>%
  select(-red_exists, -grn_exists, -both_exist)

cat("\nFinal sample sheet:", nrow(sample_sheet_final), "samples ready for QC\n")

# --- 8. Save sample sheet ----------------------------------------------------

write.csv(sample_sheet_final, SAMPLE_SHEET_BASELINE, row.names = FALSE)
cat("Sample sheet saved to:", SAMPLE_SHEET_BASELINE, "\n")
# =============================================================================
# PPMI Project 140 Methylation QC Pipeline
# Author: Melise Edwards
# Date: April 23, 2026
# Description: Unzip raw IDAT files, run QC, and generate QC report
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)
library(tidyverse)

# Define paths
RAW_DIR     <- "/mnt/output/data/ppmi/project_140/raw"
IDAT_DIR    <- "/mnt/output/data/ppmi/project_140/idat"
RESULTS_DIR <- "/mnt/output/data/ppmi/project_140/results"
META_DIR    <- "/mnt/output/data/ppmi/project_140/metadata"

# Create output directories if they don't exist
# dir.create(IDAT_DIR,    showWarnings = FALSE, recursive = TRUE)
# dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1. Unzip raw IDAT files -------------------------------------------------

# Get list of all zip files in raw directory
zip_files <- list.files(RAW_DIR, pattern = "\\.zip$", full.names = TRUE)
cat("Found", length(zip_files), "zip files to extract\n")

# Unzip each file into the idat directory
for (zip_file in zip_files) {
  cat("Unzipping:", basename(zip_file), "\n")
  unzip(zip_file, exdir = IDAT_DIR)
}

cat("All files unzipped to:", IDAT_DIR, "\n")

# Verify extraction - count idat files
idat_files <- list.files(IDAT_DIR, pattern = "\\.idat$", 
                          recursive = TRUE, full.names = TRUE)
cat("Total IDAT files found:", length(idat_files), "\n")
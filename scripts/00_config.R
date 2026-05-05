# =============================================================================
# Configuration for PPMI Project_140 methylation pipeline
# Author: Melise Edwards
# Date: April 23, 2026
# Description: Defines shared paths and parameters used across all scripts
# =============================================================================

# --- Paths -------------------------------------------------------------------

# Root data directory
PPMI_DIR    <- "/mnt/output/data/ppmi"

# Project 140 directories
P140_DIR    <- file.path(PPMI_DIR, "project_140")
RAW_DIR     <- file.path(P140_DIR, "raw")
IDAT_DIR    <- file.path(P140_DIR, "idat")
RESULTS_DIR <- file.path(P140_DIR, "results")
META_DIR    <- file.path(P140_DIR, "metadata")
SUBJ_DIR    <- file.path(META_DIR, "subject_characteristics")

# Documentation/shared metadata
DOC_DIR     <- file.path(PPMI_DIR, "documentation")

# Key metadata files
LINK_LIST          <- file.path(DOC_DIR, "ppmi_140_link_list_20210607.csv")
PARTICIPANT_STATUS <- file.path(SUBJ_DIR, "Participant_Status_23Apr2026.csv")
DEMOGRAPHICS       <- file.path(SUBJ_DIR, "Demographics_23Apr2026.csv")
AGE_AT_VISIT       <- file.path(SUBJ_DIR, "Age_at_visit_23Apr2026.csv")

# --- Cohort definitions ------------------------------------------------------
# COHORT 1 = Parkinson's Disease
# COHORT 2 = Healthy Control  
# COHORT 3 = SWEDD (Scans Without Evidence of Dopaminergic Deficit)
# COHORT 4 = Prodromal (at risk)

COHORTS_OF_INTEREST <- c(1, 2)  # PD vs Healthy Control

# --- Analysis parameters -----------------------------------------------------

# Timepoint to use for primary analysis
PRIMARY_TIMEPOINT <- "BL"  # Baseline only

# QC thresholds
DETECTION_P_THRESHOLD <- 0.01   # Max detection p-value
MIN_BEADS             <- 3      # Minimum beads per probe
FAILED_SAMPLE_CUTOFF  <- 0.1    # Max fraction of failed probes per sample (ChAMP default)

# --- Pipeline outputs --------------------------------------------------------
# These files are created by one script and read by the next

SAMPLE_SHEET_BASELINE <- file.path(RESULTS_DIR, "sample_sheet_baseline.csv")
SAMPLE_SHEET_QC       <- file.path(RESULTS_DIR, "sample_sheet_qc_passed.csv")
MSET_QC               <- file.path(RESULTS_DIR, "mSetSq_qc_passed.rds")

# Batch variable derivation note
# Batch is derived from idat folder name in Basename path
# e.g. /mnt/output/.../idat/20190604_plate1/SENTRIXID/... -> batch = "20190604_plate1"
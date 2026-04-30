# =============================================================================
# Quality Control for PPMI Project 140 Methylation Data
# Author: Melise Edwards
# Date: April 29, 2026
# Description: Reads idat files into minfi, runs QC following the minfi
#              user guide, filters failed samples and probes, normalizes
#              using preprocessFunnorm, and saves QC-passed data
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load shared configuration
source("~/methylation/scripts/00_config.R")

# Create results directory if it doesn't exist
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load sample sheet ----------------------------------------------------

cat("Loading sample sheet...\n")
targets <- read.csv(SAMPLE_SHEET_BASELINE,
                    colClasses = c(SENTRIXID = "character",
                                   PATNO     = "character",
                                   Basename  = "character")) 

cat("Samples loaded:", nrow(targets), "\n")
cat("PD:", sum(targets$COHORT == 1), "\n")
cat("HC:", sum(targets$COHORT == 2), "\n")

# --- 2. Read IDAT files ------------------------------------------------------

cat("\nReading IDAT files into minfi...\n")
rgSet <- read.metharray.exp(targets = targets,
                        verbose   = TRUE,
                        force     = TRUE)

# Verify metadata attached correctly 
cat("RGChannelSet dimensions:", dim(rgSet), "\n")
cat("Array type:", annotation(rgSet)["array"], "\n")
cat("Sample metadata columns:", ncol(pData(rgSet)), "\n")
cat("First few PATNOs:", head(pData(rgSet)$PATNO), "\n")
cat("Diagnosis breakdown:\n")
print(table(pData(rgSet)$DIAGNOSIS))

# --- 3. Initial QC -----------------------------------------------------------

cat("\nRunning initial QC...\n")

# 3a. QC plot - median methylated vs unmethylated intensity per sample
cat("Generating QC plot (median intensities)...\n")
mSet <- preprocessRaw(rgSet)
qc   <- getQC(mSet)

pdf(file.path(RESULTS_DIR, "qc_01_median_intensities.pdf"))
plotQC(qc)
dev.off()
cat("QC plot saved\n")

# 3b. Detection p-values
cat("Computing detection p-values...\n")
detP <- detectionP(rgSet)

mean_detP <- colMeans(detP)
cat("Mean detection p-value range:",
    round(min(mean_detP), 6), "to", round(max(mean_detP), 6), "\n")

# Plot mean detection p-value per sample
pdf(file.path(RESULTS_DIR, "qc_02_detection_pvalues.pdf"))
barplot(mean_detP,
        col       = ifelse(targets$COHORT == 1, "steelblue", "tomato"),
        las       = 2,
        cex.names = 0.4,
        main      = "Mean Detection P-values per Sample",
        ylab      = "Mean detection p-value")
abline(h   = DETECTION_P_THRESHOLD,
       col = "red",
       lty = 2)
legend("topright",
       legend = c("PD", "HC", "Threshold"),
       fill   = c("steelblue", "tomato", NA),
       lty    = c(NA, NA, 2),
       col    = c(NA, NA, "red"))
dev.off()
cat("Detection p-value plot saved\n")

# Identify failed samples
failed_samples <- mean_detP > DETECTION_P_THRESHOLD
cat("Samples failing detection p-value threshold:", sum(failed_samples), "\n")
if (sum(failed_samples) > 0) {
  cat("Failed samples:\n")
  print(targets$PATNO[failed_samples])
}

# 3c. Sex prediction
cat("\nPredicting sex from methylation data...\n")
mSet_mapped   <- mapToGenome(mSet)
sex_predicted <- getSex(mSet_mapped, cutoff = -2)

# Compare predicted vs reported sex (PPMI: 0=Female, 1=Male)
sex_check <- data.frame(
  PATNO        = targets$PATNO,
  reported_sex = targets$SEX,
  predicted_sex = sex_predicted$predictedSex
) %>%
  mutate(
    reported_sex_label = case_when(
      reported_sex == 0 ~ "F",
      reported_sex == 1 ~ "M",
      TRUE              ~ "Unknown"
    ),
    sex_discordant = reported_sex_label != predicted_sex &
                     reported_sex_label != "Unknown"
  )

cat("Sex discordant samples:", sum(sex_check$sex_discordant, na.rm = TRUE), "\n")
if (sum(sex_check$sex_discordant, na.rm = TRUE) > 0) {
  cat("Discordant samples:\n")
  print(sex_check %>% filter(sex_discordant))
}

# Plot sex prediction
pdf(file.path(RESULTS_DIR, "qc_03_sex_prediction.pdf"))
plot(sex_predicted$xMed, sex_predicted$yMed,
     col  = ifelse(sex_predicted$predictedSex == "F", "hotpink", "steelblue"),
     pch  = 16,
     xlab = "X chromosome median intensity",
     ylab = "Y chromosome median intensity",
     main = "Sex Prediction")
# Add text labels for discordant samples
discordant_idx <- !is.na(sex_check$sex_discordant) & sex_check$sex_discordant
text(sex_predicted$xMed[discordant_idx], 
     sex_predicted$yMed[discordant_idx],
     labels = sex_check$PATNO[discordant_idx],
     pos    = 3,
     cex    = 0.7,
     col    = "red")
legend("topright", 
       legend = c("Female", "Male", "Discordant"),
       col    = c("hotpink", "steelblue", "red"), 
       pch    = c(16, 16, 16))
dev.off()
cat("Sex prediction plot saved\n")

# --- 4. Remove failed samples ------------------------------------------------

cat("\nRemoving failed samples...\n")
samples_to_remove <- failed_samples | sex_check$sex_discordant
cat("Total samples removed:", sum(samples_to_remove), "\n")
cat("Samples remaining:", sum(!samples_to_remove), "\n")

rgSet_clean   <- rgSet[, !samples_to_remove]
detP_clean    <- detP[,  !samples_to_remove]
targets_clean <- targets[!samples_to_remove, ]

# --- 5. Normalization --------------------------------------------------------

cat("\nNormalizing with preprocessFunnorm...\n")

# Density plot BEFORE normalization (fix 2: use getBeta(preprocessRaw()))
pdf(file.path(RESULTS_DIR, "qc_04_density_before_normalization.pdf"))
densityPlot(getBeta(preprocessRaw(rgSet_clean)),
            sampGroups = targets_clean$DIAGNOSIS,
            main       = "Beta Values - Before Normalization",
            legend     = TRUE)
dev.off()

# Apply functional normalization
mSetSq <- preprocessFunnorm(rgSet_clean)
cat("Normalized object dimensions:", dim(mSetSq), "\n")

# Density plot AFTER normalization
pdf(file.path(RESULTS_DIR, "qc_05_density_after_normalization.pdf"))
densityPlot(getBeta(mSetSq),
            sampGroups = targets_clean$DIAGNOSIS,
            main       = "Beta Values - After Normalization",
            legend     = TRUE)
dev.off()
cat("Density plots saved\n")

# --- 6. Probe Filtering -------------------------------------------------------

cat("\nFiltering probes...\n")
n_probes_start <- nrow(mSetSq)

# Ensure detP rows match mSetSq probe order after normalization
detP_clean <- detP_clean[match(featureNames(mSetSq), rownames(detP_clean)), ]

# 6a. Remove failed probes
cat("Removing failed probes...\n")
failed_probes <- rowSums(detP_clean > DETECTION_P_THRESHOLD) >
                 (FAILED_SAMPLE_CUTOFF * ncol(detP_clean))
mSetSq <- mSetSq[!failed_probes, ]
n_after_failed <- nrow(mSetSq)  # fix 3: track counts progressively
cat("Probes removed (failed detection):", sum(failed_probes), "\n")

# 6b. Remove SNP-overlapping probes
cat("Removing SNP-overlapping probes...\n")
mSetSq <- dropLociWithSnps(mSetSq)
n_after_snp <- nrow(mSetSq)  # fix 3
cat("Probes removed (SNP-overlapping):", n_after_failed - n_after_snp, "\n")

# 6c. Remove sex chromosome probes
cat("Removing sex chromosome probes...\n")
ann_epic   <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
sex_probes <- ann_epic$Name[ann_epic$chr %in% c("chrX", "chrY")]
mSetSq     <- mSetSq[!rownames(mSetSq) %in% sex_probes, ]
n_after_sex <- nrow(mSetSq)  # fix 3
cat("Probes removed (sex chromosomes):", n_after_snp - n_after_sex, "\n")

cat("Total probes removed:", n_probes_start - n_after_sex, "\n")
cat("Final probe count:", n_after_sex, "\n")

# --- 7. Calculate M/Beta values and Save QC-passed data -------------------------

cat("\nSaving QC-passed data...\n")

cat("\nCalculating M and Beta values...\n")
mVals <- getM(mSetSq)
bVals <- getBeta(mSetSq)

cat("M values dimensions:", dim(mVals), "\n")
cat("Beta values dimensions:", dim(bVals), "\n")

# Save as RDS (more efficient than CSV for large matrices)
saveRDS(mVals, file.path(RESULTS_DIR, "mVals.rds"))
saveRDS(bVals, file.path(RESULTS_DIR, "bVals.rds"))
cat("M and Beta values saved\n")

# Save normalized filtered GenomicRatioSet
saveRDS(mSetSq, MSET_QC)

# Save cleaned sample sheet
write.csv(targets_clean, SAMPLE_SHEET_QC, row.names = FALSE)

# Save QC summary (fix 3: use progressive counts)
qc_summary <- data.frame(
  step = c("Input samples", "Failed detection p-value",
            "Sex discordant", "Final samples",
            "Input probes", "Failed probes removed",
            "SNP probes removed", "Sex chromosome probes removed",
            "Final probes"),
  n    = c(nrow(targets), sum(failed_samples),
            sum(sex_check$sex_discordant, na.rm = TRUE),
            nrow(targets_clean),
            n_probes_start,
            n_probes_start - n_after_failed,
            n_after_failed - n_after_snp,
            n_after_snp - n_after_sex,
            n_after_sex)
)

write.csv(qc_summary,
          file.path(RESULTS_DIR, "qc_summary.csv"),
          row.names = FALSE)

cat("QC summary saved\n")
cat("\nQC pipeline complete!\n")
cat("Results saved to:", RESULTS_DIR, "\n")
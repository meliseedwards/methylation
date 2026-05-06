# =============================================================================
# Batch Correction and Cell Type Deconvolution for PPMI Project 140
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: May 6, 2026
# Description: Loads QC-passed methylation data, performs PCA exploration,
#              confounder analysis, cell type deconvolution, and batch
#              correction via ComBat. Saves corrected M values for use
#              in differential methylation analysis.
# =============================================================================

# --- 0. Setup ----------------------------------------------------------------

library(minfi)
library(tidyverse)
library(sva)                      # ComBat
library(FlowSorted.Blood.EPIC)    # cell type deconvolution
library(ExperimentHub)            # required by FlowSorted.Blood.EPIC
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load shared configuration
source("~/methylation/scripts/00_config.R")

# --- 1. Load QC-passed data --------------------------------------------------

cat("Loading QC-passed data from Script 02...\n")

# Load the GenomicRatioSet (normalized, filtered)
mSetSq <- readRDS(MSET_QC)
cat("GenomicRatioSet dimensions:", dim(mSetSq), "\n")

# Load M values
mVals <- readRDS(file.path(RESULTS_DIR, "mVals.rds"))
cat("M values dimensions:", dim(mVals), "\n")

# Load QC-passed sample sheet
targets_clean <- read.csv(SAMPLE_SHEET_QC,
                          colClasses = c(SENTRIXID = "character",
                                         PATNO     = "character",
                                         Basename  = "character"))
cat("Sample sheet rows:", nrow(targets_clean), "\n")

# Sample order consistency check 
# colnames of mVals should match the sample order in targets_clean
cat("\n--- Sample order consistency check ---\n")
cat("mVals columns match Basename order:", 
    all(colnames(mVals) == basename(targets_clean$Basename)), "\n")

# If not aligned, reorder targets_clean to match mVals columns
if (!all(colnames(mVals) == basename(targets_clean$Basename))) {
  cat("Reordering sample sheet to match M values column order...\n")
  targets_clean <- targets_clean[match(colnames(mVals), 
                                       basename(targets_clean$Basename)), ]
}

# Derive Batch variable from Basename path
# Path format: /mnt/.../idat/<BATCH_FOLDER>/<SENTRIXID>/<SENTRIXID_POSITION>
targets_clean <- targets_clean %>%
  mutate(Batch = basename(dirname(dirname(Basename))))

cat("\nBatch breakdown:\n")
print(table(targets_clean$Batch))
cat("Number of unique batches:", length(unique(targets_clean$Batch)), "\n")

# Flag small batches for downstream awareness
small_batches <- table(targets_clean$Batch) %>%
  as.data.frame() %>%
  filter(Freq < 5)
if (nrow(small_batches) > 0) {
  cat("\nWARNING: Small batches (<5 samples):\n")
  print(small_batches)
}

# --- 2. PCA exploration BEFORE batch correction ------------------------------

cat("\n--- PCA exploration (before correction) ---\n")

# Use top 1000 most variable probes for PCA (standard EWAS practice)
n_top_probes <- 1000
probe_vars <- rowVars(mVals)
top_var_probes <- order(probe_vars, decreasing = TRUE)[1:n_top_probes]

# PCA on transposed matrix (samples as rows for prcomp)
pca <- prcomp(t(mVals[top_var_probes, ]), scale. = TRUE, center = TRUE)

# Calculate variance explained
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
cat("Variance explained by PC1:", round(var_explained[1], 2), "%\n")
cat("Variance explained by PC2:", round(var_explained[2], 2), "%\n")

# Build a dataframe for plotting
pca_df <- data.frame(
  PC1       = pca$x[, 1],
  PC2       = pca$x[, 2],
  Batch     = targets_clean$Batch,
  DIAGNOSIS = targets_clean$DIAGNOSIS,
  SEX_LABEL = targets_clean$SEX_LABEL,
  AGE       = targets_clean$ENROLL_AGE
)

# 2a. PCA colored by Batch
pdf(file.path(RESULTS_DIR, "qc_06_pca_before_by_batch.pdf"),
    width = 10, height = 8)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA Before Batch Correction - Colored by Batch",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 7))
dev.off()

# 2b. PCA colored by Diagnosis
pdf(file.path(RESULTS_DIR, "qc_07_pca_before_by_diagnosis.pdf"),
    width = 10, height = 8)
ggplot(pca_df, aes(x = PC1, y = PC2, color = DIAGNOSIS)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("HC" = "tomato", "PD" = "steelblue")) +
  labs(
    title = "PCA Before Batch Correction - Colored by Diagnosis",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw()
dev.off()

# 2c. PCA colored by Sex
pdf(file.path(RESULTS_DIR, "qc_08_pca_before_by_sex.pdf"),
    width = 10, height = 8)
ggplot(pca_df, aes(x = PC1, y = PC2, color = SEX_LABEL)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Female" = "hotpink", "Male" = "steelblue")) +
  labs(
    title = "PCA Before Batch Correction - Colored by Sex",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_bw()
dev.off()

# 2d. PCA colored by Age (continuous)
pdf(file.path(RESULTS_DIR, "qc_09_pca_before_by_age.pdf"),
    width = 10, height = 8)
ggplot(pca_df, aes(x = PC1, y = PC2, color = AGE)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_viridis_c() +
  labs(
    title = "PCA Before Batch Correction - Colored by Age",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)"),
    color = "Age at\nEnrollment"
  ) +
  theme_bw()
dev.off()

cat("PCA plots saved\n")

# --- 3. Statistical confounder analysis --------------------------------------

cat("\n--- Statistical confounder analysis ---\n")

# 3a. Age by diagnosis (t-test)
age_test <- t.test(ENROLL_AGE ~ DIAGNOSIS, data = targets_clean)
cat("\nAge by Diagnosis (t-test):\n")
cat("  PD mean age:", 
    round(mean(targets_clean$ENROLL_AGE[targets_clean$DIAGNOSIS == "PD"], 
               na.rm = TRUE), 2), "\n")
cat("  HC mean age:", 
    round(mean(targets_clean$ENROLL_AGE[targets_clean$DIAGNOSIS == "HC"], 
               na.rm = TRUE), 2), "\n")
cat("  t-statistic:", round(age_test$statistic, 3), "\n")
cat("  p-value:", format(age_test$p.value, scientific = TRUE, digits = 3), "\n")

# 3b. Sex by diagnosis (chi-square)
sex_table <- table(targets_clean$DIAGNOSIS, targets_clean$SEX_LABEL)
sex_test  <- chisq.test(sex_table)
cat("\nSex by Diagnosis (chi-square):\n")
print(sex_table)
cat("  chi-square:", round(sex_test$statistic, 3), "\n")
cat("  p-value:", format(sex_test$p.value, scientific = TRUE, digits = 3), "\n")

# 3c. Batch by diagnosis (chi-square)
batch_table <- table(targets_clean$DIAGNOSIS, targets_clean$Batch)
batch_test  <- chisq.test(batch_table)
cat("\nBatch by Diagnosis (chi-square):\n")
cat("  chi-square:", round(batch_test$statistic, 3), "\n")
cat("  p-value:", format(batch_test$p.value, scientific = TRUE, digits = 3), "\n")

# Save confounder analysis summary
confounder_summary <- data.frame(
  variable    = c("Age (continuous)", "Sex (categorical)", "Batch (categorical)"),
  test        = c("t-test", "chi-square", "chi-square"),
  statistic   = c(age_test$statistic, sex_test$statistic, batch_test$statistic),
  p_value     = c(age_test$p.value, sex_test$p.value, batch_test$p.value),
  significant = c(age_test$p.value, sex_test$p.value, batch_test$p.value) < 0.05
)
write.csv(confounder_summary,
          file.path(RESULTS_DIR, "confounder_summary.csv"),
          row.names = FALSE)
cat("\nConfounder analysis saved to confounder_summary.csv\n")

# --- 4. Cell type deconvolution ----------------------------------------------

cat("\n--- Cell type deconvolution ---\n")
cat("Using FlowSorted.Blood.EPIC (Salas et al. 2018 IDOL probes)\n")

# We need to load the original RGChannelSet for cell type estimation
# estimateCellCounts2 requires the raw RGSet, not the normalized GenomicRatioSet
# We will need to reload the rgSet from idat files
# This is a known limitation - cell type deconvolution requires raw intensities

# NOTE: For this step we need to read idat files again
# Alternative: save rgSet from Script 02 (would require modification)
# For now, we'll re-read idat files

cat("Re-reading idat files for cell type estimation...\n")
rgSet <- read.metharray.exp(targets  = targets_clean,
                            verbose  = TRUE,
                            force    = TRUE,
                            extended = TRUE)

cat("Estimating cell counts using IDOL probes...\n")
cell_counts <- estimateCellCounts2(
  rgSet,
  compositeCellType  = "Blood",
  processMethod      = "preprocessNoob",
  probeSelect        = "IDOL",
  cellTypes          = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
  referencePlatform  = "IlluminaHumanMethylationEPIC",
  returnAll          = FALSE
)

# estimateCellCounts2 returns a list; we want the counts matrix
cell_proportions <- as.data.frame(cell_counts$counts)
cell_proportions$Basename_short <- rownames(cell_proportions)

cat("Cell type proportions estimated for", 
    nrow(cell_proportions), "samples\n")
cat("Cell type means across all samples:\n")
print(round(colMeans(cell_proportions[, 1:6]), 4))

# Add cell type proportions to targets_clean
# Sample order consistency check
cat("\nVerifying cell counts sample order matches targets_clean...\n")
cell_proportions <- cell_proportions[
  match(basename(targets_clean$Basename), cell_proportions$Basename_short), ]

stopifnot(all(cell_proportions$Basename_short == 
              basename(targets_clean$Basename)))
cat("Sample order verified\n")

# Bind cell counts to sample sheet
targets_clean <- bind_cols(targets_clean,
                            cell_proportions[, c("CD8T", "CD4T", "NK", 
                                                 "Bcell", "Mono", "Neu")])

# Plot cell type proportions by diagnosis
cell_long <- targets_clean %>%
  select(PATNO, DIAGNOSIS, CD8T, CD4T, NK, Bcell, Mono, Neu) %>%
  pivot_longer(cols = c(CD8T, CD4T, NK, Bcell, Mono, Neu),
               names_to  = "cell_type",
               values_to = "proportion")

pdf(file.path(RESULTS_DIR, "qc_10_cell_proportions.pdf"),
    width = 10, height = 7)
ggplot(cell_long, aes(x = cell_type, y = proportion, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("HC" = "tomato", "PD" = "steelblue")) +
  labs(
    title = "Cell Type Proportions by Diagnosis",
    x     = "Cell Type",
    y     = "Estimated Proportion"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
cat("Cell proportions plot saved\n")

# --- 5. ComBat batch correction ----------------------------------------------

cat("\n--- ComBat batch correction ---\n")
cat("Note: This is a first-pass batch correction.\n")
cat("With Project 120 data added (~500 more samples),\n")
cat("we will revisit this with a more balanced design.\n")

# Sample order consistency check before ComBat
cat("\nFinal sample order verification before ComBat:\n")
cat("colnames(mVals) match Basename in targets_clean:", 
    all(colnames(mVals) == basename(targets_clean$Basename)), "\n")

# Build model matrix protecting biological variable (DIAGNOSIS)
# This tells ComBat to preserve DIAGNOSIS-related variation
mod <- model.matrix(~ DIAGNOSIS, data = targets_clean)

# Apply ComBat to M values
# Following Gonzalez-Latapi et al. 2023 PPMI methodology
cat("Running ComBat...\n")
combat_mVals <- ComBat(
  dat         = as.matrix(mVals),
  batch       = targets_clean$Batch,
  mod         = mod,
  par.prior   = TRUE,
  prior.plots = FALSE
)

cat("ComBat complete\n")
cat("Corrected M values dimensions:", dim(combat_mVals), "\n")

# --- 6. PCA exploration AFTER batch correction -------------------------------

cat("\n--- PCA exploration (after correction) ---\n")

# Recompute PCA on corrected M values
top_var_probes_post <- order(rowVars(combat_mVals), 
                              decreasing = TRUE)[1:n_top_probes]
pca_post <- prcomp(t(combat_mVals[top_var_probes_post, ]), 
                    scale. = TRUE, center = TRUE)
var_explained_post <- (pca_post$sdev^2 / sum(pca_post$sdev^2)) * 100

cat("Variance explained by PC1 (after):", 
    round(var_explained_post[1], 2), "%\n")
cat("Variance explained by PC2 (after):", 
    round(var_explained_post[2], 2), "%\n")

pca_df_post <- data.frame(
  PC1       = pca_post$x[, 1],
  PC2       = pca_post$x[, 2],
  Batch     = targets_clean$Batch,
  DIAGNOSIS = targets_clean$DIAGNOSIS,
  SEX_LABEL = targets_clean$SEX_LABEL,
  AGE       = targets_clean$ENROLL_AGE
)

# 6a. PCA after - by Batch (should show LESS clustering by batch)
pdf(file.path(RESULTS_DIR, "qc_11_pca_after_by_batch.pdf"),
    width = 10, height = 8)
ggplot(pca_df_post, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA After Batch Correction - Colored by Batch",
    subtitle = "Successful correction = reduced clustering by batch",
    x     = paste0("PC1 (", round(var_explained_post[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained_post[2], 2), "%)")
  ) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 7))
dev.off()

# 6b. PCA after - by Diagnosis (should preserve biology)
pdf(file.path(RESULTS_DIR, "qc_12_pca_after_by_diagnosis.pdf"),
    width = 10, height = 8)
ggplot(pca_df_post, aes(x = PC1, y = PC2, color = DIAGNOSIS)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("HC" = "tomato", "PD" = "steelblue")) +
  labs(
    title = "PCA After Batch Correction - Colored by Diagnosis",
    subtitle = "Successful correction = biological signal preserved",
    x     = paste0("PC1 (", round(var_explained_post[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained_post[2], 2), "%)")
  ) +
  theme_bw()
dev.off()

cat("Post-correction PCA plots saved\n")

# --- 7. Save outputs ---------------------------------------------------------

cat("\n--- Saving outputs ---\n")

# Save batch-corrected M values
saveRDS(combat_mVals, file.path(RESULTS_DIR, "combat_mVals.rds"))
cat("Batch-corrected M values saved\n")

# Save updated sample sheet with cell counts and batch
write.csv(targets_clean,
          file.path(RESULTS_DIR, "sample_sheet_final.csv"),
          row.names = FALSE)
cat("Final sample sheet saved (with cell counts and batch)\n")

# Save analysis summary
analysis_summary <- data.frame(
  metric = c("Final samples",
             "Final probes (autosomal)",
             "Number of batches",
             "Smallest batch size",
             "Age difference (PD vs HC) p-value",
             "Sex distribution p-value",
             "Batch distribution p-value",
             "PC1 variance explained (before)",
             "PC1 variance explained (after)"),
  value  = c(nrow(targets_clean),
             nrow(combat_mVals),
             length(unique(targets_clean$Batch)),
             min(table(targets_clean$Batch)),
             format(age_test$p.value, scientific = TRUE, digits = 3),
             format(sex_test$p.value, scientific = TRUE, digits = 3),
             format(batch_test$p.value, scientific = TRUE, digits = 3),
             round(var_explained[1], 2),
             round(var_explained_post[1], 2))
)
write.csv(analysis_summary,
          file.path(RESULTS_DIR, "batch_correction_summary.csv"),
          row.names = FALSE)

cat("\nScript 03 complete!\n")
cat("Outputs saved to:", RESULTS_DIR, "\n")
cat("\nKey files for Script 04:\n")
cat("  - combat_mVals.rds: batch-corrected M values\n")
cat("  - sample_sheet_final.csv: sample sheet with cell counts and batch\n")
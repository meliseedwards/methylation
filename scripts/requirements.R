# =============================================================================
# requirements.R - R/Bioconductor packages not managed by conda
# Author: GP2 Subtypes and Mechanisms - M.E.
# Date: April 2026
# =============================================================================

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "minfi",
    "minfiData",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "IlluminaHumanMethylationEPICmanifest",
    "DMRcate",
    "limma",
    "bumphunter"
), ask = FALSE)

# Install preprocessCore without threading (required for Linux/GCP environments)
BiocManager::install("preprocessCore",
                     configure.args = "--disable-threading",
                     force = TRUE)

# CRAN packages
install.packages("remotes")

# GitHub package for cross-reactive probe removal - depends on minfiData
remotes::install_github("markgene/maxprobes")
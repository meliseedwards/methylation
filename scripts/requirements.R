# =============================================================================
# requirements.R - R/Bioconductor packages not managed by conda
# Author: Melise Edwards
# Date: April 2026
# Usage: Rscript requirements.R
# =============================================================================

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "minfi",
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

# Install maxprobes from GitHub for cross-reactive probe removal
if (!require("devtools")) install.packages("devtools")
devtools::install_github("markgene/maxprobes")
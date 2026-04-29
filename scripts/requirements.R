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

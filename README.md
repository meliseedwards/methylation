# methylation
methylation data and scripts 

## Pipeline Order

Run scripts in this order:

1. `scripts/01_build_sample_sheet.R` — builds sample sheet from metadata
2. `scripts/02_qc.R` — QC, normalization, probe filtering
3. `scripts/03_normalize.R` — batch correction, cell type deconvolution
4. `scripts/04_differential_methylation.R` — DMPs and DMRs

Each script reads input files produced by the previous script.
All file paths are defined in `scripts/00_config.R`.

## Setup
```bash
conda env create -f environment.yml
conda activate methylation
Rscript scripts/requirements.R
```

## Troubleshooting

### preprocessFunnorm pthread_create error on Linux/GCP VMs

If you encounter this error during normalization:

Reinstall `preprocessCore` with threading disabled:
```r
BiocManager::install("preprocessCore", 
                     configure.args = "--disable-threading",
                     force = TRUE)
```

This is a known issue with `preprocessCore` in conda R environments on Linux. The reinstallation only needs to be done once per environment.
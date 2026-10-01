# =============================================================================
# run_all.R
# Reproduces every statistical analysis and data figure of:
#   Ribes Martinez E, Viana AES, et al. "A biomimetic 3D bovine endometrial gland–stroma
#   model responds to ovarian hormones and embryonic signals" (revised manuscript).
#
# Usage: open this folder in RStudio (or setwd() to it) and run
#   source("run_all.R")
# Outputs are written to output/ (tables, figures, key_results.txt, sessionInfo.txt).
# =============================================================================
rm(list = ls())
source("R/00_setup.R")
source("R/01_rnaseq.R")
source("R/02_enrichment.R")
source("R/03_nanoindentation.R")
source("R/04_gland_area.R")
source("R/05_embryo.R")
source("R/06_assemble_figures.R")
source("R/07_proteomics.R")
writeLines(capture.output(sessionInfo()), file.path(dir_out, "sessionInfo.txt"))
message("Done. See output/key_results.txt and output/figures/.")

# =============================================================================
# 00_setup.R — packages, paths, plotting defaults and helpers
# =============================================================================

cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "tibble", "lme4", "lmerTest",
               "emmeans", "ggrepel", "openxlsx", "patchwork", "scales", "png")
bioc_pkgs <- c("DESeq2", "clusterProfiler", "org.Bt.eg.db", "org.Hs.eg.db", "ReactomePA", "AnnotationDbi")

install_if_missing <- function() {
  miss_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_cran)) install.packages(miss_cran)
  miss_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(miss_bioc, update = FALSE, ask = FALSE)
  }
}
install_if_missing()

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(tibble)
  library(lme4); library(lmerTest); library(emmeans); library(ggrepel)
  library(openxlsx); library(patchwork); library(scales)
  library(DESeq2)
})

set.seed(20260926)
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 5000)

# ---- paths ------------------------------------------------------------------
if (!file.exists("run_all.R")) stop("Please set the working directory to the repository root (the folder containing run_all.R).")
dir_data  <- "data"
dir_out   <- "output"
dir_tab   <- file.path(dir_out, "tables")
dir_fig   <- file.path(dir_out, "figures")
dir_panel <- file.path(dir_fig, "panels")
for (d in c(dir_out, dir_tab, dir_fig, dir_panel)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- results log (every number quoted in the manuscript is written here) -----
results_file <- file.path(dir_out, "key_results.txt")
if (!exists("results_log_started")) {
  writeLines(c("Key results for the manuscript", paste("Generated:", Sys.time()), ""), results_file)
  results_log_started <- TRUE
}
log_result <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(txt, "\n", file = results_file, append = TRUE)
}
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), formatC(p, format = "f", digits = 3)))

# ---- plotting defaults --------------------------------------------------------
theme_paper <- function(base = 8) {
  theme_classic(base_size = base) +
    theme(axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(face = "bold", size = base + 1),
          plot.tag = element_text(face = "bold", size = base + 3))
}
theme_set(theme_paper())

col_group <- c("No cells" = "#9E9E9E", "Cells" = "#F8766D", "Cells + hormones" = "#00BFC4",
               "Native tissue" = "#7B3294")
col_hormone <- c("No hormone" = "#F8766D", "Hormone" = "#00BFC4")
col_embryo  <- c("No embryo" = "#F8766D", "Embryo" = "#3B7FB6")
col_hatch   <- c("2D well" = "#9E9E9E", "Cell-free hydrogel" = "#A6CEE3", "Co-culture" = "#F8766D")

save_panel <- function(p, name, w = 90, h = 70) {
  ggsave(file.path(dir_panel, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(dir_panel, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600)
  invisible(p)
}

p_label <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))

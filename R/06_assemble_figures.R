# =============================================================================
# 06_assemble_figures.R — multi-panel Figures 3, 4, 5 and Supplementary Figure 4
# Panel letters follow the published figure legends.
# =============================================================================

img_panel <- function(path) {
  if (!file.exists(path)) return(ggplot() + annotate("text", x = 0, y = 0, label = paste("missing:", basename(path))) + theme_void())
  wrap_elements(full = grid::rasterGrob(png::readPNG(path), interpolate = TRUE))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(dir_fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600)
  ggsave(file.path(dir_fig, paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
}

# Figure 3: a timeline (schematic), b Young's modulus, c–e viscoelasticity
fig3 <- (img_panel(file.path(dir_data, "figure_assets", "Fig3a_timeline.png")) | p_ym) /
        (p_ve[["Cells"]] | p_ve[["Cells + hormones"]] | p_ve[["No cells"]]) +
        plot_layout(heights = c(1, 0.9), guides = "collect") + plot_annotation(tag_levels = "a")
save_fig(fig3, "Figure3", 180, 150)

# Figure 4: a gland area, b PCA, c MA, d volcano, e heatmap
left4 <- p_area_h / (plots_hormone$pca | plots_hormone$ma) / plots_hormone$volcano + plot_layout(heights = c(1, 0.9, 1))
fig4 <- (left4 | plots_hormone$heatmap) + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig4, "Figure4", 180, 210)

# Figure 5: a gland area, b hatching, c PCA, d MA, e volcano, f heatmap
left5 <- (p_area_e | p_hat) / (plots_embryo$pca | plots_embryo$ma) / plots_embryo$volcano + plot_layout(heights = c(1, 0.9, 1))
fig5 <- (left5 | plots_embryo$heatmap) + plot_layout(widths = c(2, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig5, "Figure5", 180, 210)

# Supplementary Figure 4: all nine libraries (incl. excluded)
sfig4 <- (p_pca_all | p_cor_all) + plot_layout(widths = c(1, 1.1)) + plot_annotation(tag_levels = "a")
save_fig(sfig4, "Supplementary_Figure4_RNAseq_all_libraries", 180, 85)

# Supplementary Figures 2 and 3 (Reactome) are written by 02_enrichment.R to output/figures/panels
for (nm in c("SuppFig2_Reactome_hormone", "SuppFig3_Reactome_embryo")) {
  f <- file.path(dir_panel, paste0(nm, ".pdf"))
  if (file.exists(f)) file.copy(f, file.path(dir_fig, paste0(nm, ".pdf")), overwrite = TRUE)
  f <- file.path(dir_panel, paste0(nm, ".png"))
  if (file.exists(f)) file.copy(f, file.path(dir_fig, paste0(nm, ".png")), overwrite = TRUE)
}

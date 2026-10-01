# =============================================================================
# 07_proteomics.R — secondary re-analysis of the EndoECM proteome (Figure 1, Supplementary Figure 1)
#
# Data: DIA-MS MaxLFQ intensities of decellularised endometrium hydrogels (ref. 28; PRIDE PXD053248),
# 10 samples: 5 follicular-phase (EF*) and 5 luteal-phase (EL*) EndoECM preparations.
# Proteins detected (non-missing) in at least 6 of the 10 samples are retained; remaining
# missing values are set to 0; proteins are classified with the Bos taurus matrisome list
# (Matrisome Project). Relative abundance = share of the summed mean intensity within the
# matrisome (Figure 1a) or within each matrisome category (Figure 1b, Supplementary Figure 1).
# =============================================================================

prot <- read_tsv(file.path(dir_data, "proteomics", "endoECM_proteomics_maxLFQ.tsv"), show_col_types = FALSE)
matri <- read_csv(file.path(dir_data, "proteomics", "Bt_Matrisome.csv"), show_col_types = FALSE) |> dplyr::distinct(Genes, .keep_all = TRUE)
samp_cols <- grep("^MaxLFQIntensity\\.E", names(prot), value = TRUE)
min_detect <- 6

prot[samp_cols] <- lapply(prot[samp_cols], function(v) suppressWarnings(as.numeric(v)))
det <- rowSums(!is.na(as.matrix(prot[samp_cols])))
kept <- prot[det >= min_detect, ]
kept[samp_cols] <- lapply(kept[samp_cols], function(v) replace(v, is.na(v), 0))
kept$mean_intensity <- rowMeans(as.matrix(kept[samp_cols]))
mat <- dplyr::inner_join(kept, matri, by = "Genes")

log_result("")
log_result("== EndoECM proteomics (", length(samp_cols), " samples: ", sum(grepl("\\.EF", samp_cols)), " follicular, ",
           sum(grepl("\\.EL", samp_cols)), " luteal; detected in >= ", min_detect, " samples)")
log_result("  proteins retained: ", nrow(kept), "; matrisome proteins: ", nrow(mat), " (",
           sum(mat$Matrisome_Division == "Core matrisome"), " core matrisome, ",
           sum(mat$Matrisome_Division == "Matrisome-associated"), " matrisome-associated)")
cats <- mat |> dplyr::count(Matrisome_Category, name = "proteins")
comp <- mat |> dplyr::group_by(Matrisome_Category) |> dplyr::summarise(intensity = sum(mean_intensity), .groups = "drop") |>
  dplyr::mutate(percent = 100 * intensity / sum(intensity)) |> dplyr::left_join(cats, by = "Matrisome_Category") |> dplyr::arrange(desc(percent))
write_csv(comp, file.path(dir_tab, "proteomics_matrisome_composition.csv"))
for (i in seq_len(nrow(comp))) log_result(sprintf("  %-18s %3d proteins, %.1f%% of matrisome intensity",
                                                  comp$Matrisome_Category[i], comp$proteins[i], comp$percent[i]))
within <- mat |> dplyr::group_by(Matrisome_Category) |> dplyr::mutate(percent_within = 100 * mean_intensity / sum(mean_intensity)) |>
  ungroup() |> dplyr::arrange(Matrisome_Category, desc(percent_within))
write_csv(dplyr::select(within, Matrisome_Division, Matrisome_Category, Genes, mean_intensity, percent_within),
          file.path(dir_tab, "proteomics_relative_abundance_within_category.csv"))
top_txt <- within |> dplyr::group_by(Matrisome_Category) |> dplyr::slice_head(n = 3) |>
  dplyr::summarise(txt = paste0(Genes, " ", sprintf("%.1f%%", percent_within), collapse = ", "), .groups = "drop")
for (i in seq_len(nrow(top_txt))) log_result("  top in ", top_txt$Matrisome_Category[i], ": ", top_txt$txt[i])
col1 <- sum(within$percent_within[within$Genes %in% c("COL1A1", "COL1A2")])
col6 <- sum(within$percent_within[grepl("^COL6", within$Genes)])
log_result(sprintf("  collagen I (COL1A1 + COL1A2) = %.1f%% of collagens; collagen VI = %.1f%%; laminin isoforms detected: %d",
                   col1, col6, sum(grepl("^LAM", mat$Genes))))

# Figure 1a — composition
cat_cols <- c("Collagens" = "#66C2A5", "ECM Glycoproteins" = "#FC8D62", "Proteoglycans" = "#8DA0CB",
              "ECM Regulators" = "#E78AC3", "ECM-affiliated" = "#A6D854", "Secreted Factors" = "#FFD92F")
comp$Matrisome_Category <- factor(comp$Matrisome_Category, levels = names(cat_cols))
p_comp <- ggplot(comp, aes(x = 1, y = percent, fill = Matrisome_Category)) +
  geom_col(width = 1, colour = "white") + coord_polar(theta = "y") +
  geom_text(data = dplyr::filter(comp, percent >= 5), aes(label = sprintf("%.1f%%", percent)),
            position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = cat_cols, labels = setNames(sprintf("%s (%.1f%%)", comp$Matrisome_Category, comp$percent), comp$Matrisome_Category)) +
  labs(fill = "Matrisome category", title = "Matrisome composition (%)") +
  theme_void(base_size = 8) + theme(legend.key.size = unit(3, "mm"), plot.tag = element_text(face = "bold", size = 11))
save_panel(p_comp, "Fig1a_matrisome_composition", 90, 65)

bar_panel <- function(cat, title) {
  d <- dplyr::filter(within, Matrisome_Category == cat) |>
    dplyr::group_by(Genes) |> dplyr::summarise(percent_within = sum(percent_within), .groups = "drop") |>   # protein groups sharing a gene symbol are summed
    dplyr::arrange(dplyr::desc(percent_within)) |> dplyr::mutate(Genes = factor(Genes, levels = rev(Genes)))
  ggplot(d, aes(percent_within, Genes)) + geom_col(fill = cat_cols[[cat]]) +
    geom_text(aes(label = sprintf("%.2f%%", percent_within)), hjust = -0.1, size = 1.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Relative abundance (%)", y = NULL, title = title) +
    theme(axis.text.y = element_text(size = if (nrow(d) > 40) 3.5 else 5.5))
}
p_col <- bar_panel("Collagens", "Collagens")
save_panel(p_col, "Fig1b_collagens", 90, 80)
fig1 <- (p_comp | p_col) + plot_annotation(tag_levels = "a")
save_fig <- function(p, name, w, h) {
  ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(dir_fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600)
  ggsave(file.path(dir_fig, paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw")
}
save_fig(fig1, "Figure1", 180, 80)
sfig1 <- ((bar_panel("ECM Glycoproteins", "ECM glycoproteins") / bar_panel("Proteoglycans", "Proteoglycans") /
           bar_panel("ECM Regulators", "ECM regulators") + plot_layout(heights = c(1.6, 0.6, 1.2))) |
         (bar_panel("ECM-affiliated", "ECM-affiliated proteins") / bar_panel("Secreted Factors", "Secreted factors") / plot_spacer() +
           plot_layout(heights = c(0.9, 0.7, 1.8)))) + plot_annotation(tag_levels = "a")
save_fig(sfig1, "Supplementary_Figure1_matrisome_categories", 180, 240)

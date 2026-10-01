# =============================================================================
# 01_rnaseq.R — differential expression (DESeq2), PCA, MA, volcano, heatmaps
#
# Settings are kept identical to the analysis reported in the revised manuscript:
#   * genes whose identifier starts with LOC or NON are removed
#   * genes with < 10 total counts across the libraries of a contrast are removed
#   * design ~ condition, Wald test, Benjamini–Hochberg, padj < 0.05 and |log2FC| > 1
#   * Hormone vs No hormone: Hydrog.Exp1.B3 and Hydrog.Exp1.B4 excluded (low depth)
#   * Embryo vs No embryo: Hydrog.Exp1.B3 excluded; the no-embryo controls are the
#     hormone-treated (B3) libraries, because embryos were cultured on hormone-
#     stimulated constructs.
# A sensitivity analysis adds experiment as a covariate (~ experiment + condition).
# =============================================================================

padj_cut <- 0.05
lfc_cut  <- 1
min_libsize <- 10e6   # gene-assigned counts; libraries below this are excluded

counts_all <- read_csv(file.path(dir_data, "rnaseq", "counts_gene_level.csv"), show_col_types = FALSE) |>
  column_to_rownames("gene") |> as.matrix()
samples <- read_csv(file.path(dir_data, "rnaseq", "sample_metadata.csv"), show_col_types = FALSE) |>
  as.data.frame()
rownames(samples) <- samples$library
stopifnot(all(samples$library %in% colnames(counts_all)))
counts_all <- counts_all[, samples$library]

# ---- library QC ---------------------------------------------------------------
samples$gene_assigned_counts <- colSums(counts_all)
samples$genes_detected       <- colSums(counts_all > 0)
samples$passes_depth         <- samples$gene_assigned_counts >= min_libsize
write_csv(samples, file.path(dir_tab, "rnaseq_library_qc.csv"))
log_result("RNA-seq library QC (threshold ", min_libsize / 1e6, " million gene-assigned counts):")
for (i in seq_len(nrow(samples))) {
  log_result(sprintf("  %s  %-10s  %s counts  %s genes  %s", samples$library[i], samples$condition[i],
                     format(samples$gene_assigned_counts[i], big.mark = ","),
                     format(samples$genes_detected[i], big.mark = ","),
                     ifelse(samples$passes_depth[i], "retained", "EXCLUDED (low depth)")))
}

keep_gene <- !grepl("^(LOC|NON)", rownames(counts_all))
counts_annot <- counts_all[keep_gene, ]

# ---- global view of ALL nine libraries (incl. excluded) -------------------------
dds_all <- DESeqDataSetFromMatrix(round(counts_annot[rowSums(counts_annot) >= 10, ]),
                                  colData = samples, design = ~ 1)
vsd_all <- vst(dds_all, blind = TRUE)
pca_df <- function(vsd, ntop = 500) {
  m <- assay(vsd)
  top <- head(order(apply(m, 1, var), decreasing = TRUE), ntop)
  pc <- prcomp(t(m[top, ]))
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  list(df = cbind(as.data.frame(colData(vsd)), PC1 = pc$x[, 1], PC2 = pc$x[, 2]), ve = ve)
}
pa <- pca_df(vsd_all)
pa$df$status <- ifelse(pa$df$passes_depth, "Retained", "Excluded (low depth)")
p_pca_all <- ggplot(pa$df, aes(PC1, PC2, colour = condition, shape = status)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = sub("Hydrog.", "", library)), size = 2.2, show.legend = FALSE) +
  scale_shape_manual(values = c("Retained" = 16, "Excluded (low depth)" = 4)) +
  scale_colour_manual(values = c("Embryo" = "#3B7FB6", "Hormone" = "#00BFC4", "No hormone" = "#F8766D")) +
  labs(x = paste0("PC1 (", pa$ve[1], "%)"), y = paste0("PC2 (", pa$ve[2], "%)"), colour = NULL, shape = NULL)
cor_all <- cor(assay(vsd_all), method = "pearson")
cor_long <- as.data.frame(as.table(cor_all)) |> setNames(c("a", "b", "r")) |>
  dplyr::mutate(a = factor(sub("Hydrog.", "", a), levels = sub("Hydrog.", "", samples$library)),
         b = factor(sub("Hydrog.", "", b), levels = rev(sub("Hydrog.", "", samples$library))))
p_cor_all <- ggplot(cor_long, aes(a, b, fill = r)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 1.8) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Pearson r") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
write_csv(rownames_to_column(as.data.frame(cor_all), "library"), file.path(dir_tab, "rnaseq_sample_correlation_all_libraries.csv"))
save_panel(p_pca_all, "SuppFig4a_PCA_all_libraries", 90, 75)
save_panel(p_cor_all, "SuppFig4b_correlation_all_libraries", 100, 85)

# ---- contrast helper ----------------------------------------------------------------
run_contrast <- function(libs, cond_levels, numerator, denominator, label) {
  cts <- counts_annot[, libs]
  cts <- cts[rowSums(cts) >= 10, ]
  cd  <- samples[libs, ]
  cd$cond_label <- factor(cd$condition, levels = cond_levels)            # display names
  cd$condition  <- factor(gsub(" ", "_", cd$condition), levels = gsub(" ", "_", cond_levels))
  cd$experiment <- factor(cd$experiment)
  dds <- DESeqDataSetFromMatrix(round(cts), colData = cd, design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", gsub(" ", "_", numerator), gsub(" ", "_", denominator)), alpha = padj_cut)
  tab <- as.data.frame(res) |> rownames_to_column("gene") |> dplyr::arrange(padj)
  tab$direction <- dplyr::case_when(
    !is.na(tab$padj) & tab$padj < padj_cut & tab$log2FoldChange >  lfc_cut ~ paste("Higher in", numerator),
    !is.na(tab$padj) & tab$padj < padj_cut & tab$log2FoldChange < -lfc_cut ~ paste("Lower in", numerator),
    TRUE ~ "Not significant")
  deg <- dplyr::filter(tab, direction != "Not significant")

  # sensitivity: experiment as covariate
  dds_b <- DESeqDataSetFromMatrix(round(cts), colData = cd, design = ~ experiment + condition)
  dds_b <- DESeq(dds_b, quiet = TRUE)
  res_b <- results(dds_b, contrast = c("condition", gsub(" ", "_", numerator), gsub(" ", "_", denominator)), alpha = padj_cut)
  tab_b <- as.data.frame(res_b) |> rownames_to_column("gene")
  deg_b <- tab_b$gene[!is.na(tab_b$padj) & tab_b$padj < padj_cut & abs(tab_b$log2FoldChange) > lfc_cut]
  same_sign <- sign(tab$log2FoldChange[match(deg$gene, tab$gene)]) ==
               sign(tab_b$log2FoldChange[match(deg$gene, tab_b$gene)])

  vsd <- vst(dds, blind = TRUE)
  norm <- counts(dds, normalized = TRUE)

  log_result("")
  log_result("== ", label, " (", numerator, " vs ", denominator, "; libraries: ", paste(libs, collapse = ", "), ")")
  log_result("  genes tested after filtering: ", nrow(cts), "; genes with non-NA padj: ", sum(!is.na(res$padj)))
  log_result("  DEGs: ", nrow(deg), " (", sum(deg$log2FoldChange > 0), " higher, ",
             sum(deg$log2FoldChange < 0), " lower in ", numerator, ")")
  log_result("  Sensitivity (~ experiment + condition): ", length(deg_b), " DEGs; ",
             sum(deg$gene %in% deg_b), " of the ", nrow(deg), " primary DEGs also significant; ",
             sum(same_sign, na.rm = TRUE), " of ", nrow(deg), " have the same direction of change")

  list(label = label, dds = dds, res = tab, deg = deg, res_batch = tab_b, deg_batch = deg_b,
       vsd = vsd, norm = norm, numerator = numerator, denominator = denominator, libs = libs)
}

hormone <- run_contrast(libs = c("Hydrog.Exp2.B3", "Hydrog.Exp3.B3", "Hydrog.Exp2.B4", "Hydrog.Exp3.B4"),
                        cond_levels = c("No hormone", "Hormone"),
                        numerator = "Hormone", denominator = "No hormone", label = "Hormone vs No hormone")

# embryo contrast: B3 libraries relabelled as "No embryo"
samples$condition_embryo <- ifelse(samples$embryo == "yes", "Embryo", "No embryo")
emb_libs <- c("Hydrog.Exp1.B2", "Hydrog.Exp2.B2", "Hydrog.Exp3.B2", "Hydrog.Exp2.B3", "Hydrog.Exp3.B3")
samples_bak <- samples
samples$condition[samples$library %in% emb_libs] <- samples$condition_embryo[samples$library %in% emb_libs]
embryo <- run_contrast(libs = emb_libs, cond_levels = c("No embryo", "Embryo"),
                       numerator = "Embryo", denominator = "No embryo", label = "Embryo vs No embryo")
samples <- samples_bak

# key genes quoted in the text
quote_genes <- list(
  hormone = c("GREB1", "SELENOP", "DHRS9", "HLTF", "TRPM3", "EIF4A2", "PTGS1", "WNT7B", "RGS3", "NR4A1", "IL11",
              "SLC6A12", "CD247", "RANBP3L", "LY6G6C", "LY6G6E", "LY6G6F", "POSTN", "LUM", "DCN", "FMOD", "EMB",
              "LOX", "HBA1", "HBA", "HBA2", "PADI1", "EPYC"),
  embryo  = c("RSAD2", "CMPK2", "ISG15", "IFI6", "IFI44", "IFI44L", "IFIT1", "IFIT2", "OAS1Z", "MX2", "STAT1",
              "TNFSF10", "ZBP1", "PLAC8", "PSMA8", "SFI1", "NRXN1", "KEH36_p01"))
for (nm in names(quote_genes)) {
  obj <- if (nm == "hormone") hormone else embryo
  q <- obj$res[match(quote_genes[[nm]], obj$res$gene), c("gene", "log2FoldChange", "padj", "direction")]
  q$gene <- quote_genes[[nm]]
  write_csv(q, file.path(dir_tab, paste0("rnaseq_genes_quoted_in_text_", nm, ".csv")))
  log_result("  Genes quoted in text (", nm, "): ",
             paste0(q$gene, " ", ifelse(is.na(q$log2FoldChange), "not tested",
                    paste0(sprintf("%+.2f", q$log2FoldChange), ifelse(q$direction == "Not significant", " (ns)", ""))),
                    collapse = "; "))
}

# ---- plots per contrast -------------------------------------------------------------
plot_contrast <- function(obj, cols, prefix) {
  num <- obj$numerator; den <- obj$denominator
  cols_dir <- setNames(c(cols[[num]], cols[[den]], "grey75"),
                       c(paste("Higher in", num), paste("Lower in", num), "Not significant"))
  # PCA
  pc <- pca_df(obj$vsd)
  pc$df$cond_label <- factor(pc$df$cond_label, levels = c(den, num))
  p_pca <- ggplot(pc$df, aes(PC1, PC2, colour = cond_label)) + geom_point(size = 2.5) +
    geom_text_repel(aes(label = paste0("Exp", sub(".*Exp(\\d).*", "\\1", library))), size = 2.2, box.padding = 0.5, point.padding = 0.4, min.segment.length = 0, show.legend = FALSE) +
    scale_colour_manual(values = cols) +
    labs(x = paste0("PC1 (", pc$ve[1], "%)"), y = paste0("PC2 (", pc$ve[2], "%)"), colour = NULL)
  # MA
  r <- obj$res |> dplyr::filter(!is.na(log2FoldChange))
  p_ma <- ggplot(r, aes(baseMean, log2FoldChange, colour = direction)) +
    geom_point(data = dplyr::filter(r, direction == "Not significant"), size = 0.3) +
    geom_point(data = dplyr::filter(r, direction != "Not significant"), size = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_x_log10(breaks = c(1e1, 1e3, 1e5), labels = c(expression(10^1), expression(10^3), expression(10^5))) + scale_colour_manual(values = cols_dir) +
    labs(x = "Mean normalised count", y = paste0("log2 fold change (", num, " / ", den, ")"), colour = NULL)
  # volcano
  rv <- dplyr::filter(r, !is.na(padj))
  lab <- head(dplyr::filter(rv, direction != "Not significant") |> dplyr::arrange(padj), 20)
  p_vol <- ggplot(rv, aes(log2FoldChange, -log10(padj), colour = direction)) +
    geom_point(data = dplyr::filter(rv, direction == "Not significant"), size = 0.3) +
    geom_point(data = dplyr::filter(rv, direction != "Not significant"), size = 0.8) +
    geom_hline(yintercept = -log10(padj_cut), linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, linewidth = 0.3) +
    geom_text_repel(data = lab, aes(label = gene), size = 2, max.overlaps = 30, fontface = "italic", show.legend = FALSE) +
    scale_colour_manual(values = cols_dir) +
    labs(x = paste0("log2 fold change (", num, " / ", den, ")"), y = expression(-log[10]~adjusted~p), colour = NULL)
  # heatmap of up to 100 DEGs (ranked by padj), z-scored VST, genes clustered (Euclidean, Ward.D2)
  g <- head(obj$deg$gene, 100)
  m <- assay(obj$vsd)[g, , drop = FALSE]
  z <- t(scale(t(m)))
  ord <- if (nrow(z) > 2) hclust(dist(z), method = "ward.D2")$order else seq_len(nrow(z))
  cd <- as.data.frame(colData(obj$vsd))
  sample_order <- rownames(cd)[order(factor(cd$cond_label, levels = c(num, den)))]
  hl <- data.frame(gene = rep(rownames(z), times = ncol(z)),
                   library = rep(colnames(z), each = nrow(z)),
                   z = as.vector(z), stringsAsFactors = FALSE)
  hl$gene <- factor(hl$gene, levels = rownames(z)[ord])
  hl$library <- factor(hl$library, levels = sample_order)
  p_hm <- ggplot(hl, aes(library, gene, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "Row z-score") +
    scale_x_discrete(labels = function(x) paste0(cd$cond_label[match(x, rownames(cd))], "\nExp", sub(".*Exp(\\d).*", "\\1", x))) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_text(size = if (length(g) > 60) 3.2 else 4.5, face = "italic"),
          axis.text.x = element_text(size = 6, angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
  save_panel(p_pca, paste0(prefix, "_PCA"), 70, 60)
  save_panel(p_ma, paste0(prefix, "_MA"), 90, 60)
  save_panel(p_vol, paste0(prefix, "_volcano"), 90, 70)
  save_panel(p_hm, paste0(prefix, "_heatmap"), 70, 150)
  list(pca = p_pca, ma = p_ma, volcano = p_vol, heatmap = p_hm, pca_var = pc$ve)
}
plots_hormone <- plot_contrast(hormone, col_hormone, "Fig4_hormone")
plots_embryo  <- plot_contrast(embryo,  col_embryo,  "Fig5_embryo")
log_result("  PCA variance explained — hormone: PC1 ", plots_hormone$pca_var[1], "%, PC2 ", plots_hormone$pca_var[2],
           "%; embryo: PC1 ", plots_embryo$pca_var[1], "%, PC2 ", plots_embryo$pca_var[2], "%")

# ---- Supplementary Data 1 and 2 (Excel) -----------------------------------------------
write_supp_data <- function(obj, file, title) {
  wb <- createWorkbook()
  addWorksheet(wb, "README")
  writeData(wb, "README", data.frame(Item = c(
    title,
    paste0("Comparison: ", obj$numerator, " vs ", obj$denominator, " (", obj$numerator, " is the numerator)"),
    paste0("Positive log2FoldChange = higher expression in ", obj$numerator, "; negative = lower expression in ", obj$numerator),
    paste0("Libraries: ", paste(obj$libs, collapse = ", ")),
    "DESeq2 design: ~ condition; Wald test; Benjamini-Hochberg adjusted p-values (padj)",
    "Genes with identifiers starting with LOC or NON removed; genes with < 10 total counts across these libraries removed",
    "Differentially expressed genes (DEGs): padj < 0.05 and |log2FoldChange| > 1",
    "Sheet 'Sensitivity_experiment' reports the same contrast with design ~ experiment + condition")))
  addWorksheet(wb, "DEGs");                   writeData(wb, "DEGs", obj$deg)
  addWorksheet(wb, "All_genes");              writeData(wb, "All_genes", obj$res)
  addWorksheet(wb, "Sensitivity_experiment"); writeData(wb, "Sensitivity_experiment", dplyr::arrange(obj$res_batch, padj))
  addWorksheet(wb, "Normalised_counts")
  writeData(wb, "Normalised_counts", rownames_to_column(as.data.frame(round(obj$norm, 2)), "gene"))
  saveWorkbook(wb, file, overwrite = TRUE)
}
write_supp_data(hormone, file.path(dir_tab, "Supplementary_Data_1_DEGs_Hormone_vs_NoHormone.xlsx"),
                "Supplementary Data 1. Differential expression, hormone-treated vs untreated co-cultures")
write_supp_data(embryo,  file.path(dir_tab, "Supplementary_Data_2_DEGs_Embryo_vs_NoEmbryo.xlsx"),
                "Supplementary Data 2. Differential expression, embryo-exposed vs non-exposed co-cultures")
saveRDS(list(hormone = hormone[c("res", "deg", "res_batch")], embryo = embryo[c("res", "deg", "res_batch")]),
        file.path(dir_tab, "rnaseq_results.rds"))

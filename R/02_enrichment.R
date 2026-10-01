# =============================================================================
# 02_enrichment.R — Gene Ontology (biological process) over-representation of the
# genes with higher and lower expression in each contrast (clusterProfiler, Bos taurus).
# Universe = all genes tested in the respective contrast.
# =============================================================================
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Bt.eg.db) })

run_go <- function(obj, prefix) {
  to_entrez <- function(sym) {
    m <- suppressWarnings(bitr(unique(sym), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Bt.eg.db))
    unique(m$ENTREZID)
  }
  universe <- to_entrez(obj$res$gene[!is.na(obj$res$padj)])
  sets <- list(higher = obj$deg$gene[obj$deg$log2FoldChange > 0],
               lower  = obj$deg$gene[obj$deg$log2FoldChange < 0])
  out <- list()
  for (s in names(sets)) {
    ids <- to_entrez(sets[[s]])
    ego <- if (length(ids) >= 3) enrichGO(ids, universe = universe, OrgDb = org.Bt.eg.db, ont = "BP",
                                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                                          minGSSize = 10, maxGSSize = 500, readable = TRUE) else NULL
    tab <- if (is.null(ego)) data.frame() else as.data.frame(ego)
    lab <- paste0(ifelse(s == "higher", "Higher", "Lower"), " in ", obj$numerator)
    write_csv(tab, file.path(dir_tab, paste0(prefix, "_GO_BP_", s, "_in_", gsub(" ", "", obj$numerator), ".csv")))
    log_result("  GO BP (", obj$label, ", ", lab, "): ", length(sets[[s]]), " genes, ", length(ids),
               " mapped to Entrez; ", nrow(tab), " terms with padj < 0.05",
               if (nrow(tab)) paste0(". Top: ", paste(head(tab$Description, 8), collapse = "; ")) else "")
    if (nrow(tab)) out[[lab]] <- head(tab, 15) |> dplyr::mutate(set = lab)
  }
  if (!length(out)) return(NULL)
  d <- dplyr::bind_rows(out) |>
    dplyr::mutate(GeneRatioNum = vapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)),
           Description = factor(Description, levels = unique(rev(Description))))
  p <- ggplot(d, aes(GeneRatioNum, Description, size = Count, colour = p.adjust)) + geom_point() +
    facet_grid(set ~ ., scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#B2182B", high = "#2166AC", name = "Adjusted p") +
    labs(x = "Gene ratio", y = NULL, size = "Genes") + theme(axis.text.y = element_text(size = 6))
  save_panel(p, prefix, 160, max(60, 6 * nrow(d) + 25))
  p
}
log_result("")
log_result("== Gene Ontology enrichment")
go_hormone <- run_go(hormone, "SuppTable_GO_hormone")
go_embryo  <- run_go(embryo,  "SuppTable_GO_embryo")

# =============================================================================
# Reactome pathway over-representation (as in the original submission).
# Reactome has no bovine annotation, so bovine gene symbols are mapped to their human
# orthologues by gene symbol (org.Hs.eg.db) and tested with ReactomePA (human pathways).
# Universe = all genes tested in the respective contrast.
# =============================================================================
suppressPackageStartupMessages({ library(ReactomePA); library(org.Hs.eg.db) })
run_reactome <- function(obj, prefix) {
  to_hs <- function(sym) {
    m <- suppressWarnings(bitr(unique(sym), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    unique(m$ENTREZID)
  }
  universe <- to_hs(obj$res$gene[!is.na(obj$res$padj)])
  sets <- list(higher = obj$deg$gene[obj$deg$log2FoldChange > 0],
               lower  = obj$deg$gene[obj$deg$log2FoldChange < 0])
  out <- list()
  for (s in names(sets)) {
    ids <- to_hs(sets[[s]])
    er <- if (length(ids) >= 3) enrichPathway(ids, organism = "human", universe = universe, pAdjustMethod = "BH",
                                              pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 10,
                                              maxGSSize = 500, readable = TRUE) else NULL
    tab <- if (is.null(er)) data.frame() else as.data.frame(er)
    lab <- paste0(ifelse(s == "higher", "Higher", "Lower"), " in ", obj$numerator)
    write_csv(tab, file.path(dir_tab, paste0(prefix, "_Reactome_", s, "_in_", gsub(" ", "", obj$numerator), ".csv")))
    log_result("  Reactome (", obj$label, ", ", lab, "): ", length(sets[[s]]), " genes, ", length(ids),
               " mapped to human Entrez; ", nrow(tab), " pathways with padj < 0.05",
               if (nrow(tab)) paste0(". Top: ", paste(head(tab$Description, 8), collapse = "; ")) else "")
    if (nrow(tab)) out[[lab]] <- head(tab, 15) |> dplyr::mutate(set = lab)
  }
  if (!length(out)) return(NULL)
  d <- dplyr::bind_rows(out) |>
    dplyr::mutate(GeneRatioNum = vapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)),
                  Description = factor(Description, levels = unique(rev(Description))))
  p <- ggplot(d, aes(GeneRatioNum, Description, size = Count, colour = p.adjust)) + geom_point() +
    facet_grid(set ~ ., scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#B2182B", high = "#2166AC", name = "Adjusted p") +
    labs(x = "Gene ratio", y = NULL, size = "Genes") + theme(axis.text.y = element_text(size = 6))
  save_panel(p, prefix, 160, max(60, 6 * nrow(d) + 25))
  p
}
log_result("")
log_result("== Reactome pathway enrichment (human orthologues by gene symbol)")
re_hormone <- run_reactome(hormone, "SuppFig2_Reactome_hormone")
re_embryo  <- run_reactome(embryo,  "SuppFig3_Reactome_embryo")

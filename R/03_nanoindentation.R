# =============================================================================
# 03_nanoindentation.R — Young's modulus (Figure 3b) and viscoelasticity (Figure 3c–e)
#
# Unit structure: indentation sites are nested in hydrogels (Batch x Replicate x Drop x
# condition), hydrogels are nested in independent preparations (EndoECM batch).
# Models: linear mixed-effects models on log-transformed moduli with preparation and
# hydrogel as random intercepts; effects are reported as ratios (back-transformed).
# =============================================================================

ym_raw <- read_csv(file.path(dir_data, "nanoindentation", "young_modulus.csv"), show_col_types = FALSE)
assign_group <- function(cells, horm) dplyr::case_when(
  horm == "yTissue" ~ "Native tissue",
  cells == "no" ~ "No cells",
  cells == "yes" & horm == "no" ~ "Cells",
  cells == "yes" & horm == "yes" ~ "Cells + hormones")
ym <- ym_raw |>
  dplyr::mutate(group = factor(assign_group(Cells, Hormonal_treatment), levels = names(col_group)),
         Day = factor(Day, levels = c(-3, 1, 7)),
         prep = factor(Batch),
         gel = factor(paste(Batch, Replicate, Drop, group, sep = "_")))
tissue_ym <- dplyr::filter(ym, group == "Native tissue")
gels <- dplyr::filter(ym, group != "Native tissue") |> droplevels()
gels$cond <- interaction(gels$group, gels$Day, sep = "|", drop = TRUE)

log_result("")
log_result("== Young's modulus")
log_result("  ", nrow(gels), " indentation sites from ", dplyr::n_distinct(gels$gel), " hydrogels in ",
           dplyr::n_distinct(gels$prep), " independent preparations; native tissue: ", nrow(tissue_ym), " sites")
cnt <- gels |> dplyr::group_by(Day, group) |> dplyr::summarise(sites = n(), hydrogels = dplyr::n_distinct(gel),
                                                 preparations = dplyr::n_distinct(prep), mean_YM = mean(YM),
                                                 median_YM = median(YM), .groups = "drop")
write_csv(cnt, file.path(dir_tab, "young_modulus_design_and_means.csv"))

m_ym <- lmer(log(YM) ~ cond + (1 | prep) + (1 | gel), data = gels)
emm <- emmeans(m_ym, ~ cond)
build_contrasts <- function(dat) {
  lv <- levels(dat$cond)
  mk <- function(a, b) { v <- setNames(rep(0, length(lv)), lv); v[a] <- 1; v[b] <- -1; v }
  out <- list()
  for (dy in levels(dat$Day)) {
    present <- intersect(c("No cells", "Cells", "Cells + hormones"), as.character(unique(dat$group[dat$Day == dy])))
    if (length(present) < 2) next
    for (pr in combn(present, 2, simplify = FALSE))
      out[[paste0("Day ", dy, ": ", pr[2], " / ", pr[1])]] <- mk(paste0(pr[2], "|", dy), paste0(pr[1], "|", dy))
  }
  for (g in c("No cells", "Cells", "Cells + hormones")) {
    dys <- as.character(unique(dat$Day[dat$group == g]))
    if (all(c("-3", "7") %in% dys)) out[[paste0(g, ": Day 7 / Day -3")]] <- mk(paste0(g, "|7"), paste0(g, "|-3"))
    if (all(c("1", "7") %in% dys))  out[[paste0(g, ": Day 7 / Day 1")]]  <- mk(paste0(g, "|7"), paste0(g, "|1"))
  }
  out
}
con_list <- build_contrasts(gels)
ct <- summary(contrast(emm, con_list, adjust = "holm"), infer = TRUE, type = "response")
ct <- as.data.frame(ct)
write_csv(ct, file.path(dir_tab, "young_modulus_contrasts.csv"))
for (i in seq_len(nrow(ct))) log_result(sprintf("  %-40s ratio %.2f (95%% CI %.2f–%.2f), Holm p = %s",
  as.character(ct$contrast[i]), ct$ratio[i], ct$lower.CL[i], ct$upper.CL[i], fmt_p(ct$p.value[i])))

# cell x day interaction (no-hormone hydrogels only), likelihood-ratio test
nh <- dplyr::filter(gels, group %in% c("No cells", "Cells")) |> droplevels()
m1 <- lmer(log(YM) ~ group * Day + (1 | prep) + (1 | gel), data = nh, REML = FALSE)
m0 <- lmer(log(YM) ~ group + Day + (1 | prep) + (1 | gel), data = nh, REML = FALSE)
lrt <- anova(m0, m1)
log_result(sprintf("  Cells x Day interaction (no-hormone hydrogels): chi2(%d) = %.2f, p = %s",
                   lrt$Df[2], lrt$Chisq[2], fmt_p(lrt$`Pr(>Chisq)`[2])))
log_result(sprintf("  Native tissue: mean %.0f Pa, median %.0f Pa (n = %d sites)", mean(tissue_ym$YM), median(tissue_ym$YM), nrow(tissue_ym)))

# ---- Figure 3b ----------------------------------------------------------------------
plot_ym <- dplyr::bind_rows(gels, tissue_ym) |> dplyr::mutate(Day = factor(paste("Day", Day), levels = paste("Day", c(-3, 1, 7))))
sig <- ct |> dplyr::filter(p.value < 0.05, grepl("^Day", contrast)) |>
  dplyr::mutate(Day = paste("Day", sub("^Day (-?\\d+):.*", "\\1", contrast)),
         g2 = sub("^Day -?\\d+: (.*) / .*$", "\\1", contrast), g1 = sub("^.* / ", "", contrast))
means <- plot_ym |> dplyr::group_by(Day, group) |> dplyr::summarise(m = mean(YM), .groups = "drop")
p_ym <- ggplot(plot_ym, aes(group, YM / 1000, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.3, width = 0.7) +
  geom_point(aes(colour = group), position = position_jitter(width = 0.15, height = 0), size = 0.5, alpha = 0.7) +
  geom_text(data = means, aes(y = m / 1000, label = sprintf("%.2f", m / 1000)), size = 2, vjust = -0.6, fontface = "bold") +
  facet_grid(. ~ Day, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = col_group) + scale_colour_manual(values = col_group) +
  labs(x = NULL, y = "Young's modulus (kPa)", fill = NULL, colour = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
if (nrow(sig)) {
  ymax <- max(plot_ym$YM) / 1000
  sig <- sig |> dplyr::group_by(Day) |> dplyr::mutate(y = ymax * (1.05 + 0.08 * (row_number() - 1))) |> ungroup()
  p_ym <- p_ym + geom_segment(data = sig, aes(x = g1, xend = g2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = sig, aes(x = g1, y = y, label = p_label(p.value)), inherit.aes = FALSE, size = 2, vjust = -0.3, hjust = -0.1)
}
save_panel(p_ym, "Fig3b_young_modulus", 120, 75)

# ---- Viscoelasticity (E', E'') -----------------------------------------------------------
ve_raw <- read_csv(file.path(dir_data, "nanoindentation", "viscoelasticity.csv"), show_col_types = FALSE)
key <- c("Batch", "Replicate", "Drop", "Date", "Day", "Cells", "Hormonal_treatment")
ve <- ve_raw |>
  dplyr::group_by(across(all_of(c(key, "Module", "Hz")))) |> dplyr::mutate(site = row_number()) |> ungroup() |>
  pivot_wider(id_cols = all_of(c(key, "site", "Hz")), names_from = Module, values_from = Indentation) |>
  dplyr::mutate(group = factor(assign_group(Cells, Hormonal_treatment), levels = names(col_group)),
         Day = factor(Day, levels = c(-3, 1, 7)), prep = factor(Batch),
         gel = factor(paste(Batch, Replicate, Drop, group, sep = "_")))

# tan(delta) per site = mean over frequencies of E''/E'; sites with non-positive moduli excluded
tan_site <- ve |> dplyr::filter(E1 > 0, E2 > 0) |>
  dplyr::group_by(across(all_of(c(key, "site", "group", "prep", "gel")))) |>
  dplyr::summarise(tan_delta = mean(E2 / E1), n_hz = n(), .groups = "drop") |> dplyr::filter(n_hz == 4) |> droplevels()
tan_site$cond <- interaction(tan_site$group, tan_site$Day, sep = "|", drop = TRUE)
m_tan <- lmer(log(tan_delta) ~ cond + (1 | prep) + (1 | gel), data = tan_site)
emm_t <- emmeans(m_tan, ~ cond)
ct_t <- as.data.frame(summary(contrast(emm_t, build_contrasts(tan_site), adjust = "holm"), infer = TRUE, type = "response"))
write_csv(ct_t, file.path(dir_tab, "tan_delta_contrasts.csv"))
tan_means <- tan_site |> dplyr::group_by(Day, group) |> dplyr::summarise(mean_tan_delta = mean(tan_delta), sites = n(),
                                                          hydrogels = dplyr::n_distinct(gel), .groups = "drop")
write_csv(tan_means, file.path(dir_tab, "tan_delta_means.csv"))
log_result("")
log_result("== Viscoelasticity: tan(delta) = E''/E' (mean over 1–20 Hz); all conditions < 1: ",
           all(tan_means$mean_tan_delta < 1))
for (i in seq_len(nrow(tan_means))) log_result(sprintf("  Day %s %-17s mean tan(delta) %.3f (%d sites, %d hydrogels)",
  as.character(tan_means$Day[i]), as.character(tan_means$group[i]), tan_means$mean_tan_delta[i], tan_means$sites[i], tan_means$hydrogels[i]))
for (i in seq_len(nrow(ct_t))) log_result(sprintf("  %-40s ratio %.2f (95%% CI %.2f–%.2f), Holm p = %s",
  as.character(ct_t$contrast[i]), ct_t$ratio[i], ct_t$lower.CL[i], ct_t$upper.CL[i], fmt_p(ct_t$p.value[i])))

# Figure 3c–e: E' and E'' vs frequency; mean ± 95% CI across hydrogel means
gel_means <- ve |> pivot_longer(c(E1, E2), names_to = "Module", values_to = "value") |>
  dplyr::group_by(group, Day, prep, gel, Hz, Module) |> dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
curve <- gel_means |> dplyr::group_by(group, Day, Hz, Module) |>
  dplyr::summarise(n = n(), mean = mean(value), se = sd(value) / sqrt(n),
            lo = mean - qt(0.975, pmax(n - 1, 1)) * se, hi = mean + qt(0.975, pmax(n - 1, 1)) * se, .groups = "drop") |>
  dplyr::mutate(Module = dplyr::recode(Module, E1 = "E'", E2 = "E''"), Day = factor(paste("Day", Day), levels = paste("Day", c(-3, 1, 7))))
tissue_file <- file.path(dir_data, "nanoindentation", "native_tissue_viscoelasticity.csv")
if (file.exists(tissue_file)) {
  tv <- read_csv(tissue_file, show_col_types = FALSE) |>   # native endometrium (ref. 15): Animal_ID, Estrus_cycle, Module, Hz, Indentation
    dplyr::group_by(Animal_ID, Hz, Module) |> dplyr::summarise(value = mean(Indentation, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(Hz, Module) |> dplyr::summarise(n = n(), mean = mean(value), se = sd(value) / sqrt(n),
                                       lo = mean - qt(0.975, n - 1) * se, hi = mean + qt(0.975, n - 1) * se, .groups = "drop") |>
    dplyr::mutate(Module = dplyr::recode(Module, E1 = "E'", E2 = "E''"), Day = "Native tissue")
} else tv <- NULL
ve_panel <- function(g) {
  d <- dplyr::filter(curve, group == g)
  if (!is.null(tv)) d <- dplyr::bind_rows(d, tv)
  d$Day <- factor(d$Day, levels = c(paste("Day", c(-3, 1, 7)), "Native tissue"))
  ggplot(d, aes(Hz, mean / 1000, colour = Day, fill = Day, linetype = Module)) +
    geom_ribbon(aes(ymin = lo / 1000, ymax = hi / 1000), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.5) + geom_point(size = 0.8) +
    scale_colour_manual(values = c("Day -3" = "#E7298A", "Day 1" = "#A6761D", "Day 7" = "#1B9E77", "Native tissue" = "#7B3294")) +
    scale_fill_manual(values = c("Day -3" = "#E7298A", "Day 1" = "#A6761D", "Day 7" = "#1B9E77", "Native tissue" = "#7B3294")) +
    scale_x_continuous(breaks = c(1, 5, 10, 20)) +
    labs(title = g, x = "Frequency (Hz)", y = "E' / E'' (kPa)", colour = NULL, fill = NULL, linetype = NULL)
}
p_ve <- lapply(c("No cells", "Cells", "Cells + hormones"), ve_panel)
names(p_ve) <- c("No cells", "Cells", "Cells + hormones")
for (g in names(p_ve)) save_panel(p_ve[[g]], paste0("Fig3_viscoelastic_", gsub("[^A-Za-z]", "", g)), 80, 60)
if (is.null(tv)) log_result("  NOTE: native-tissue E'/E'' curve not plotted (data/nanoindentation/native_tissue_viscoelasticity.csv not supplied)")

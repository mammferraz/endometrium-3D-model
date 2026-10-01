# =============================================================================
# 04_gland_area.R — gland area under the hormone regimen (Figure 4a) and with embryos (Figure 5a)
#
# Individual glands are subsamples nested in hydrogels (wells); hydrogels are nested in
# independent co-culture experiments (Replicate). Linear mixed-effects models on log(area)
# with experiment and hydrogel as random intercepts; effects reported as ratios.
# =============================================================================

ga <- read_csv(file.path(dir_data, "gland_area", "gland_area.csv"), show_col_types = FALSE) |>
  dplyr::mutate(experiment = factor(Replicate), hydrogel = factor(paste(Replicate, Well_group, sep = "_")))

# ---- Hormone vs no hormone (embryo-free wells), all imaging days ---------------------------
hd <- ga |> dplyr::filter(Embryos == "No") |>
  dplyr::mutate(treatment = factor(ifelse(Hormones == "Yes", "Hormone", "No hormone"), levels = c("No hormone", "Hormone")),
         Day = factor(Day_exp, levels = sort(unique(Day_exp))))
log_result("")
log_result("== Gland area, hormone vs no hormone")
log_result("  ", nrow(hd), " gland measurements from ", dplyr::n_distinct(hd$hydrogel), " hydrogels in ",
           dplyr::n_distinct(hd$experiment), " independent experiments; imaging days: ", paste(levels(hd$Day), collapse = ", "))
write_csv(hd |> dplyr::count(experiment, Day, treatment, name = "glands"), file.path(dir_tab, "gland_area_hormone_design.csv"))

mh1 <- lmer(log(Area_um2) ~ treatment * Day + (1 | experiment) + (1 | hydrogel), data = hd, REML = FALSE)
mh0 <- lmer(log(Area_um2) ~ treatment + Day + (1 | experiment) + (1 | hydrogel), data = hd, REML = FALSE)
lrt <- anova(mh0, mh1)
log_result(sprintf("  Treatment x day interaction: chi2(%d) = %.2f, p = %s", lrt$Df[2], lrt$Chisq[2], fmt_p(lrt$`Pr(>Chisq)`[2])))
mh <- lmer(log(Area_um2) ~ treatment * Day + (1 | experiment) + (1 | hydrogel), data = hd)
ct_h <- as.data.frame(summary(contrast(emmeans(mh, ~ treatment | Day), method = "revpairwise", adjust = "none"),
                              infer = TRUE, type = "response"))
ct_h$p.holm <- p.adjust(ct_h$p.value, "holm")
write_csv(ct_h, file.path(dir_tab, "gland_area_hormone_contrasts_by_day.csv"))
for (i in seq_len(nrow(ct_h))) log_result(sprintf("  Day %s: Hormone / No hormone area ratio %.3f (95%% CI %.3f–%.3f), p = %s, Holm p = %s",
  as.character(ct_h$Day[i]), ct_h$ratio[i], ct_h$lower.CL[i], ct_h$upper.CL[i], fmt_p(ct_h$p.value[i]), fmt_p(ct_h$p.holm[i])))

p_area_h <- ggplot(hd, aes(Day, Area_um2 / 1000, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.3, position = position_dodge(0.8), width = 0.7) +
  geom_point(aes(colour = treatment), position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.8),
             size = 0.3, alpha = 0.5) +
  scale_fill_manual(values = col_hormone) + scale_colour_manual(values = col_hormone) +
  labs(x = "Day of culture", y = expression("Gland area (" * 10^3 ~ mu * m^2 * ")"), fill = NULL, colour = NULL) +
  theme(legend.position = "top")
sig_h <- dplyr::filter(ct_h, p.holm < 0.05)
if (nrow(sig_h)) {
  ytop <- quantile(hd$Area_um2, 0.995) / 1000
  p_area_h <- p_area_h + geom_text(data = sig_h, aes(x = Day, y = ytop, label = p_label(p.holm)), inherit.aes = FALSE, size = 2)
}
save_panel(p_area_h, "Fig4a_gland_area_hormone", 120, 70)

# ---- Embryo vs no embryo (hormone-treated wells), Days 8–10 --------------------------------
ed <- ga |> dplyr::filter(Hormones == "Yes", Day_exp %in% c(8, 9, 10)) |>
  dplyr::mutate(treatment = factor(ifelse(Embryos == "Yes", "Embryo", "No embryo"), levels = c("No embryo", "Embryo")),
         Day = factor(Day_exp))
log_result("")
log_result("== Gland area, embryo vs no embryo (hormone-treated co-cultures)")
log_result("  ", nrow(ed), " gland measurements from ", dplyr::n_distinct(ed$hydrogel), " hydrogels in ", dplyr::n_distinct(ed$experiment), " experiments")
me1 <- lmer(log(Area_um2) ~ treatment * Day + (1 | experiment) + (1 | hydrogel), data = ed, REML = FALSE)
me0 <- lmer(log(Area_um2) ~ treatment + Day + (1 | experiment) + (1 | hydrogel), data = ed, REML = FALSE)
me00 <- lmer(log(Area_um2) ~ Day + (1 | experiment) + (1 | hydrogel), data = ed, REML = FALSE)
lrt_i <- anova(me0, me1); lrt_m <- anova(me00, me0)
log_result(sprintf("  Embryo x day interaction: chi2(%d) = %.2f, p = %s; embryo main effect: chi2(%d) = %.2f, p = %s",
                   lrt_i$Df[2], lrt_i$Chisq[2], fmt_p(lrt_i$`Pr(>Chisq)`[2]), lrt_m$Df[2], lrt_m$Chisq[2], fmt_p(lrt_m$`Pr(>Chisq)`[2])))
me <- lmer(log(Area_um2) ~ treatment * Day + (1 | experiment) + (1 | hydrogel), data = ed)
ct_e <- as.data.frame(summary(contrast(emmeans(me, ~ treatment | Day), method = "revpairwise", adjust = "none"),
                              infer = TRUE, type = "response"))
ct_e$p.holm <- p.adjust(ct_e$p.value, "holm")
write_csv(ct_e, file.path(dir_tab, "gland_area_embryo_contrasts_by_day.csv"))
for (i in seq_len(nrow(ct_e))) log_result(sprintf("  Day %s: Embryo / No embryo area ratio %.3f (95%% CI %.3f–%.3f), Holm p = %s",
  as.character(ct_e$Day[i]), ct_e$ratio[i], ct_e$lower.CL[i], ct_e$upper.CL[i], fmt_p(ct_e$p.holm[i])))

p_area_e <- ggplot(ed, aes(Day, Area_um2 / 1000, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.3, position = position_dodge(0.8), width = 0.7) +
  geom_point(aes(colour = treatment), position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.8),
             size = 0.3, alpha = 0.5) +
  scale_fill_manual(values = col_embryo) + scale_colour_manual(values = col_embryo) +
  labs(x = "Day of culture", y = expression("Gland area (" * 10^3 ~ mu * m^2 * ")"), fill = NULL, colour = NULL) +
  theme(legend.position = "top")
save_panel(p_area_e, "Fig5a_gland_area_embryo", 80, 70)

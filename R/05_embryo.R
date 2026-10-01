# =============================================================================
# 05_embryo.R — blastocyst hatching (Figure 5b)
#
# Hatching is a binomial outcome (hatched / allocated embryos per group and IVF replicate).
# Primary model: binomial GLMM with IVF replicate as random intercept; if the replicate
# variance is not estimable (singular fit), IVF replicate is entered as a fixed blocking factor.
# Pairwise odds ratios with Wald 95% CI, Holm-adjusted.
# =============================================================================

hat <- read_csv(file.path(dir_data, "embryo", "hatching.csv"), show_col_types = FALSE) |>
  dplyr::mutate(group = factor(dplyr::recode(Group, Cell = "Co-culture", No_Cell = "Cell-free hydrogel", Well_2D = "2D well"),
                        levels = c("2D well", "Cell-free hydrogel", "Co-culture")),
         replicate = factor(Replicate), not_hatched = Total - Hatched, pct = 100 * Hatched / Total)
log_result("")
log_result("== Blastocyst hatching (hatched/total per IVF replicate)")
for (g in levels(hat$group)) {
  d <- dplyr::filter(hat, group == g)
  log_result(sprintf("  %-19s %s; total %d/%d (%.1f%%)", g, paste0(d$Hatched, "/", d$Total, collapse = ", "),
                     sum(d$Hatched), sum(d$Total), 100 * sum(d$Hatched) / sum(d$Total)))
}
glmm <- tryCatch(glmer(cbind(Hatched, not_hatched) ~ group + (1 | replicate), data = hat, family = binomial), error = function(e) NULL)
use_glmm <- !is.null(glmm) && !isSingular(glmm)
if (use_glmm) {
  m_hat <- glmm
  m_hat0 <- glmer(cbind(Hatched, not_hatched) ~ 1 + (1 | replicate), data = hat, family = binomial)
  lrt <- anova(m_hat0, m_hat); chi <- lrt$Chisq[2]; dfc <- lrt$Df[2]; pc <- lrt$`Pr(>Chisq)`[2]
  model_txt <- "binomial GLMM, IVF replicate as random intercept"
} else {
  m_hat  <- glm(cbind(Hatched, not_hatched) ~ group + replicate, data = hat, family = binomial)
  m_hat0 <- glm(cbind(Hatched, not_hatched) ~ replicate, data = hat, family = binomial)
  lrt <- anova(m_hat0, m_hat, test = "LRT"); chi <- lrt$Deviance[2]; dfc <- lrt$Df[2]; pc <- lrt$`Pr(>Chi)`[2]
  model_txt <- "binomial GLM, IVF replicate as fixed blocking factor (random-effect variance not estimable)"
}
log_result("  Model: ", model_txt)
log_result(sprintf("  Treatment effect: likelihood-ratio chi2(%d) = %.2f, p = %s", dfc, chi, fmt_p(pc)))
or <- as.data.frame(summary(contrast(emmeans(m_hat, ~ group), method = "revpairwise", adjust = "holm"),
                            infer = TRUE, type = "response"))
write_csv(or, file.path(dir_tab, "hatching_odds_ratios.csv"))
for (i in seq_len(nrow(or))) log_result(sprintf("  %-36s OR %.2f (95%% CI %.2f–%.2f), Holm p = %s",
  as.character(or$contrast[i]), or$odds.ratio[i], or$asymp.LCL[i], or$asymp.UCL[i], fmt_p(or$p.value[i])))

sig_hat <- or |> dplyr::filter(p.value < 0.05) |>
  dplyr::mutate(contrast = gsub("[()]", "", as.character(contrast)), g2 = sub(" / .*", "", contrast), g1 = sub(".* / ", "", contrast), y = 100 + 8 * (row_number() - 1))
p_hat <- ggplot(hat, aes(group, pct, fill = group)) +
  stat_summary(fun = mean, geom = "col", alpha = 0.6, width = 0.6) +
  geom_point(aes(shape = replicate), size = 1.8, position = position_dodge(width = 0.35)) +
  scale_fill_manual(values = col_hatch, guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), name = "IVF replicate") +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 20)) +
  labs(x = NULL, y = "Hatched blastocysts (%)") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
if (nrow(sig_hat)) p_hat <- p_hat +
  geom_segment(data = sig_hat, aes(x = g1, xend = g2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
  geom_text(data = sig_hat, aes(x = g1, y = y, label = p_label(p.value)), inherit.aes = FALSE, size = 2, vjust = -0.3, hjust = -0.2)
save_panel(p_hat, "Fig5b_hatching", 75, 70)


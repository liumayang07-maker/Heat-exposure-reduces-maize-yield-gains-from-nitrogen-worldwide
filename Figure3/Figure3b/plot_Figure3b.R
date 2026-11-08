#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table); library(metafor); library(ggplot2)})
entry <- grep("^--file=", commandArgs(), value = TRUE)
out <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))

# 1. 从当前拟合结果重算环境等权预测；不在区域实际热暴露范围外外推。
bundle <- readRDS(file.path(out, "model_REML.rds"))
fit <- bundle$fit
form <- bundle$formula
d <- fread(file.path(out, "model_input.csv"))
support <- fread(file.path(out, "regional_support.csv"))
regional_tests <- fread(file.path(out, "regional_heat_N_joint_tests_REML.csv"))
transformations <- fread(file.path(out, "transformation_parameters.csv"))
region_levels <- c(
  "North America_Irrigated", "North America_Rainfed", "Europe", "WANA",
  "Northwest China", "Northeast China", "South Asia", "North China Plain",
  "Southwest China", "East Africa", "West Africa", "South America"
)
d[, region := relevel(factor(region, levels = region_levels), "North China Plain")]
d[, irrigation := factor(irrigation, levels = c("Rainfed", "Irrigation"))]
d[, N_strategy := factor(N_strategy, levels = c("S", "SV", "SVR", "SR", "VR", "V", "R"))]
N_center <- transformations[target == "N_c", center]
stopifnot(length(N_center) == 1L, transformations[target == "N_c", divisor] == 100)
beta <- coef(fit)
covariance <- vcov(fit)
predictions <- list()
for (reg in support$region) {
  ds <- copy(d[as.character(region) == reg])
  env_values <- unique(ds[, .(env_ID, heat_log)])
  heat_grid <- seq(min(env_values$heat_log), max(env_values$heat_log), length.out = 101)
  n_values <- support[region == reg, c(N_q25, N_q50, N_q75)]
  ds[, rows_per_environment := .N, by = env_ID]
  ds[, marginal_weight := 1 / rows_per_environment]
  ds[, marginal_weight := marginal_weight / sum(marginal_weight)]
  for (scenario in seq_along(n_values)) {
    n_value <- as.numeric(n_values[scenario])
    scenario_N_c <- (n_value - N_center) / 100
    for (heat_value in heat_grid) {
      nd <- copy(ds)
      nd[, `:=`(heat_wr = heat_value - heat_region_mean, N_c = scenario_N_c, N_c_sq = scenario_N_c^2)]
      X <- model.matrix(form, nd)
      colnames(X)[colnames(X) == "(Intercept)"] <- "intrcpt"
      stopifnot(identical(colnames(X), names(beta)))
      L <- colSums(X * nd$marginal_weight)
      estimate <- sum(L * beta)
      se <- sqrt(as.numeric(t(L) %*% covariance %*% L))
      predictions[[length(predictions) + 1L]] <- data.table(
        region = reg,
        region_abbr = support[region == reg, region_abbr],
        N_scenario = c("Low N", "Medium N", "High N")[scenario],
        N_rate = n_value,
        log1p_EDH30 = heat_value,
        EDH30_total = expm1(heat_value),
        predicted_lnRR = estimate,
        se = se,
        ci_lb_lnRR = estimate - qnorm(.975) * se,
        ci_ub_lnRR = estimate + qnorm(.975) * se,
        effect_size_pct = 100 * expm1(estimate),
        ci_lb_pct = 100 * expm1(estimate - qnorm(.975) * se),
        ci_ub_pct = 100 * expm1(estimate + qnorm(.975) * se)
      )
    }
  }
}
predictions <- rbindlist(predictions)
stopifnot(nrow(predictions) == 3636L, all(is.finite(as.matrix(predictions[, 4:13]))),
          all(predictions$ci_lb_pct <= predictions$effect_size_pct & predictions$effect_size_pct <= predictions$ci_ub_pct))
fwrite(predictions, file.path(out, "regional_N_EDH30_predictions_REML.csv"))

# 2. 沿用原版三行四列、区域顺序、Q25/Q50/Q75 配色和 BH 校正标注。
plot_order <- c("NAMR", "EUR", "NWC", "NEC", "NAMI", "WAS", "SWC", "NCP", "SAM", "WAF", "EAF", "SAS")
scenario_cols <- c("Low N" = "#6C9BD2", "Medium N" = "#FDD35C", "High N" = "#3071B9")
predictions[, region_abbr := factor(region_abbr, levels = plot_order)]
predictions[, N_scenario := factor(N_scenario, levels = names(scenario_cols))]
labels <- merge(support, regional_tests[, .(region_abbr, p_BH)], by = "region_abbr")
labels[, region_abbr := factor(region_abbr, levels = plot_order)]
labels[, panel_label := paste0(
  region_abbr, "\nK=", format(effects, big.mark = ",", trim = TRUE), ", n=", environments,
  "\nN=", round(N_q25), "/", round(N_q50), "/", round(N_q75),
  "\n", fifelse(p_BH < .001, "P<0.001", paste0("P=", format.pval(p_BH, digits = 2, eps = .001)))
)]
p <- ggplot(predictions, aes(log1p_EDH30, effect_size_pct, colour = N_scenario, fill = N_scenario)) +
  geom_ribbon(aes(ymin = ci_lb_pct, ymax = ci_ub_pct), alpha = .18, colour = NA) +
  geom_line(linewidth = .55) +
  geom_text(data = labels, aes(x = Inf, y = Inf, label = panel_label), inherit.aes = FALSE,
            hjust = 1.08, vjust = 1.08, family = "Arial", size = 2.1, lineheight = .92, colour = "black") +
  facet_wrap(~ region_abbr, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = scenario_cols, breaks = names(scenario_cols), labels = c("Q25", "Q50", "Q75"), name = NULL) +
  scale_fill_manual(values = scenario_cols, guide = "none") +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = function(limits) seq(limits[1], limits[2], length.out = 4),
                     labels = scales::label_number(accuracy = 1), expand = expansion(mult = c(.06, .08))) +
  labs(x = "Within-region ln(1+EDH30)", y = "Model-predicted effect size (%)") +
  theme_classic(base_size = 6, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.text = element_text(colour = "black", size = 5.5),
    axis.title = element_text(colour = "black", size = 6.5),
    axis.ticks = element_line(linewidth = .4, colour = "black"),
    axis.ticks.length = grid::unit(1, "mm"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = .4),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "top",
    legend.key.width = grid::unit(5, "mm"),
    panel.spacing.x = grid::unit(3.5, "mm"),
    panel.spacing.y = grid::unit(1.8, "mm")
  )
ggsave(file.path(out, "Figure3b.pdf"), p, width = 235, height = 102, units = "mm", device = cairo_pdf)
# 预览仅在用户指定外部路径时生成，不增加投稿目录文件。
preview <- Sys.getenv("FIGURE3B_PREVIEW", "")
if (nzchar(preview)) ggsave(preview, p, width = 235, height = 102, units = "mm", dpi = 250)
message("Completed Figure3b: 3636 predictions; 12 regional panels.")

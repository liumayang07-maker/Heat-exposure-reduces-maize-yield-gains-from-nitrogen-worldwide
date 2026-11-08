invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(metafor)
  library(patchwork)
})

# 1. 读取本目录的新拟合结果，沿用原版组内百分位预测。
entry <- grep("^--file=", commandArgs(), value = TRUE)
module <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
table_dir <- module
bundle <- readRDS(file.path(module, "model_REML.rds"))
fit <- bundle$fit
model_formula <- bundle$formula
d <- fread(file.path(module, "model_input.csv"))
joint <- fread(file.path(module, "interaction_joint_Wald.csv"))

groups <- c("0_100", "100_200", "200_300", "300_400", "400_600")
group_labels <- c("0-100", "100-200", "200-300", "300-400", "400-600")
d[, N_group := factor(N_group, levels = groups)]
d[, country := relevel(factor(country), "China")]
d[, irrigation := factor(irrigation, levels = c("Rainfed", "Irrigation"))]
d[, N_strategy := factor(N_strategy, levels = c("S", "SV", "SVR", "SR", "VR", "V", "R"))]
stopifnot(!anyNA(d$plant_density_c), !anyNA(d$N_strategy))
for (variable in c("reference_ID", "ref_site_year_ID", "effect_id")) {
  set(d, j = variable, value = factor(d[[variable]]))
}

beta <- as.numeric(coef(fit))
names(beta) <- names(coef(fit))
coef_vcov <- vcov(fit)
percentile_grid <- 0:100

predict_group_curve <- function(group_name) {
  group_data <- copy(d[N_group == group_name])
  if (!nrow(group_data)) stop("No observations for N group: ", group_name)

  # 热暴露百分位基于组内唯一环境，避免重复处理行改变分位数。
  environment_heat <- unique(group_data[, .(env_ID, EDH30_total, heat_c)])
  if (environment_heat[, anyDuplicated(env_ID)]) {
    stop("EDH30 is not unique within env_ID for N group: ", group_name)
  }
  heat_quantiles <- as.numeric(quantile(
    environment_heat$heat_c, probs = percentile_grid / 100,
    names = FALSE, type = 7
  ))
  raw_log_quantiles <- as.numeric(quantile(
    log1p(environment_heat$EDH30_total), probs = percentile_grid / 100,
    names = FALSE, type = 7
  ))

  # 对同一环境内处理行均分权重，使各环境对边际预测的总权重相同。
  group_data[, environment_rows := .N, by = env_ID]
  group_data[, marginal_weight := 1 / environment_rows]
  group_data[, marginal_weight := marginal_weight / sum(marginal_weight)]

  rbindlist(lapply(seq_along(percentile_grid), function(index) {
    prediction_data <- copy(group_data)
    prediction_data[, heat_c := heat_quantiles[index]]
    X <- model.matrix(model_formula, prediction_data)
    if (colnames(X)[1] == "(Intercept)" && names(beta)[1] == "intrcpt") {
      colnames(X)[1] <- "intrcpt"
    }
    stopifnot(identical(colnames(X), names(beta)))
    L <- colSums(X * prediction_data$marginal_weight)
    mu <- sum(L * beta)
    se <- sqrt(as.numeric(t(L) %*% coef_vcov %*% L))
    data.table(
      N_group = group_name,
      percentile = percentile_grid[index],
      log1p_EDH30 = raw_log_quantiles[index],
      EDH30_total = expm1(raw_log_quantiles[index]),
      predicted_lnRR = mu,
      se = se,
      ci_lb_lnRR = mu - qnorm(0.975) * se,
      ci_ub_lnRR = mu + qnorm(0.975) * se,
      effect_size_pct = 100 * expm1(mu),
      ci_lb_pct = 100 * expm1(mu - qnorm(0.975) * se),
      ci_ub_pct = 100 * expm1(mu + qnorm(0.975) * se),
      effects = nrow(group_data),
      environments = uniqueN(group_data$env_ID)
    )
  }))
}

predictions <- rbindlist(lapply(groups, predict_group_curve))
predictions[, N_group := factor(N_group, levels = groups, labels = group_labels)]
fwrite(predictions, file.path(table_dir, "N_group_EDH30_percentile_curves_REML.csv"))

coverage <- predictions[, .(
  minimum = min(percentile), maximum = max(percentile), points = .N
), by = N_group]
validation <- data.table(
  check = c("five_groups_present", "percentiles_0_to_100", "finite_predictions", "ordered_intervals", "effect_size_transform"),
  passed = c(
    uniqueN(predictions$N_group) == 5L,
    all(coverage$minimum == 0 & coverage$maximum == 100 & coverage$points == 101),
    all(is.finite(predictions$effect_size_pct)),
    all(predictions$ci_lb_pct <= predictions$effect_size_pct & predictions$effect_size_pct <= predictions$ci_ub_pct),
    max(abs(predictions$effect_size_pct - 100 * expm1(predictions$predicted_lnRR))) < 1e-12
  )
)
stopifnot(all(validation$passed))

palette <- c(
  "0-100" = "#0072B2", "100-200" = "#56B4E9", "200-300" = "#009E73",
  "300-400" = "#E69F00", "400-600" = "#D55E00"
)

exposure_distribution <- unique(d[, .(N_group, env_ID, EDH30_total)])
exposure_distribution[, N_group := factor(N_group, levels = groups, labels = group_labels)]
exposure_distribution[, log1p_EDH30 := log1p(EDH30_total)]
plot_density <- ggplot(
  exposure_distribution,
  aes(x = log1p_EDH30, colour = N_group, fill = N_group)
) +
  geom_density(alpha = 0.08, linewidth = 0.45, adjust = 0.9, trim = TRUE) +
  scale_colour_manual(values = palette, guide = "none") +
  scale_fill_manual(values = palette, guide = "none") +
  scale_x_continuous(position = "top", breaks = c(0, 2, 4, 6, 8), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(
    breaks = c(0, 0.15, 0.30),
    limits = c(0, 0.32),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Observed log1p(EDH30)", y = "Density") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.text = element_text(colour = "black", size = 6),
    axis.title = element_text(colour = "black", size = 7),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = grid::unit(1, "mm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    plot.margin = margin(3, 3, 0, 3)
  )

plot_curves <- ggplot(
  predictions,
  aes(x = percentile, y = effect_size_pct, colour = N_group, fill = N_group)
) +
  geom_ribbon(aes(ymin = ci_lb_pct, ymax = ci_ub_pct), alpha = 0.30, colour = NA) +
  geom_line(linewidth = 0.50) +
  annotate(
    "text", x = 30, y = max(predictions$ci_ub_pct) * 0.99,
    label = paste0(
      "Heat × N-group\nQM = ", sprintf("%.2f", joint$QM),
      ", P = ", format.pval(joint$p, digits = 3, eps = 1e-05),
      "\nK = ", format(nrow(d), big.mark = ",")
    ),
    hjust = 0, vjust = 1, size = 2.2, family = "Arial", colour = "black"
  ) +
  scale_colour_manual(values = palette, name = "N rate (kg N/ha)") +
  scale_fill_manual(values = palette, name = "N rate (kg N/ha)") +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100), expand = expansion(mult = 0)) +
  labs(
    x = "log1p(EDH30) percentile within N-rate group",
    y = "Model-predicted effect size (%)"
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black", size = 7),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = grid::unit(1, "mm"),
    legend.position = "inside",
    legend.position.inside = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(5, "mm"),
    legend.title = element_text(size = 6.5, colour = "black"),
    legend.text = element_text(size = 6, colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  )

plot_combined <- (plot_density / plot_curves) + plot_layout(heights = c(1, 4))

# 2. 保存实际密度曲线坐标与最终面板；不导出额外探索图。
density_layer <- as.data.table(ggplot_build(plot_density)$data[[1]])
density_points <- density_layer[, .(
  N_group = group_labels[group],
  log1p_EDH30 = x,
  density = density
)]
stopifnot(all(is.finite(density_points$log1p_EDH30)), all(is.finite(density_points$density)))
fwrite(density_points, file.path(module, "exposure_density_curves.csv"))
ggsave(
  file.path(module, "Figure2b.pdf"),
  plot_combined,
  width = 80,
  height = 90,
  units = "mm",
  device = cairo_pdf
)
message(format(Sys.time(), "%F %T"), " | Figure2b plot completed")

#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# 1. 从当前模型重新计算共同氮量范围内的环境等权边际预测。
entry <- grep("^--file=", commandArgs(), value = TRUE)
module <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
bundle <- readRDS(file.path(module, "model_REML.rds"))
fit <- bundle$fit
model_formula <- bundle$formula
d <- fread(file.path(module, "model_input.csv"))
parameters <- fread(file.path(module, "model_parameters.csv"))
transformations <- fread(file.path(module, "transformation_parameters.csv"))
n_center <- transformations[target == "N_c", center]
n_divisor <- transformations[target == "N_c", divisor]
heat_levels <- c(
  "Low heat exposure", "Early heat exposure",
  "Middle heat exposure", "Late heat exposure"
)
d[, Env_Type := factor(Env_Type, levels = heat_levels)]
d[, country := relevel(factor(country), "China")]
d[, irrigation := factor(irrigation, levels = c("Rainfed", "Irrigation"))]
d[, N_strategy := factor(N_strategy, levels = c("S", "SV", "SVR", "SR", "VR", "V", "R"))]
d[, standard_weight := 1 / .N, by = env_ID]
w <- d$standard_weight / sum(d$standard_weight)
mean_design <- function(heat_level, n_rate) {
  nd <- copy(d)
  nd[, Env_Type := factor(heat_level, levels = heat_levels)]
  nd[, N_c := (n_rate - n_center) / n_divisor]
  nd[, N_c_sq := N_c^2]
  X <- model.matrix(model_formula, nd)
  colnames(X)[colnames(X) == "(Intercept)"] <- "intrcpt"
  stopifnot(identical(colnames(X), names(coef(fit))))
  as.numeric(crossprod(w, X))
}
n_grid <- seq(parameters$common_N_min, parameters$common_N_max, length.out = 251)
prediction_grid <- CJ(Env_Type = heat_levels, N_rate = n_grid, unique = TRUE)
M <- do.call(rbind, Map(mean_design, prediction_grid$Env_Type, prediction_grid$N_rate))
estimate <- drop(M %*% as.numeric(coef(fit)))
se <- sqrt(pmax(0, rowSums((M %*% vcov(fit)) * M)))
z <- cbind(prediction_grid, data.table(
  estimate = estimate,
  se = se,
  ci_lb = estimate - qnorm(0.975) * se,
  ci_ub = estimate + qnorm(0.975) * se
))
z[, effect_size_pct := 100 * expm1(estimate)]
z[, effect_size_ci_lb := 100 * expm1(ci_lb)]
z[, effect_size_ci_ub := 100 * expm1(ci_ub)]
stopifnot(
  nrow(z) == 1004L,
  all(is.finite(z$effect_size_pct)),
  all(z$effect_size_ci_lb <= z$effect_size_pct),
  all(z$effect_size_ci_ub >= z$effect_size_pct)
)
fwrite(z, file.path(module, "adjusted_N_response_curves_REML.csv"))
wald <- fread(file.path(module, "interaction_Wald_tests.csv"))[
  test == "All environment by N interactions"
]
k <- parameters$effects
stopifnot(nrow(wald) == 1L)
p_text <- if (wald$p < 1e-05) {
  "P < 1e-05"
} else {
  paste0("P = ", format.pval(wald$p, digits = 3))
}

cols <- c(
  "Low heat exposure" = "#66A5ED",
  "Early heat exposure" = "#F6C45C",
  "Middle heat exposure" = "#DB3E8C",
  "Late heat exposure" = "#7758A3"
)
labels <- c(
  "Low heat exposure" = "Low",
  "Early heat exposure" = "Early",
  "Middle heat exposure" = "Middle",
  "Late heat exposure" = "Late"
)

# 2. 密度图保留原版的比较行计数口径，仅显示共同范围上限以内的数据。
n_input <- copy(d[, .(N_rate, Env_Type)])
n_input[, N_rate := suppressWarnings(as.numeric(N_rate))]
n_input <- n_input[is.finite(N_rate) & N_rate >= 0 & Env_Type %in% names(labels)]
stopifnot(nrow(n_input) > 0L)

z[, Env_Type := factor(Env_Type, levels = names(labels))]

curve_max <- max(z$N_rate)
n_input <- n_input[N_rate <= curve_max]
stopifnot(nrow(n_input) > 0L)
n_input[, Env_Type := factor(Env_Type, levels = names(labels))]

p <- ggplot(z, aes(N_rate, effect_size_pct, colour = Env_Type, fill = Env_Type)) +
  geom_ribbon(aes(ymin = effect_size_ci_lb, ymax = effect_size_ci_ub),
              alpha = 0.30, colour = NA) +
  geom_line(linewidth = 0.50) +
  annotate(
    "text",
    x = 200,
    y = min(z$effect_size_ci_lb) + 0.20 * diff(range(c(z$effect_size_ci_lb, z$effect_size_ci_ub))),
    label = paste0(
      "Environment x N rate\nQM = ", sprintf("%.2f", wald$QM),
      ", ", p_text,
      "\nK = ", format(k, big.mark = ",")
    ),
    hjust = 0, vjust = 1, size = 2.2, family = "Arial", colour = "black"
  ) +
  scale_colour_manual(values = cols, labels = labels, name = "Environmental type") +
  scale_fill_manual(values = cols, labels = labels, guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 400, 100), limits = c(0, curve_max),
    expand = expansion(mult = 0)
  ) +
  labs(x = "N rate (kg N ha-1)", y = "Model-estimated effect size (%)") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black", size = 7),
    legend.position = "inside",
    legend.position.inside = c(0.03, 0.97),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.80), colour = NA),
    legend.margin = margin(2, 3, 2, 3, unit = "pt"),
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(5, "mm"),
    legend.text = element_text(size = 6, colour = "black"),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = grid::unit(1, "mm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank()
  )

p_density <- ggplot(n_input, aes(N_rate, colour = Env_Type, fill = Env_Type)) +
  geom_density(alpha = 0.08, linewidth = 0.45, adjust = 0.9,
               trim = TRUE) +
  scale_x_continuous(
    breaks = seq(0, 400, 100), limits = c(0, curve_max),
    position = "top", expand = expansion(mult = 0)
  ) +
  labs(x = "N rate (kg N ha-1)", y = "Density") +
  scale_colour_manual(values = cols, labels = labels, guide = "none") +
  scale_fill_manual(values = cols, labels = labels, guide = "none") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = grid::unit(1, "mm"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    plot.margin = margin(3, 3, 0, 3)
  )

combined <- p_density / p + plot_layout(heights = c(1, 4))

# 3. 保存实际密度坐标和最终矢量面板，不保留其他探索性图表。
density_layer <- as.data.table(ggplot_build(p_density)$data[[1]])
density_points <- density_layer[, .(
  Env_Type = names(labels)[group],
  N_rate = x,
  density = density
)]
stopifnot(all(is.finite(density_points$N_rate)), all(is.finite(density_points$density)))
fwrite(density_points, file.path(module, "N_rate_density_curves.csv"))
ggsave(
  file.path(module, "Figure2c.pdf"),
  combined,
  width = 80,
  height = 90,
  units = "mm",
  device = cairo_pdf
)
message(format(Sys.time(), "%F %T"), " | Figure2c plot completed")

#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
  library(ggplot2)
  library(patchwork)
})

# 1. 读取当前模型及输入；沿用原版比较行等权的协变量标准化。
entry <- grep("^--file=", commandArgs(), value = TRUE)
module <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
tab <- function(x, name) fwrite(x, file.path(module, paste0(name, ".csv")))
bundle <- readRDS(file.path(module, "model_REML.rds"))
fit <- bundle$fit
form <- bundle$formula
d <- fread(file.path(module, "model_input.csv"))
patterns <- c("Early heat exposure", "Middle heat exposure", "Late heat exposure")
d[, Env_Type := factor(Env_Type, levels = patterns)]
d[, country := factor(country)]
d[, irrigation := factor(irrigation, levels = c("Rainfed", "Irrigation"))]
d[, N_strategy := factor(N_strategy, levels = c("S", "SV", "SVR", "SR", "VR", "V", "R"))]
transformations <- fread(file.path(module, "transformation_parameters.csv"))
heat_center <- transformations[target == "heat_c", center]
n_center <- transformations[target == "N_c", center]
coef_names <- names(coef(fit))
wald_from_L <- function(L, label) {
  L <- as.matrix(L)
  est <- as.numeric(L %*% coef(fit))
  S <- L %*% vcov(fit) %*% t(L)
  eig <- eigen((S + t(S)) / 2, symmetric = TRUE)
  tol <- max(dim(S)) * max(abs(eig$values)) * .Machine$double.eps
  keep <- eig$values > tol
  if (!any(keep)) stop("Wald contrast covariance has zero rank: ", label)
  Sinv <- eig$vectors[, keep, drop = FALSE] %*%
    diag(1 / eig$values[keep], nrow = sum(keep)) %*%
    t(eig$vectors[, keep, drop = FALSE])
  QM <- as.numeric(crossprod(est, Sinv %*% est))
  df <- sum(keep)
  data.table(test = label, QM = QM, df = df,
             p = pchisq(QM, df = df, lower.tail = FALSE))
}

# 2. 各热模式使用唯一环境的百分位网格；其他协变量按比较行分布标准化。
env_support <- unique(d[, .(env_ID, Env_Type, EDH30_total)])
stopifnot(!anyDuplicated(env_support$env_ID))
env_support[, log1p_EDH30 := log1p(EDH30_total)]
probs <- seq(0, 1, by = 0.01)
grid <- rbindlist(lapply(patterns, function(pattern) {
  x <- env_support[Env_Type == pattern, log1p_EDH30]
  data.table(
    Env_Type = pattern,
    percentile = probs * 100,
    log1p_EDH30 = as.numeric(quantile(x, probs = probs, names = FALSE))
  )
}))
grid[, Env_Type := factor(Env_Type, levels = patterns)]

X_obs <- model.matrix(form, data = d)
colnames(X_obs)[colnames(X_obs) == "(Intercept)"] <- "intrcpt"
base <- as.numeric(crossprod(rep(1 / nrow(X_obs), nrow(X_obs)), X_obs))
names(base) <- colnames(X_obs)

design_row <- function(pattern, heat_value, n_rate) {
  nd <- copy(d[1])
  nd[, Env_Type := factor(pattern, levels = patterns)]
  nd[, heat_c := heat_value]
  nd[, N_c := (n_rate - n_center) / 100]
  nd[, N_c_sq := N_c^2]
  x <- model.matrix(form, data = nd)
  colnames(x)[colnames(x) == "(Intercept)"] <- "intrcpt"
  x[, coef_names, drop = FALSE]
}

prediction_row <- function(pattern, heat_value, n_rate) {
  x <- design_row(pattern, heat_value - heat_center, n_rate)
  focus <- grepl("^Env_Type|^heat_c$|^N_c$|^N_c_sq$|:", names(base))
  out <- base[coef_names]
  out[focus] <- x[1, names(out)[focus]]
  matrix(out, nrow = 1)
}

fixed_N <- c(100, 200, 300)
predict_grid <- function(n_rate) {
  out <- copy(grid)
  M <- do.call(rbind, Map(
    prediction_row, as.character(out$Env_Type), out$log1p_EDH30, n_rate
  ))
  b <- as.numeric(coef(fit)); C <- vcov(fit)
  out[, estimate := as.numeric(M %*% b)]
  out[, se := sqrt(pmax(0, rowSums((M %*% C) * M)))]
  out[, `:=`(
    ci_lb = estimate - qnorm(0.975) * se,
    ci_ub = estimate + qnorm(0.975) * se,
    effect_size_pct = 100 * expm1(estimate),
    effect_size_ci_lb = 100 * expm1(estimate - qnorm(0.975) * se),
    effect_size_ci_ub = 100 * expm1(estimate + qnorm(0.975) * se),
    N_control_kg_ha = n_rate
  )]
  out
}

fixed_predictions <- rbindlist(lapply(fixed_N, predict_grid))
tab(fixed_predictions, "fixed_N_percentile_curves_REML")

# 3. 固定施氮量下，联合检验三种模式曲线的水平及热暴露斜率是否相同。
curve_tests_at_N <- function(n_rate) {
  x0 <- lapply(patterns, design_row, heat_value = 0, n_rate = n_rate)
  x1 <- lapply(patterns, design_row, heat_value = 1, n_rate = n_rate)
  names(x0) <- names(x1) <- patterns
  slope <- Map(function(a, b) b - a, x0, x1)

  L_overall <- rbind(
    middle_level = x0[["Middle heat exposure"]] - x0[["Early heat exposure"]],
    middle_slope = slope[["Middle heat exposure"]] - slope[["Early heat exposure"]],
    late_level = x0[["Late heat exposure"]] - x0[["Early heat exposure"]],
    late_slope = slope[["Late heat exposure"]] - slope[["Early heat exposure"]]
  )
  overall <- wald_from_L(L_overall, "Overall equality of the three heat-pattern curves")
  overall[, `:=`(
    N_control_kg_ha = n_rate,
    comparison = "Early = Middle = Late",
    scope = "overall"
  )]
  overall
}

fixed_N_tests <- rbindlist(lapply(fixed_N, curve_tests_at_N), use.names = TRUE, fill = TRUE)
tab(fixed_N_tests, "fixed_N_curve_Wald_tests_REML")

# 4. 仅输出投稿面板，不导出三个子图或其他分析图。
palette <- c(
  "Early heat exposure" = "#EEBF5D",
  "Middle heat exposure" = "#CC3F84",
  "Late heat exposure" = "#71559B"
)
common_ylim <- c(0, 150)
common_y_breaks <- c(0, 75, 150)

make_plot <- function(n_rate, show_x = TRUE, show_legend = TRUE) {
  pred <- fixed_predictions[N_control_kg_ha == n_rate]
  test_row <- fixed_N_tests[N_control_kg_ha == n_rate & scope == "overall"]
  p_label <- if (test_row$p < 0.001) "P < 0.001" else sprintf("P = %.3f", test_row$p)
  annotation <- sprintf(
    "QM = %.2f, %s, K = %s",
    test_row$QM, p_label, format(nrow(d), big.mark = ",")
  )
  n_label <- paste0("N = ", n_rate, " kg N ha-1")
  ggplot(pred, aes(percentile, effect_size_pct, colour = Env_Type, fill = Env_Type)) +
    geom_ribbon(aes(ymin = effect_size_ci_lb, ymax = effect_size_ci_ub),
                alpha = 0.14, colour = NA) +
    geom_line(linewidth = 0.65) +
    annotate("text", x = 2, y = 3, label = annotation,
             hjust = 0, vjust = 0, size = 2.1, family = "Arial") +
    annotate("text", x = 98, y = 147, label = n_label,
             hjust = 1, vjust = 1, size = 2.1, family = "Arial") +
    scale_colour_manual(values = palette, labels = c("Early", "Middle", "Late"), name = NULL) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(breaks = common_y_breaks, expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = common_ylim) +
    labs(
      x = "log1p(EDH30) percentile within heat environment",
      y = "Model-estimated effect size (%)"
    ) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      legend.position = if (show_legend) "inside" else "none",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      legend.background = element_blank(),
      legend.key.width = grid::unit(4, "mm"),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.text = element_text(size = 5.5, family = "Arial"),
      legend.spacing.x = grid::unit(1, "mm"),
      axis.text = element_text(size = 6, family = "Arial", colour = "black"),
      axis.title = element_text(size = 7, family = "Arial"),
      axis.text.x = if (show_x) element_text(size = 6, family = "Arial", colour = "black") else element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.5) else element_blank(),
      axis.title.x = if (show_x) element_text(size = 7, family = "Arial") else element_blank(),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = grid::unit(1, "mm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(1.2, 1.5, 1.2, 1.5, unit = "mm")
    )
}

combined_plot <-
  make_plot(100, show_x = FALSE, show_legend = TRUE) /
  make_plot(200, show_x = FALSE, show_legend = FALSE) /
  make_plot(300, show_x = TRUE, show_legend = FALSE) +
  plot_layout(ncol = 1, heights = c(1, 1, 1), axis_titles = "collect_y")

stopifnot(
  nrow(fixed_predictions) == 909L,
  all(is.finite(fixed_predictions$effect_size_pct)),
  all(fixed_predictions$effect_size_ci_lb <= fixed_predictions$effect_size_pct),
  all(fixed_predictions$effect_size_ci_ub >= fixed_predictions$effect_size_pct),
  nrow(fixed_N_tests) == 3L
)
ggsave(
  file.path(module, "Figure2d.pdf"),
  combined_plot,
  width = 80,
  height = 90,
  units = "mm",
  device = cairo_pdf
)
message(format(Sys.time(), "%F %T"), " | Figure2d plot completed")

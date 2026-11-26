#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(metafor)
  library(ggplot2)
  library(patchwork)
})

if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}

# 1. 读取 Figure 2a 投稿版输入和 SD 来源，按 data_ID 对齐。
entry <- grep("^--file=", commandArgs(), value = TRUE)
module_dir <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
data_path <- file.path(module_dir, "..", "Figure2a", "comparison_input.csv")
output_path <- file.path(module_dir, "publication_bias_assessment_abc.pdf")
figure_width_mm <- 288
figure_height_mm <- 80
input_columns <- c(
  "data_ID", "reference_ID", "site_year_ID", "Label_of_control_treatment",
  "mean_c", "sd_c", "n_c", "mean_t", "sd_t", "n_t", "lnRR"
)
message("Reading submission data")
dt <- fread(
  data_path,
  encoding = "UTF-8",
  na.strings = c("", "NA"),
  select = input_columns
)
stopifnot(
  nrow(dt) == 4613L,
  !anyNA(dt),
  !anyDuplicated(dt$data_ID),
  uniqueN(dt$reference_ID) == 281L,
  all(dt$Label_of_control_treatment == "N0"),
  all(is.finite(dt$lnRR)),
  all(dt$mean_c > 0 & dt$mean_t > 0),
  all(dt$sd_c >= 0 & dt$sd_t >= 0),
  all(dt$n_c > 0 & dt$n_t > 0)
)
sd_input <- fread(
  file.path(module_dir, "..", "Supplementary_Table_SD_imputation", "current_SD_input.csv"),
  select = c("data_ID", "sd_c_lab", "sd_t_lab", "sd_c", "sd_t")
)
stopifnot(!anyDuplicated(sd_input$data_ID), setequal(sd_input$data_ID, dt$data_ID))
sd_input <- sd_input[match(dt$data_ID, data_ID)]
# origin 包括原文报告或可换算的 SD。核对保留值与投稿版一致。
keep <- sd_input$sd_c_lab %in% "origin" & sd_input$sd_t_lab %in% "origin"
stopifnot(
  sum(keep) == 1326L,
  all(dt$sd_c[keep] == sd_input$sd_c[keep]),
  all(dt$sd_t[keep] == sd_input$sd_t[keep])
)

# 2. 分别构建全数据和未插补 SD 子集的共享对照采样协方差矩阵。
prepare_comparisons <- function(dt) {
  dt[, row_id := .I]
  dt[, treatment_var_lnRR := sd_t^2 / (n_t * mean_t^2)]
  dt[, control_var_lnRR := sd_c^2 / (n_c * mean_c^2)]
  dt[, vi_lnRR := treatment_var_lnRR + control_var_lnRR]
  dt[, sei_lnRR := sqrt(vi_lnRR)]
  stopifnot(
    all(is.finite(dt$sei_lnRR) & dt$sei_lnRR > 0),
    max(abs(dt$lnRR - log(dt$mean_t / dt$mean_c))) < 1e-6
  )

  # 2. 按 Figure 2a 规则识别共享对照；同块比较的协方差等于对照方差项。
  message("Building shared-control sampling covariance matrix")
  dt[, ref_site_year_ID := paste(reference_ID, site_year_ID, sep = "__")]
  dt[, control_id := paste(
    reference_ID,
    site_year_ID,
    Label_of_control_treatment,
    sprintf("%.10g", mean_c),
    sprintf("%.10g", sd_c),
    sprintf("%.10g", n_c),
    sep = "||"
  )]
  setorder(dt, control_id, row_id)
  dt[, model_row := .I]
  dt[, effect_id := factor(model_row)]
  dt[, reference_ID := factor(as.character(reference_ID))]
  dt[, ref_site_year_ID := factor(ref_site_year_ID)]

  block_summary <- dt[, .(treatment_arms = .N), by = control_id]
  shared_ids <- block_summary[treatment_arms > 1, control_id]
  shared_pairs <- dt[control_id %in% shared_ids, {
    pairs <- utils::combn(model_row, 2)
    .(
      i = c(pairs[1, ], pairs[2, ]),
      j = c(pairs[2, ], pairs[1, ]),
      x = rep(control_var_lnRR[1], 2 * ncol(pairs))
    )
  }, by = control_id]
  V_sparse <- sparseMatrix(
    i = c(dt$model_row, shared_pairs$i),
    j = c(dt$model_row, shared_pairs$j),
    x = c(dt$vi_lnRR, shared_pairs$x),
    dims = c(nrow(dt), nrow(dt))
  )
  stopifnot(isSymmetric(V_sparse), all(diag(V_sparse) == dt$vi_lnRR))
  list(data = dt, covariance = V_sparse)
}

all_comparisons <- prepare_comparisons(copy(dt))
observed_comparisons <- prepare_comparisons(dt[keep])

# 3. 用总体 REML 模型确定中心线，用多层 Egger 型回归评估不对称。
make_egger_plot <- function(comparisons, panel_label) {
  dt <- comparisons$data
  V_sparse <- comparisons$covariance
  random_effects <- list(
    ~ 1 | reference_ID,
    ~ 1 | ref_site_year_ID,
    ~ 1 | effect_id
  )

  message("Fitting intercept-only REML model")
  fit_overall <- rma.mv(
    yi = lnRR,
    V = V_sparse,
    mods = ~ 1,
    random = random_effects,
    data = dt,
    method = "REML",
    sparse = TRUE,
    control = list(optimizer = "nlminb")
  )
  overall_lnRR <- unname(coef(fit_overall)[1])

  message("Fitting multilevel Egger-type REML model")
  fit_egger <- rma.mv(
    yi = lnRR,
    V = V_sparse,
    mods = ~ sei_lnRR,
    random = random_effects,
    data = dt,
    method = "REML",
    sparse = TRUE,
    control = list(optimizer = "nlminb")
  )
  egger_coefficients <- as.matrix(coef(summary(fit_egger)))
  se_coefficient <- egger_coefficients["sei_lnRR", ]
  stopifnot(
    fit_overall$k == nrow(dt),
    fit_egger$k == nrow(dt),
    all(is.finite(se_coefficient))
  )

  plot_data <- dt[, .(data_ID, lnRR, sei_lnRR)]
  funnel_bounds <- data.table(
    sei_lnRR = seq(0, max(plot_data$sei_lnRR), length.out = 200)
  )
  funnel_bounds[, upper := overall_lnRR + stats::qnorm(0.975) * sei_lnRR]
  funnel_bounds[, lower := overall_lnRR - stats::qnorm(0.975) * sei_lnRR]

  egger_p <- unname(se_coefficient["pval"])
  egger_label <- if (egger_p < 0.00001) {
    'italic(P[Egger])*" < 0.00001"'
  } else {
    sprintf('italic(P[Egger])*" = %.5f"', egger_p)
  }

  funnel_plot <- ggplot(plot_data, aes(x = lnRR, y = sei_lnRR)) +
    geom_point(
      shape = 21,
      size = 0.85,
      stroke = 0.125,
      color = "grey25",
      fill = "grey65",
      alpha = 0.75
    ) +
    geom_line(
      data = funnel_bounds,
      aes(x = lower, y = sei_lnRR),
      inherit.aes = FALSE,
      color = "grey45",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_line(
      data = funnel_bounds,
      aes(x = upper, y = sei_lnRR),
      inherit.aes = FALSE,
      color = "grey45",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = overall_lnRR,
      color = "grey45",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    scale_y_reverse(expand = expansion(mult = c(0.04, 0.04))) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = egger_label,
      parse = TRUE,
      family = "Arial",
      size = 9 / ggplot2::.pt,
      hjust = 1.04,
      vjust = -1.1
    ) +
    labs(x = "lnRR", y = "Standard error") +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.ticks.length = grid::unit(1.5, "mm")
    )
  message(sprintf(
    "Panel %s: %d comparisons, %d publications; lnRR = %.12f; Egger slope = %.12f; P = %.12f",
    panel_label, nrow(dt), uniqueN(dt$reference_ID), overall_lnRR,
    se_coefficient[["estimate"]], egger_p
  ))
  list(
    plot = funnel_plot,
    points = plot_data,
    overall_model = fit_overall,
    egger_model = fit_egger,
    egger_result = data.table(
      panel = panel_label,
      comparisons = nrow(dt),
      publications = uniqueN(dt$reference_ID),
      overall_lnRR = overall_lnRR,
      slope = se_coefficient[["estimate"]],
      se = se_coefficient[["se"]],
      z = se_coefficient[["zval"]],
      p = egger_p,
      ci_lower = se_coefficient[["ci.lb"]],
      ci_upper = se_coefficient[["ci.ub"]]
    )
  )
}

all_results <- make_egger_plot(all_comparisons, "a")
observed_results <- make_egger_plot(observed_comparisons, "b")

# 4. 用 GLS 合并同一篇文献内的比较，再按文献执行 L0 剪补法。
dt <- all_comparisons$data
V_sparse <- all_comparisons$covariance
overall_model <- fread(file.path(module_dir, "..", "Figure2a", "overall_REML_results.csv"))
stopifnot(
  nrow(overall_model) == 1L,
  overall_model$k == nrow(dt),
  overall_model$unique_references == uniqueN(dt$reference_ID)
)

# 2. 用广义最小二乘计算每篇文献的合并 lnRR 和方差。
# 文献内协方差包含共享对照、站点年和比较层变异。文献间方差由 rma.uni 重新估计。
message("Combining ", nrow(dt), " comparisons within ", uniqueN(dt$reference_ID), " publications")
publication_effects <- dt[, {
  site_year_match <- outer(ref_site_year_ID, ref_site_year_ID, "==")
  within_covariance <- as.matrix(V_sparse[model_row, model_row]) +
    overall_model$sigma2_ref_site_year_ID * site_year_match
  diag(within_covariance) <- diag(within_covariance) + overall_model$sigma2_effect_id
  precision_weights <- solve(within_covariance, rep(1, .N))
  total_precision <- sum(precision_weights)
  .(
    comparisons = .N,
    lnRR = sum(precision_weights * lnRR) / total_precision,
    sampling_variance = 1 / total_precision
  )
}, by = reference_ID]
stopifnot(
  sum(publication_effects$comparisons) == nrow(dt),
  all(is.finite(publication_effects$lnRR)),
  all(publication_effects$sampling_variance > 0)
)

# 3. 按文献拟合随机效应 REML，再用 L0 剪补法自动选择补齐侧。
unadjusted <- rma.uni(
  yi = lnRR,
  vi = sampling_variance,
  slab = as.character(reference_ID),
  data = publication_effects,
  method = "REML"
)
baseline_lnRR <- unname(coef(unadjusted)[1])
stopifnot(abs(baseline_lnRR - overall_model$estimate_lnRR) < 1e-6)

tf <- trimfill(unadjusted, estimator = "L0")
side_fit <- regtest(unadjusted, model = "rma", predictor = "sei", ret.fit = TRUE)$fit
side_slope <- unname(coef(side_fit)[2])
adjusted_lnRR <- unname(coef(tf)[1])
is_filled <- tf$fill
stopifnot(
  sum(is_filled) == tf$k0,
  tf$k == nrow(publication_effects) + tf$k0,
  identical(tf$slab[!is_filled], as.character(publication_effects$reference_ID)),
  isTRUE(all.equal(unname(tf$yi[!is_filled]), publication_effects$lnRR)),
  isTRUE(all.equal(unname(tf$vi[!is_filled]), publication_effects$sampling_variance)),
  tf$side == if (side_slope < 0) "right" else "left"
)

point_labels <- c("Observed publications", "Imputed studies")
plot_data <- data.table(
  study_label = tf$slab,
  is_filled = is_filled,
  lnRR = tf$yi,
  sei_lnRR = sqrt(tf$vi)
)
message(sprintf(
  "Panel c: %d publications; %d imputed on the %s; lnRR before = %.12f; after = %.12f",
  unadjusted$k, tf$k0, tf$side, baseline_lnRR, adjusted_lnRR
))
# 5. 空心蓝圈表示算法补齐的效应。虚线使用补齐后的合并 lnRR 为中心。
plot_data[, status := factor(ifelse(is_filled, point_labels[2], point_labels[1]), levels = point_labels)]
funnel_se <- c(0, max(plot_data$sei_lnRR))
funnel_bounds <- data.table(
  sei_lnRR = rep(funnel_se, 2),
  lnRR = c(
    adjusted_lnRR - stats::qnorm(0.975) * funnel_se,
    adjusted_lnRR + stats::qnorm(0.975) * funnel_se
  ),
  side = rep(c("left", "right"), each = 2)
)
trimfill_plot <- ggplot(plot_data, aes(x = lnRR, y = sei_lnRR)) +
  geom_line(
    data = funnel_bounds,
    aes(group = side),
    colour = "grey45",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = adjusted_lnRR,
    colour = "grey45",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_point(
    aes(colour = status, fill = status),
    shape = 21,
    size = 1.5,
    stroke = 0.38
  ) +
  scale_colour_manual(values = setNames(c("grey35", "#23628A"), point_labels)) +
  scale_fill_manual(values = setNames(c("grey65", "white"), point_labels)) +
  scale_y_reverse(expand = expansion(mult = c(0.04, 0.04))) +
  labs(x = "lnRR", y = "Standard error", colour = NULL, fill = NULL) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.ticks.length = grid::unit(1.5, "mm"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    legend.margin = margin(0, 0, 0, 0)
  )

# 5. 对齐三幅图的绘图区，并用 Cairo 导出矢量图。
tag_theme <- theme(plot.tag = element_text(family = "Arial", face = "bold", size = 11))
combined_plot <- (
  wrap_plots(all_results$plot, observed_results$plot, trimfill_plot, nrow = 1) +
    plot_annotation(tag_levels = "a")
) & tag_theme

cairo_path <- tempfile(fileext = ".pdf")
save_funnel_pdf <- function(plot, filename) {
  grDevices::cairo_pdf(
    filename = filename,
    width = figure_width_mm / 25.4,
    height = figure_height_mm / 25.4,
    family = "Arial"
  )
  on.exit(grDevices::dev.off())

  # 让 Cairo 分别写入相邻刻度和图例标签，避免 PDF 提取时将它们视为同一个文本对象。
  split_pdf_labels <- function(grob) {
    if (inherits(grob, "text")) {
      grob$gp$alpha <- rep_len(c(1, 254 / 255), length(grob$label))
      if (length(grob$label) == 1L) {
        grob$gp$alpha <- 254 / 255
      }
    }
    if (length(grob$children)) {
      grob$children[] <- lapply(grob$children, split_pdf_labels)
    }
    if (length(grob$grobs)) {
      grob$grobs <- lapply(grob$grobs, split_pdf_labels)
    }
    grob
  }
  layout <- patchworkGrob(plot)
  label_indices <- grep("^(axis-[lb]|ylab-l|guide-box-top)-[0-9]+$", layout$layout$name)
  for (index in label_indices) {
    layout$grobs[[index]] <- split_pdf_labels(layout$grobs[[index]])
  }
  grid::grid.draw(layout)
}
save_funnel_pdf(combined_plot, cairo_path)

# 6. Cairo 将纸张尺寸取整到整数 pt。仅校正 PDF 页框，保留矢量内容。
pdf_page_box_code <- r"(
from pathlib import Path
import sys

import pymupdf

source_path, output_path, width_mm, height_mm = sys.argv[1:]
width = float(width_mm) * 72 / 25.4
height = float(height_mm) * 72 / 25.4
output = Path(output_path)
staged_output = output.with_suffix('.tmp.pdf')

with pymupdf.open(source_path) as pdf:
    assert pdf.page_count == 1
    pdf[0].set_mediabox(pymupdf.Rect(0, 0, width, height))
    pdf.save(staged_output, deflate=True)

with pymupdf.open(staged_output) as pdf:
    page = pdf[0]
    assert abs(page.rect.width - width) < 0.002
    assert abs(page.rect.height - height) < 0.002
    lines = page.get_text().splitlines()
    assert all(label in lines for label in ('a', 'b', 'c'))
    assert lines.count('lnRR') == 3 and lines.count('Standard error') == 3

staged_output.replace(output)
)"
# Windows ???? Python Launcher??????? python3?
use_launcher <- .Platform$OS.type == "windows" && nzchar(Sys.which("py"))
python <- if (use_launcher) "py" else "python3"
python_args <- if (use_launcher) c("-3", "-") else "-"
status <- system2(
  python,
  args = c(
    python_args,
    shQuote(cairo_path),
    shQuote(output_path),
    figure_width_mm,
    figure_height_mm
  ),
  input = strsplit(pdf_page_box_code, "\n", fixed = TRUE)[[1]]
)
unlink(cairo_path)
if (status != 0L) {
  stop("PDF page-box adjustment failed; check Python 3 and PyMuPDF.")
}
message("Saved ", output_path)


# 7. 保存模型对象、摘要、关键结果和三面板绘图数据。
models <- list(
  all_overall = all_results$overall_model,
  all_egger = all_results$egger_model,
  observed_sd_overall = observed_results$overall_model,
  observed_sd_egger = observed_results$egger_model,
  publication_before_fill = unadjusted,
  publication_side_selection = side_fit,
  publication_after_fill = tf
)
saveRDS(models, file.path(module_dir, "publication_bias_models.rds"))

model_summaries <- unlist(lapply(names(models), function(model_name) {
  c(model_name, capture.output(print(summary(models[[model_name]]))), "")
}), use.names = FALSE)
writeLines(
  c(model_summaries, capture.output(sessionInfo())),
  file.path(module_dir, "publication_bias_model_summary.txt")
)

fwrite(
  rbindlist(list(all_results$egger_result, observed_results$egger_result)),
  file.path(module_dir, "egger_results.csv")
)
trimfill_results <- data.table(
  stage = c("Before fill", "After fill"),
  observed_publications = unadjusted$k,
  added_effects = c(0L, tf$k0),
  imputation_side = c(NA_character_, tf$side),
  side_selection_slope = c(NA_real_, side_slope),
  lnRR = c(baseline_lnRR, adjusted_lnRR),
  ci_lower = c(unadjusted$ci.lb, tf$ci.lb),
  ci_upper = c(unadjusted$ci.ub, tf$ci.ub)
)
trimfill_results[, yield_change_percent := 100 * expm1(lnRR)]
fwrite(trimfill_results, file.path(module_dir, "trim_and_fill_results.csv"))

funnel_points <- rbindlist(
  list(a = all_results$points, b = observed_results$points, c = plot_data),
  use.names = TRUE,
  fill = TRUE,
  idcol = "panel"
)
fwrite(funnel_points, file.path(module_dir, "funnel_plot_input.csv"))

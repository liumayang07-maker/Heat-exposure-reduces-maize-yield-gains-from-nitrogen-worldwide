#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# 1. 读取实际比较值及已拟合结果，不在绘图时重新拟合。
entry <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript plot_histogram.R [OUTPUT_DIRECTORY]")
out <- if (length(args)) args[[1]] else root
dir.create(out, recursive = TRUE, showWarnings = FALSE)
lnrr_valid <- fread(file = file.path(root, "comparison_input.csv"))
lnrr_valid[, lnRR_num := as.numeric(lnRR)]
overall_summary <- fread(file = file.path(out, "overall_REML_results.csv"))
stopifnot(nrow(lnrr_valid) == overall_summary$k, all(is.finite(lnrr_valid$lnRR_num)))
est <- overall_summary$estimate_lnRR
# 2. 沿用原直方图的分箱、线型和标注。
hist_plot <- ggplot(lnrr_valid, aes(x = lnRR_num)) +
  geom_histogram(bins = 50, fill = "#67BE91", alpha = 0.7, color = "white", linewidth = 0.25) +
  geom_vline(xintercept = quantile(lnrr_valid$lnRR_num, c(.05, .95), na.rm = TRUE),
             color = "#C62828", linewidth = 0.25, linetype = "dashed") +
  geom_vline(xintercept = est, color = "#D55E00", linewidth = 0.25) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3,
           label = sprintf("Overall REML\nestimate: %.2f%%\n95%% CI:\n%.2f%% to %.2f%%",
                           overall_summary$response_change_percent,
                           overall_summary$ci_lb_response_change_percent,
                           overall_summary$ci_ub_response_change_percent),
           family = "Helvetica", size = 7 / (72.27 / 25.4), lineheight = 1) +
  labs(
    x = "lnRR",
    y = "Frequency"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 8) +
  theme(
    line = element_line(linewidth = 0.25),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = grid::unit(1.5, "mm"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

q <- quantile(lnrr_valid$lnRR_num, c(.05, .95))
peak <- max(ggplot_build(hist_plot)$data[[1]]$count)
hist_plot <- hist_plot + annotate(
  "text",
  x = as.numeric(q),
  y = peak * .94,
  label = c("5th", "95th"),
  hjust = c(-.15, 1.15),
  colour = "#C62828",
  size = 7 / (72.27 / 25.4),
  family = "Helvetica"
)
# 3. 导出原始分图及直接用于直方图的分箱计数。
pdf(
  file.path(out, "Figure2a_histogram.pdf"),
  width = 67 / 25.4,
  height = 46 / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)
print(hist_plot)
invisible(dev.off())
bins <- as.data.table(ggplot_build(hist_plot)$data[[1]])
fwrite(bins[, .(bin_left = xmin, bin_right = xmax, bin_center = x, count)],
       file.path(out, "histogram_bins.csv"))

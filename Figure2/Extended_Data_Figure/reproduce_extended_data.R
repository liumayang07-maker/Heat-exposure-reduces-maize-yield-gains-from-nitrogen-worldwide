#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

# 1. 只读取本目录输入；分类由前置 K=3 分析确定，此处不重新聚类。
entry <- grep("^--file=", commandArgs(), value = TRUE)
out <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
environments <- fread(file.path(out, "environment_input.csv"))
profiles <- fread(file.path(out, "profiles_input.csv"))
r1_input <- fread(file.path(out, "R1_input.csv"))
smooth_half_window <- 2L
heatmap_cap_probability <- 0.99
class_names <- c("Low heat exposure", "Early heat exposure", "Middle heat exposure", "Late heat exposure")
class_colors <- c("#CC99CC", "#0072B2", "#E69F00", "#009E73")
pattern_names <- c("Early pattern", "Middle pattern", "Late pattern")
pattern_colors <- setNames(class_colors[2:4], 1:3)
stopifnot(
  nrow(environments) == 856L,
  uniqueN(environments$env_ID) == 856L,
  !anyNA(environments),
  all(environments$heat_exposure_code %in% 0:3)
)
environments[, heat_exposure_class := class_names[heat_exposure_code + 1L]]
counts <- environments[, .(environment_n = .N), by = .(heat_exposure_code, heat_exposure_class)]
setorder(counts, heat_exposure_code)
counts[, percentage := 100 * environment_n / sum(environment_n)]
stopifnot(identical(counts$environment_n, c(246L, 270L, 201L, 139L)))
assignment <- environments[heat_exposure_code > 0L, .(env_ID, cluster = heat_exposure_code)]
stopifnot(
  nrow(profiles) == 61000L,
  !anyDuplicated(profiles[, .(env_ID, bin)]),
  setequal(profiles$env_ID, assignment$env_ID),
  all(is.finite(profiles$EDH30_per_hour)),
  all(profiles$EDH30_per_hour >= 0)
)
profile_check <- profiles[, .(complete = identical(sort(bin), 1:100)), by = env_ID]
stopifnot(all(profile_check$complete))

# 2. 沿用原版 R1 规则：每个环境先取中位数，再计算环境间 type-7 分位数。
stopifnot(
  all(r1_input$env_ID %in% assignment$env_ID),
  !anyNA(r1_input),
  all(is.finite(r1_input$pre_silking_total_GDD)),
  all(is.finite(r1_input$authoritative_growing_season_GDD10_30_C_d)),
  all(r1_input$authoritative_growing_season_GDD10_30_C_d > 0)
)
r1_valid <- unique(r1_input[, .(
  env_ID,
  R1_date,
  r1_relative_GDD_pct = 100 * pre_silking_total_GDD / authoritative_growing_season_GDD10_30_C_d
)])
stopifnot(all(r1_valid$r1_relative_GDD_pct >= 0 & r1_valid$r1_relative_GDD_pct <= 100))
r1_by_env <- r1_valid[, .(
  r1_relative_GDD_pct = median(r1_relative_GDD_pct),
  n_R1_dates = uniqueN(R1_date)
), by = env_ID]
r1_quantiles <- data.table(
  quantile = c("Q25", "Q50", "Q75"),
  probability = c(0.25, 0.5, 0.75),
  r1_relative_GDD_pct = as.numeric(quantile(
    r1_by_env$r1_relative_GDD_pct,
    probs = c(0.25, 0.5, 0.75),
    type = 7
  )),
  valid_environment_n = nrow(r1_by_env)
)
stopifnot(nrow(r1_by_env) == 236L)
r1_guides <- r1_quantiles[c(1, 3), r1_relative_GDD_pct]

# 3. 保留原始强度绘图；仅在每类内部以 Pearson 距离、average linkage 排序。
wide <- dcast(profiles, env_ID ~ bin, value.var = "EDH30_per_hour")
bin_cols <- as.character(1:100)
setcolorder(wide, c("env_ID", bin_cols))
x_raw <- as.matrix(wide[, ..bin_cols])
rownames(x_raw) <- wide$env_ID
x_raw <- x_raw[assignment$env_ID, , drop = FALSE]
order_within_cluster <- function(ids) {
  if (length(ids) <= 2L) return(ids)
  correlation <- cor(t(x_raw[ids, , drop = FALSE]), method = "pearson")
  stopifnot(all(is.finite(correlation)))
  diag(correlation) <- 1
  distance <- 1 - correlation
  distance[distance < 0] <- 0
  fit <- hclust(as.dist(distance), method = "average")
  fit$labels[fit$order]
}
ordered_ids <- unlist(lapply(1:3, function(group) {
  order_within_cluster(assignment[cluster == group, env_ID])
}), use.names = FALSE)
stopifnot(length(ordered_ids) == 610L, !anyDuplicated(ordered_ids))
order_dt <- data.table(env_ID = ordered_ids, display_order = seq_along(ordered_ids))
order_dt <- merge(order_dt, assignment, by = "env_ID", sort = FALSE)
setorder(order_dt, display_order)
heat <- merge(profiles, assignment, by = "env_ID", sort = FALSE)
heat[, env_factor := factor(env_ID, levels = rev(ordered_ids))]
color_cap <- as.numeric(quantile(heat[EDH30_per_hour > 0, EDH30_per_hour], heatmap_cap_probability))
heat[, value_cap := pmin(EDH30_per_hour, color_cap)]

# 4. 每个环境内居中五 bin 平滑，边缘使用可用的三或四 bin；阴影为 Q25-Q75。
smooth_centered_5bin <- function(values) {
  vapply(seq_along(values), function(index) {
    bounds <- pmax(1L, index - smooth_half_window):pmin(length(values), index + smooth_half_window)
    mean(values[bounds])
  }, numeric(1))
}
setorder(heat, env_ID, bin)
heat[, value_5bin := smooth_centered_5bin(EDH30_per_hour), by = env_ID]
curves <- heat[, .(
  mean_value = mean(value_5bin),
  q25_value = as.numeric(quantile(value_5bin, 0.25, type = 7)),
  q75_value = as.numeric(quantile(value_5bin, 0.75, type = 7)),
  environment_n = uniqueN(env_ID)
), by = .(cluster, bin)]
setorder(curves, cluster, bin)

# 5. 输出三个独立矢量面板，不复制整版 TIFF 或生成重复 PNG。
pdf(file.path(out, "Extended_Data_a.pdf"),
    width = 76 / 25.4, height = 34 / 25.4,
    pointsize = 7, family = "Helvetica", useDingbats = FALSE)
par(fig = c(0, 34 / 76, 0, 1), mar = rep(0, 4), family = "sans", xaxs = "i", yaxs = "i")
fractions <- counts$environment_n / sum(counts$environment_n)
pie(fractions,
    labels = rep("", 4), col = class_colors, border = "white", lwd = 0.7,
    clockwise = FALSE, init.angle = 180, radius = 0.99)
mid <- 180 + 360 * (cumsum(fractions) - fractions / 2)
text(0.62 * cos(mid * pi / 180), 0.62 * sin(mid * pi / 180),
     labels = sprintf("%.1f%%\n(n=%d)", counts$percentage, counts$environment_n),
     col = c("black", "white", "black", "white"), cex = 0.88)
par(fig = c(34 / 76, 1, 0, 1), mar = rep(0, 4), new = TRUE)
plot.new()
legend("center", legend = class_names, fill = class_colors, border = NA, bty = "n",
       cex = 1, x.intersp = 0.6, y.intersp = 1.05)
invisible(dev.off())

theme_panel <- theme_bw(base_size = 8) + theme(
  panel.grid = element_blank(),
  panel.border = element_rect(linewidth = 0.4),
  axis.text = element_text(color = "black", size = 7),
  legend.title = element_blank(),
  plot.margin = margin(2, 2, 2, 2, unit = "mm")
)
p_curve <- ggplot(curves, aes(x = bin, y = mean_value, color = factor(cluster), fill = factor(cluster))) +
  geom_ribbon(aes(ymin = q25_value, ymax = q75_value), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.65) +
  geom_vline(xintercept = r1_guides, linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = pattern_colors, labels = pattern_names) +
  scale_fill_manual(values = pattern_colors, labels = pattern_names) +
  scale_x_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Relative thermal progress (%)",
       y = expression(atop("Mean EDH"[30]~"intensity", "("*degree*C %.% h~"per hour)"))) +
  guides(fill = "none") + theme_panel +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.background = element_blank(),
        legend.text = element_text(size = 7), legend.key.height = unit(3, "mm"))
ggsave(file.path(out, "Extended_Data_b.pdf"), p_curve,
       width = 90, height = 62, units = "mm", device = cairo_pdf)

cluster_n <- assignment[, .N, by = cluster][order(cluster)]$N
cluster_centers <- c(610 - cluster_n[1] / 2, cluster_n[3] + cluster_n[2] / 2, cluster_n[3] / 2)
cluster_boundaries <- c(cluster_n[3] + 0.5, cluster_n[3] + cluster_n[2] + 0.5)
strip_dt <- unique(heat[, .(env_factor, cluster)])
p_strip <- ggplot(strip_dt, aes(x = 1, y = env_factor, fill = factor(cluster))) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = pattern_colors, guide = "none") +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
p_heat <- ggplot(heat, aes(x = bin, y = env_factor, fill = value_cap)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = r1_guides, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = cluster_boundaries, linewidth = 0.3) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 100.5), clip = "on") +
  scale_fill_gradient(low = "white", high = "firebrick3", limits = c(0, color_cap),
                      oob = squish, labels = label_number(accuracy = 0.1),
                      name = expression(atop("EDH"[30]~"intensity", "("*degree*C %.% h~"per hour)"))) +
  guides(fill = guide_colorbar(
    title.position = "right",
    barheight = unit(20, "mm"),
    barwidth = unit(2, "mm"),
    theme = theme(legend.title = element_text(angle = 270, size = 7))
  )) +
  annotate("text", x = c(75, 20, 27), y = cluster_centers,
           label = sprintf("%s\n(n=%d)", pattern_names, cluster_n), size = 2.4, lineheight = 0.95) +
  annotate("text", x = r1_guides + c(-1, 1), y = 6,
           label = sprintf("%.1f%%", r1_guides), hjust = c(1, 0), vjust = 0, size = 2.2) +
  labs(x = "Relative thermal progress (%)", y = NULL) +
  theme_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.title = element_text(size = 7), legend.text = element_text(size = 6),
        legend.key.height = unit(5, "mm"), legend.key.width = unit(2, "mm"),
        plot.margin = margin(0, 1, 0, 0, unit = "mm"))
# 色条与热图统一绘图区，避免侧条延伸到坐标标签区域。
p_heatmap <- p_strip + p_heat + plot_layout(widths = c(0.025, 0.975))
ggsave(file.path(out, "Extended_Data_c.pdf"), p_heatmap,
       width = 100, height = 90, units = "mm", device = cairo_pdf)

# 6. 英文 Source Data 只记录直接作图输入与计算结果；热图保留未截色的原值。
heatmap_data <- merge(order_dt, wide, by = "env_ID", sort = FALSE)
setorder(heatmap_data, display_order)
setnames(heatmap_data, bin_cols, sprintf("bin_%03d", 1:100))
source_sheets <- list(
  Panel_a_counts = counts,
  Panel_b_curves = curves,
  Panel_c_heatmap = heatmap_data,
  Environment_input = environments,
  R1_input = r1_input,
  R1_by_environment = r1_by_env,
  R1_quantiles = r1_quantiles
)
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
header <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold", fgFill = "#E8EDF2")
for (name in names(source_sheets)) {
  data <- source_sheets[[name]]
  text_values <- unlist(lapply(data, function(x) if (is.character(x)) x else character()))
  stopifnot(!any(grepl("[^\\x01-\\x7F]", text_values, perl = TRUE), na.rm = TRUE))
  addWorksheet(wb, name)
  writeData(wb, name, data, headerStyle = header, withFilter = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = 16)
  if ("env_ID" %in% names(data)) setColWidths(wb, name, cols = match("env_ID", names(data)), widths = 43)
}
saveWorkbook(wb, file.path(out, "Extended_Data_Source_Data.xlsx"), overwrite = TRUE)
message("Completed: 856 environments; 610 heat profiles; 236 R1 environments.")
print(counts)
print(r1_quantiles)

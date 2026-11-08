#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages(library(ggplot2))

# 1. 使用主图的完整站点十二年累计值；不改动归档数据及原图。
entry <- grep("^--file=", commandArgs(), value = TRUE)
out <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
regional <- read.csv(file.path(out, "regional_12year_cumulative_mean_sd.csv"))
sites <- read.csv(file.path(out, "site_12year_cumulative.csv"))
region_names <- c(
  NCP = "North China Plain",
  NEC = "Northeast China",
  NWC = "Northwest China",
  SWC = "Southwest China",
  NAMI = "North America (irrigated region)",
  NAMR = "North America (rainfed region)",
  EUR = "Europe",
  SAM = "South America",
  SAS = "South Asia",
  WAS = "West Asia and North Africa",
  EAF = "East Africa",
  WAF = "West Africa"
)
resolution <- 300
main_width_mm <- 297
height_mm <- 210
width_mm <- main_width_mm / 3
stem <- "Supplementary_N2O_warming_management"

# 2. 从站点累计值独立核对区域均值，分别保留两种情景的有效站点数。
checked <- aggregate(cumulative_N2O ~ region + mode + climate, sites, mean)
checked <- merge(regional, checked, by = c("region", "mode", "climate"))
stopifnot(nrow(checked) == 72L,
          max(abs(checked$mean - checked$cumulative_N2O)) < 1e-10)
extract <- function(mode, climate, prefix) {
  selected <- regional[regional$mode == mode & regional$climate == climate,
                       c("region", "mean", "n_sites")]
  stopifnot(nrow(selected) == 12L, !anyDuplicated(selected$region))
  names(selected)[2:3] <- paste0(prefix, c("_kg_N_ha", "_n_sites"))
  # 全局按站点等权汇总，不平均区域百分比；两情景沿用各自完整站点集合。
  pooled <- sites[sites$mode == mode & sites$climate == climate, ]
  stopifnot(!anyDuplicated(pooled$xSite))
  global <- data.frame("Overall", mean(pooled$cumulative_N2O), nrow(pooled))
  names(global) <- names(selected)
  rbind(selected, global)
}
results <- Reduce(function(x, y) merge(x, y, by = "region"), list(
  extract("F1_mzI2", "Baseline", "historical_BI"),
  extract("F1_mzI2", "Tmax+2", "warming_BI"),
  extract("F2_mzI2", "Tmax+2", "warming_SI")
))
stopifnot(all(results$historical_BI_kg_N_ha > 0))
results$region_name <- unname(region_names[results$region])
results$region_name[results$region == "Overall"] <- "Overall (all simulated sites, equal site weights)"
results$warming_response_pct <- 100 * (
  results$warming_BI_kg_N_ha / results$historical_BI_kg_N_ha - 1
)
results$SI_net_change_pct <- 100 * (
  results$warming_SI_kg_N_ha / results$historical_BI_kg_N_ha - 1
)
results$fully_offset <- results$SI_net_change_pct <= 0
results <- results[order(results$region != "Overall", -results$warming_response_pct), ]
results$rank <- seq_len(nrow(results))
row_count <- nrow(results)
stopifnot(row_count == 13L, results$region[1] == "Overall")
stopifnot(!anyNA(results), all(is.finite(results$SI_net_change_pct)))
write.csv(results, file.path(out, paste0(stem, ".csv")), row.names = FALSE)

# 3. 两条柱共用历史 BI 参照；按增温响应降序，从上至下排列。
warming <- data.frame(
  region = results$region,
  y = row_count + 1 - results$rank + .18,
  change = results$warming_response_pct,
  category = "Warming response (BI)"
)
net <- data.frame(
  region = results$region,
  y = row_count + 1 - results$rank - .18,
  change = results$SI_net_change_pct,
  category = "Net change with SI"
)
bars <- rbind(warming, net)
categories <- c("Warming response (BI)", "Net change with SI")
bars$category <- factor(bars$category, levels = categories)
bars$label <- sprintf("%+.1f%%", bars$change)
bars$label_x <- bars$change + ifelse(bars$change >= 0, .28, -.28)
bars$label_hjust <- ifelse(bars$change >= 0, 0, 1)
colors <- c("#F6C45C", "#275D63")
names(colors) <- categories
plot <- ggplot(bars) +
  geom_vline(xintercept = 0, linewidth = .45, colour = "#555555") +
  geom_rect(aes(xmin = pmin(0, change), xmax = pmax(0, change),
                ymin = y - .145, ymax = y + .145, fill = category)) +
  geom_text(aes(x = label_x, y = y, label = label, hjust = label_hjust),
            family = "Arial", fontface = "plain", size = 10 / ggplot2::.pt, colour = "#222222") +
  scale_fill_manual(values = colors, breaks = categories, drop = FALSE) +
  scale_y_continuous(breaks = row_count + 1 - results$rank, labels = results$region,
                     limits = c(.45, row_count + .55), expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 6),
                     expand = expansion(add = c(6, 18))) +
  labs(x = "Change in cumulative N2O-N (%)\nRelative to historical BI",
       y = NULL, fill = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    text = element_text(family = "Arial", face = "plain"),
    axis.text = element_text(size = 10, colour = "#222222"),
    axis.text.y = element_text(size = 11, margin = margin(r = 7)),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = .4),
    axis.ticks.x = element_line(linewidth = .4),
    axis.ticks.length = grid::unit(1.5, "mm"),
    legend.position = "inside",
    legend.position.inside = c(.99, .025),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 10),
    legend.key.size = grid::unit(3.5, "mm"),
    # 上边界对齐左图顶排面板，底边对齐底排横轴；底部留白包含横轴标题。
    plot.margin = margin(59.04, 8, 29.08, 5)
  )

# 字号调整后仍以物理尺寸固定绘图区上下边界，避免设备字体度量改变对齐。
ragg::agg_capture(width = width_mm, height = height_mm, units = "mm", res = resolution)
plot <- ggplotGrob(plot)
panel_row <- plot$layout[plot$layout$name == "panel", ]
above <- seq_len(panel_row$t - 1)
below <- seq.int(panel_row$b + 1, length(plot$heights))
top_mm <- grid::convertHeight(sum(plot$heights[above]), "mm", valueOnly = TRUE)
bottom_mm <- grid::convertHeight(sum(plot$heights[below]), "mm", valueOnly = TRUE)
plot$heights[1] <- plot$heights[1] + grid::unit(21 - top_mm, "mm")
last_row <- length(plot$heights)
plot$heights[last_row] <- plot$heights[last_row] + grid::unit(25.2 - bottom_mm, "mm")
invisible(dev.off())

# 4. 输出可编辑矢量图和 300 dpi 位图，避免字体转成逐字路径。
ggsave(file.path(out, paste0(stem, ".pdf")), plot,
       device = cairo_pdf, width = width_mm, height = height_mm, units = "mm")
ggsave(file.path(out, paste0(stem, ".svg")), plot,
       device = svglite::svglite, width = width_mm, height = height_mm, units = "mm")
ggsave(file.path(out, paste0(stem, ".png")), plot,
       device = ragg::agg_png, width = width_mm, height = height_mm, units = "mm",
       dpi = resolution, bg = "white")
ggsave(file.path(out, paste0(stem, ".tiff")), plot,
       device = ragg::agg_tiff, width = width_mm, height = height_mm, units = "mm",
       dpi = resolution, bg = "white", compression = "lzw")
print(results[, c("region", "warming_response_pct", "SI_net_change_pct", "fully_offset")],
      row.names = FALSE, digits = 4)
cat("Verified: Overall plus 12 regions; 26 bars; percentage changes relative to historical BI.\n")

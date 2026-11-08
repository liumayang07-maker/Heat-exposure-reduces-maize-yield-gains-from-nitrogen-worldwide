# Figure1d heatmap, trimmed from original script 07; no palette/QC/cache exports.
#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))

standardize_rows <- function(values) {
  if (any(!is.finite(values))) stop("Missing or nonfinite values cannot be standardized.")
  spread <- apply(values, 1, sd)
  if (any(!is.finite(spread) | spread <= 0)) stop("Constant rows have undefined Z-scores.")
  sweep(sweep(values, 1, rowMeans(values), "-"), 1, spread, "/")
}

main <- function() {
  # 1. 读取完整的 41 年输入，保留原区域顺序用于聚类。
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: Rscript plot_heatmap.R INPUT_CSV REGIONS_CSV OUTPUT_DIR")
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(pheatmap))
  data <- fread(args[1])
  styles <- fread(args[2])[order(order)]
  output_dir <- args[3]
  if (!dir.exists(output_dir)) stop("Output directory does not exist.")
  years <- 1985:2025
  stopifnot(nrow(styles) == 15L, !anyDuplicated(styles$region_id),
            nrow(data) == 615L, !anyDuplicated(data[, .(region_id, season_year)]),
            setequal(data$region_id, styles$region_id), all(data$season_year %in% years),
            all(is.finite(data$area_weighted_mean)), all(data$area_weighted_mean >= 0),
            all(data$area_source_year == pmin(data$season_year, 2020)),
            all(data[, .N, by = region_id]$N == 41L))
  wide <- dcast(data, region_id ~ season_year, value.var = "area_weighted_mean")
  wide <- wide[match(styles$region_id, wide$region_id)]
  values <- as.matrix(wide[, as.character(years), with = FALSE])
  rownames(values) <- styles$region_id
  z <- standardize_rows(values)
  stopifnot(inherits(try(standardize_rows(matrix(1, 2, 4)), silent = TRUE), "try-error"),
            inherits(try(standardize_rows(matrix(NA_real_, 2, 4)), silent = TRUE), "try-error"),
            max(abs(rowMeans(z))) < 1e-12, max(abs(apply(z, 1, sd) - 1)) < 1e-12,
            max(abs(t(scale(t(values))) - z)) < 1e-12)
  # 2. 区域内标准化后聚类区域，年份保持时间顺序。
  tree <- hclust(
    dist(z, method = "euclidean"),
    method = "ward.D2"
  )
  ordered_ids <- rownames(z)[tree$order]
  output <- file.path(output_dir, "Figure1d_heatmap.pdf")
  if (file.exists(output)) stop("Output exists.")
  limit <- max(4, ceiling(max(abs(z))))
  colors <- colorRampPalette(c("#053061", "#08519C", "#2171B5", "#4393C3", "#73ACCF",
                              "#FFFFFF", "#E58E79", "#D6604D", "#B2182B", "#8B0A1A",
                              "#67001F"))(1001)
  breaks <- seq(-limit, limit, length.out = length(colors) + 1L)
  legend_ticks <- seq(-5, 5, by = 2.5)
  stopifnot(colors[501] == "#FFFFFF", min(z) >= min(breaks), max(z) <= max(breaks))

  row_labels <- styles$region_name
  row_labels[styles$region_id == "CHN_NE"] <- "Northeast China"
  row_labels[styles$region_id == "CHN_NCP"] <- "North China Plain"
  row_labels[styles$region_id == "CHN_SW"] <- "Southwest China"
  row_labels[styles$region_id == "CHN_NW"] <- "Northwest China"
  ordered_labels <- row_labels[match(ordered_ids, styles$region_id)]
  plot_matrix <- t(z[ordered_ids, , drop = FALSE])
  colnames(plot_matrix) <- ordered_labels
  rownames(plot_matrix) <- as.character(years)

  draw <- function() {
    # 3. 保留原版布局，只调整图形对象的位置，不改变数据。
    result <- pheatmap(plot_matrix, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                      color = colors, breaks = breaks, border_color = "#D9D9D9",
                      labels_row = ifelse(years %% 5 == 0, as.character(years), ""),
                      fontsize_row = 4.5, labels_col = ordered_labels,
                      fontsize_col = 4.5, angle_col = "90", fontsize = 5,
                      treeheight_row = 0, treeheight_col = 0,
                      legend_breaks = legend_ticks, legend_labels = as.character(legend_ticks),
                      main = NA, fontfamily = "sans", silent = TRUE)
    stopifnot(identical(result$tree_row, NA), identical(result$tree_col, NA),
              identical(colnames(plot_matrix), ordered_labels),
              identical(rownames(plot_matrix), as.character(years)))
    # Move the label column into the unused left dendrogram space.
    label_index <- which(result$gtable$layout$name == "row_names")
    label_column <- result$gtable$layout$l[label_index]
    left_width <- result$gtable$widths[1]
    result$gtable$widths[1] <- result$gtable$widths[label_column]
    result$gtable$widths[label_column] <- left_width
    result$gtable$layout$l[label_index] <- 1L
    result$gtable$layout$r[label_index] <- 1L
    result$gtable$grobs[[label_index]]$x <- grid::unit(1, "npc") - grid::unit(3, "bigpts")
    result$gtable$grobs[[label_index]]$hjust <- 1
    # Attach ticks at the exact label positions calculated by pheatmap.
    legend_index <- which(result$gtable$layout$name == "legend")
    legend_grob <- result$gtable$grobs[[legend_index]]
    # Compress only the bar and label positions, preserving font sizes.
    for (i in 1:2) {
      legend_grob$children[[i]]$y <- grid::unit(1, "npc") -
        (grid::unit(1, "npc") - legend_grob$children[[i]]$y) / 3
    }
    legend_grob$children[[1]]$height <- legend_grob$children[[1]]$height / 3
    legend_grob$children[[1]]$width <- grid::unit(5, "bigpts")
    legend_grob$children[[2]]$x <- grid::unit(9, "bigpts")
    tick_y <- legend_grob$children[[2]]$y
    stopifnot(length(tick_y) == length(legend_ticks), all(diff(legend_ticks) == 2.5))
    legend_grob <- grid::addGrob(legend_grob, grid::segmentsGrob(
      x0 = grid::unit(5, "bigpts"), x1 = grid::unit(7.5, "bigpts"),
      y0 = tick_y, y1 = tick_y, gp = grid::gpar(col = "#303030", lwd = .6)))
    result$gtable$grobs[[legend_index]] <- legend_grob
    tiles <- result$gtable$grobs[[which(result$gtable$layout$name == "matrix")]]
    stopifnot(length(tiles$children[[1]]$gp$fill) == 615L)
    tiles$children[[1]]$gp$lwd <- .25
    tiles <- grid::addGrob(tiles, grid::rectGrob(
      gp = grid::gpar(fill = NA, col = "#303030", lwd = .7)))
    result$gtable$grobs[[which(result$gtable$layout$name == "matrix")]] <- tiles
    grid::grid.newpage()
    grid::grid.draw(result$gtable)
  }
  pdf(output, width = 78 / 25.4, height = 42 / 25.4, family = "sans", bg = "white",
      useDingbats = FALSE)
  tryCatch(draw(), finally = dev.off())

  data[, z_score := z[cbind(match(region_id, rownames(z)), match(season_year, years))]]
  # 4. 返回实际热图值；主脚本写入英文 Source Data 后清理临时文件。
  data[, row_order := match(region_id, ordered_ids)]
  setorder(data, row_order, season_year)
  fwrite(data[, .(region_id, season_year, area_weighted_mean, z_score, row_order)],
         file.path(output_dir, "Heatmap_data.csv"))

}

main()

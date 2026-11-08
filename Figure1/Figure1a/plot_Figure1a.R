#!/usr/bin/env Rscript
# Figure 1a: multi-year mean annual EDH30 in maize-area cells, 1986-2025.
# Run: Rscript plot_Figure1a.R [--overwrite]
# Requires terra, sf, ggplot2 and Cairo. Inputs and output are beside this script.
# EDH30 = sum(max(hourly temperature in C - 30, 0)), in degree C hours.
# Full UTC calendar years; not growing-season or hot-day-normalized exposure.
# Mean uses available maize-masked years; per-pixel year counts are unavailable.
# No additional area weighting, interpolation, imputation or exposure calculation.
# Upstream harvested-area data are not bundled; see README.md for the source.
# Use Harvested_area_1961-2020.h5, dataset Maize; retain maize_area > 0 for masks.
# For 2021-2025 reuse 2020 area. Plotting this supplied mean needs no area data.
# Raster: mean_annual_EDH30 band, losslessly extracted from
# global_maize_edh_1980_2025/data/derived/grid_year_stage_metrics/
# gee_era5land_multiyear_mean_annual_heat_maize_area_1986_2025.tif.
# Boundaries: Natural Earth ne_110m_admin_0_countries (public domain),
# geometry-only conversion to GeoPackage; no simplification or reprojection.
# Adapted from Figure2_v3/01_data_heat_patterns/scripts/
# plot_global_EDH30_meta_sites.R --edh-only; unrelated meta-study inputs removed.
# Retains the source panel's discrete palette/layout, not the assembled JPG's
# manually edited legend and typography. Heat is raster; borders/text are vector.

main <- function() {
  if (.Platform$OS.type == 'windows') invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
  args <- commandArgs(trailingOnly = TRUE)
  if (any(!args %in% '--overwrite')) stop('Usage: Rscript plot_Figure1a.R [--overwrite]')
  entry <- grep('^--file=', commandArgs(), value = TRUE)
  if (length(entry) != 1L) stop('Run with Rscript.')
  root <- dirname(normalizePath(sub('^--file=', '', entry), winslash = '/'))
  output <- file.path(root, 'Figure1a.pdf')
  if (file.exists(output) && !'--overwrite' %in% args) stop('PDF exists; use --overwrite explicitly.')
  for (pkg in c('terra', 'sf', 'ggplot2')) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop('Missing package: ', pkg)
  }
  if (!capabilities('cairo')) stop('Cairo PDF support is required.')
  suppressPackageStartupMessages(library(ggplot2))
  raster <- terra::rast(file.path(root, 'mean_annual_EDH30_1986_2025.tif'))
  stopifnot(terra::nlyr(raster) == 1L, names(raster) == 'mean_annual_EDH30',
            isTRUE(sf::st_crs(terra::crs(raster)) == sf::st_crs(4326)))
  values <- terra::values(raster, mat = FALSE)
  stopifnot(any(is.finite(values)), !any(is.infinite(values)),
            all(values[is.finite(values)] >= 0))
  world <- sf::st_read(file.path(root, 'world_boundaries.gpkg'), quiet = TRUE)
  stopifnot(isTRUE(sf::st_crs(world) == sf::st_crs(4326)))
  labels <- c('0', '(0, 100]', '(100, 200]', '(200, 400]', '(400, 800]',
              '(800, 1,600]', '(1,600, 3,200]', '(3,200, 6,400]',
              '(6,400, 12,800]', '>12,800')
  legend_labels <- c('0', '0-\n100', '100-\n200', '200-\n400', '400-\n800',
                     '0.8-\n1.6k', '1.6-\n3.2k', '3.2-\n6.4k', '6.4-\n12.8k', '>12.8k')
  palette <- setNames(c('#68BDC4', '#BFE3C7', '#FFF7BC', '#FEE391', '#FEC44F',
                        '#FE9929', '#F03B20', '#D7191C', '#B10026', '#67000D'), labels)
  classify <- function(x) {
    ans <- rep(NA_integer_, length(x))
    ans[is.finite(x) & x == 0] <- 1L
    positive <- is.finite(x) & x > 0
    ans[positive] <- as.integer(cut(x[positive],
      c(0, 100, 200, 400, 800, 1600, 3200, 6400, 12800, Inf), right = TRUE)) + 1L
    ans
  }
  stopifnot(identical(classify(c(NA, 0, 0.1, 100, 100.1, 12800, 12801)),
                      c(NA_integer_, 1L, 2L, 2L, 3L, 9L, 10L)))
  local <- terra::crop(raster, terra::ext(-130, 155, -60, 85), snap = 'out')
  classes <- classify(terra::values(local, mat = FALSE))
  colors <- unname(palette[classes])
  colors[is.na(classes)] <- 'transparent'
  bitmap <- as.raster(matrix(colors, nrow = nrow(local), ncol = ncol(local), byrow = TRUE))
  extent <- as.vector(terra::ext(local))
  legend <- data.frame(x = 0, y = 0, class_id = factor(labels, levels = labels))
  p <- ggplot() +
    geom_sf(data = world, fill = '#EEEEEE', color = NA) +
    annotation_raster(bitmap, xmin = extent[1], xmax = extent[2],
                      ymin = extent[3], ymax = extent[4], interpolate = FALSE) +
    geom_sf(data = world, fill = NA, color = '#8E8E8E', linewidth = 0.17) +
    geom_tile(data = legend, aes(x, y, fill = class_id), width = 0, height = 0) +
    scale_fill_manual(values = palette, breaks = labels, labels = legend_labels,
                      drop = FALSE, name = 'Annual EDH30 (degree C h)') +
    scale_x_continuous(breaks = seq(-120, 120, 60)) +
    scale_y_continuous(breaks = seq(-60, 60, 30)) +
    coord_sf(crs = sf::st_crs(4326), default_crs = sf::st_crs(4326),
             xlim = c(-130, 155), ylim = c(-60, 85), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE,
      title.position = 'top', title.hjust = 0.5, label.position = 'bottom',
      label.hjust = 0.5, keywidth = grid::unit(13, 'mm'), keyheight = grid::unit(4, 'mm'))) +
    theme_minimal(base_size = 10, base_family = 'Arial') +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white', color = NA),
      panel.border = element_rect(fill = NA, color = '#777777', linewidth = 0.35),
      axis.text = element_text(size = 7.5, color = '#404040'),
      axis.ticks = element_line(color = '#404040', linewidth = 0.25),
      axis.ticks.length = grid::unit(1.5, 'mm'), legend.position = 'bottom',
      legend.box = 'vertical', legend.title = element_text(size = 9),
      legend.text = element_text(size = 7.2, lineheight = 0.85),
      legend.key.width = grid::unit(8, 'mm'), legend.key.height = grid::unit(4, 'mm'),
      legend.spacing.x = grid::unit(0, 'mm'), legend.key.spacing.x = grid::unit(0, 'mm'),
      plot.margin = margin(4, 5, 4, 5))
  grDevices::cairo_pdf(output, width = 155 / 25.4, height = 6.3, family = 'Arial', onefile = TRUE)
  tryCatch(print(p), finally = grDevices::dev.off())
  stopifnot(file.info(output)$size > 0)
  message('Created Figure1a.pdf; R ', getRversion(), '; terra ', packageVersion('terra'),
          '; sf ', packageVersion('sf'), '; ggplot2 ', packageVersion('ggplot2'))
}

main()

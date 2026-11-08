#!/usr/bin/env Rscript
# Figure1b: six independent GLS AR(1)/REML fits to unsmoothed annual means.
# Run: Rscript plot_Figure1b.R [--overwrite]. See README.md for provenance.
# Adapted from Figure1_v3/scripts/plot_continent_annual_EDH30_GLS.R.
options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  export_arg <- grep("^--source-data-dir=", args, value = TRUE)
  if (length(export_arg) > 1L || any(!args %in% "--overwrite" & !grepl("^--source-data-dir=", args))) stop("Invalid arguments.")
  entry <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(entry) != 1L) stop("Run with Rscript.")
  root <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
  input <- file.path(root, "continent_annual_EDH30_1986_2025.csv")
  required <- c("data.table", "ggplot2", "nlme")
  if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) stop("Missing R packages.")
  if (!capabilities("cairo")) stop("Cairo PDF support required.")
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ggplot2))
  colors <- c("North America" = "#0072B2", "South America" = "#009E73",
              "Europe" = "#56B4E9", "Asia" = "#D55E00",
              "Africa" = "#E69F00", "Oceania" = "#CC79A7")
  annual <- fread(input)
  stopifnot(nrow(annual) == 240L, setequal(unique(annual$continent), names(colors)),
            !anyDuplicated(annual[, .(continent, year)]),
            all(is.finite(annual$area_weighted_mean)), all(annual$area_weighted_mean >= 0),
            all(annual$area_source_year == pmin(annual$year, 2020L)),
            all(annual$valid_area_share > 0 & annual$valid_area_share <= 1))
  pdf_path <- file.path(root, "Figure1b.pdf")
  model_path <- file.path(root, "GLS_models.csv")
  if (!length(export_arg) && !"--overwrite" %in% args && any(file.exists(c(pdf_path, model_path)))) {
    stop("Outputs exist; use --overwrite explicitly.")
  }
  models <- list()
  fitted_data <- list()
  parameters <- list()
  covariance <- list()
  for (region in names(colors)) {
    d <- annual[continent == region][order(year)]
    stopifnot(identical(d$year, 1986:2025))
    d[, year_c := year - mean(year)]
    fit <- nlme::gls(area_weighted_mean ~ year_c, data = d,
                     correlation = nlme::corAR1(form = ~ year), method = "REML")
    tab <- summary(fit)$tTable
    beta <- unname(coef(fit)["year_c"])
    se <- tab["year_c", "Std.Error"]
    df <- nrow(d) - length(coef(fit))
    phi <- unname(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
    parameters[[region]] <- data.table(continent = region, year_center = mean(d$year),
      intercept = unname(coef(fit)["(Intercept)"]), slope_per_year = beta,
      residual_sigma = fit$sigma, ar1_phi = phi, df = df)
    v <- vcov(fit)
    covariance[[region]] <- data.table(continent = region,
      term_row = rep(rownames(v), times = ncol(v)),
      term_col = rep(colnames(v), each = nrow(v)), covariance = as.vector(v))
    stopifnot(is.finite(beta), is.finite(se), se > 0, abs(phi) < 1)
    d[, gls_fit := as.numeric(predict(fit, newdata = d))]
    design <- model.matrix(~ year_c, data = d)
    d[, gls_mean_se := sqrt(rowSums((design %*% vcov(fit)) * design))]
    d[, `:=`(gls_ci_low = gls_fit - qt(.975, df) * gls_mean_se,
             gls_ci_high = gls_fit + qt(.975, df) * gls_mean_se)]
    stopifnot(all(is.finite(d$gls_mean_se)), all(d$gls_ci_low <= d$gls_fit),
              all(d$gls_ci_high >= d$gls_fit))
    stopifnot(max(abs(d$gls_fit - as.numeric(fitted(fit)))) < 1e-8)
    fitted_data[[region]] <- d[, .(continent, year, area_source_year, area_weighted_mean,
                                   valid_area_share, gls_fit, gls_mean_se, gls_ci_low, gls_ci_high)]
    models[[region]] <- data.table(continent = region, model = "GLS AR(1), REML", n_years = nrow(d),
                                   slope_per_decade = beta * 10, se_per_decade = se * 10,
                                   ci_low_per_decade = (beta - qt(.975, df) * se) * 10,
                                   ci_high_per_decade = (beta + qt(.975, df) * se) * 10,
                                   p_value = tab["year_c", "p-value"], df = df, ar1_phi = phi)
  }
  model_table <- rbindlist(models)
  model_table[, p_BH_six_continents := p.adjust(p_value, method = "BH")]
  p_text <- ifelse(model_table$p_BH_six_continents < .001, "P(BH) < 0.001",
                   sprintf("P(BH) = %.3f", model_table$p_BH_six_continents))
  legend_labels <- setNames(sprintf("%s\n%+.1f; %s", model_table$continent,
                                    model_table$slope_per_decade, p_text), model_table$continent)
  plot_data <- rbindlist(fitted_data)
  if (length(export_arg)) {
    export_dir <- sub("^--source-data-dir=", "", export_arg)
    stopifnot(dir.exists(export_dir))
    tables <- list(Annual_data = annual, Plot_data = plot_data, GLS_models = model_table,
                   Model_parameters = rbindlist(parameters), Coefficient_covariance = rbindlist(covariance))
    targets <- file.path(export_dir, paste0(names(tables), ".csv"))
    if (any(file.exists(targets))) stop("Source tables already exist.")
    for (i in seq_along(tables)) fwrite(tables[[i]], targets[i])
    return(invisible(NULL))
  }
  plot_data[, continent := factor(continent, levels = names(colors))]
  p <- ggplot(plot_data, aes(year, color = continent)) +
    geom_ribbon(aes(ymin = gls_ci_low, ymax = gls_ci_high, fill = continent),
                color = NA, alpha = .14, show.legend = FALSE) +
    geom_line(aes(y = area_weighted_mean, linetype = "Annual mean"), linewidth = .65, alpha = .65) +
    geom_line(aes(y = gls_fit, linetype = "GLS trend"), linewidth = 1.05) +
    scale_color_manual(values = colors, labels = legend_labels,
                       name = expression("GLS AR(1), REML; slope: "*degree*C~h~decade^{-1}*"; shading: 95% mean CI")) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_linetype_manual(values = c("Annual mean" = "solid", "GLS trend" = "22"), name = NULL) +
    scale_x_continuous(breaks = seq(1990, 2025, 5), limits = c(1985.5, 2025.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 5000, 1000),
                       expand = expansion(mult = c(0, .04))) +
    labs(x = "Year", y = expression("Area-weighted annual EDH30 ("*degree*C~h~year^{-1}*")")) +
    guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE,
                                override.aes = list(alpha = 1, linewidth = 1, linetype = "solid")),
           linetype = guide_legend(order = 2, override.aes = list(color = "#333333", alpha = 1))) +
    theme_classic(base_size = 7, base_family = "Arial") +
    theme(legend.position = "none", legend.position.inside = c(.015, .985),
          legend.box = "vertical", legend.justification = c(0, 1), legend.box.just = "left",
          legend.key.width = grid::unit(4, "mm"), legend.key.height = grid::unit(2, "mm"),
          legend.margin = margin(0, 0, 0, 0), legend.spacing.y = grid::unit(.4, "mm"),
          legend.text = element_text(size = 5.3), legend.title = element_text(size = 5.5),
          legend.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(size = 6, color = "#303030"), panel.grid = element_blank(),
          axis.title = element_text(size = 6.5), axis.line = element_blank(),
          panel.border = element_rect(fill = NA, color = "#303030", linewidth = .3),
          plot.margin = margin(2, 8, 2, 2))
  ggsave(pdf_path, p, width = 89, height = 81, units = "mm",
         device = grDevices::cairo_pdf, bg = "white")
  fwrite(model_table, model_path)
  stopifnot(all(file.info(c(pdf_path, model_path))$size > 0))
  print(model_table)
  message("Verified: six GLS fits; 240 observations; only PDF and model CSV written.")
}

main()

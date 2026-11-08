#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(metafor)
})

# 1. 默认读取最简归档输入；传入主表时按原区域规则重新提取，不改写主表。
entry <- grep("^--file=", commandArgs(), value = TRUE)
out <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript fit_Figure3b.R [MASTER_CSV]")
region_levels <- c(
  "North America_Irrigated", "North America_Rainfed", "Europe", "WANA",
  "Northwest China", "Northeast China", "South Asia", "North China Plain",
  "Southwest China", "East Africa", "West Africa", "South America"
)
region_abbr <- setNames(
  c("NAMI", "NAMR", "EUR", "WAS", "NWC", "NEC", "SAS", "NCP", "SWC", "EAF", "WAF", "SAM"),
  region_levels
)
input_columns <- c(
  "data_ID", "env_ID", "reference_ID", "site_year_ID", "Label_of_control_treatment",
  "mean_c", "sd_c", "n_c", "mean_t", "sd_t", "n_t", "lnRR", "region",
  "country", "EDH30_total", "N_rate", "Exp_year", "GDD_total", "soil_pH",
  "SOC", "soil_total_N", "total_rainfall", "plant_density", "N_strategy", "Irrigation_decision"
)
input_file <- file.path(out, "comparison_input.csv")
if (length(args)) {
  text <- iconv(readLines(args[[1]], warn = FALSE), from = "GB18030", to = "UTF-8")
  if (anyNA(text)) stop("Cannot decode master CSV as GB18030.")
  master <- fread(text = paste(text, collapse = "\n"))
  stopifnot(nrow(master) == 4613L, !anyDuplicated(master$data_ID))
  fwrite(master[region %chin% region_levels, ..input_columns], input_file)
}
d <- fread(input_file, encoding = "UTF-8")
stopifnot(
  nrow(d) == 4541L,
  uniqueN(d$reference_ID) == 277L,
  uniqueN(d$env_ID) == 846L,
  !anyDuplicated(d$data_ID),
  !anyNA(d[, ..input_columns]),
  setequal(d$region, region_levels),
  all(d$mean_t > 0 & d$mean_c > 0),
  all(d$sd_t >= 0 & d$sd_c >= 0),
  all(d$n_t > 0 & d$n_c > 0),
  max(abs(d$lnRR - log(d$mean_t / d$mean_c))) < 1e-8
)
tab <- function(x, name) fwrite(x, file.path(out, paste0(name, ".csv")))

# 2. 用当前 SD 重建稀疏采样协方差，不读取旧 design.rds 或已拟合模型。
d[, yi := lnRR]
d[, treatment_var := sd_t^2 / (n_t * mean_t^2)]
d[, control_var := sd_c^2 / (n_c * mean_c^2)]
d[, vi := treatment_var + control_var]
d[, control_id := paste(
  reference_ID,
  site_year_ID,
  Label_of_control_treatment,
  sprintf("%.10g", mean_c),
  sprintf("%.10g", sd_c),
  sprintf("%.10g", n_c),
  env_ID,
  sep = "||"
)]
d[, model_row := .I]
blocks <- split(d$model_row, d$control_id)
shared <- blocks[lengths(blocks) > 1L]
pairs <- do.call(cbind, lapply(shared, combn, m = 2))
stopifnot(all(is.finite(d$vi) & d$vi > 0))
stopifnot(all(vapply(blocks, function(i) max(abs(d$control_var[i] - d$control_var[i[1]])) < 1e-12, logical(1))))
V <- sparseMatrix(
  i = c(d$model_row, pairs[1, ], pairs[2, ]),
  j = c(d$model_row, pairs[2, ], pairs[1, ]),
  x = c(d$vi, d$control_var[pairs[1, ]], d$control_var[pairs[1, ]]),
  dims = c(nrow(d), nrow(d))
)
stopifnot(isSymmetric(V), max(abs(diag(V) - d$vi)) < 1e-12)

# 3. 沿用原区域模型的变换、参考水平、区域内环境等权热暴露中心化。
transformations <- fread(file.path(out, "transformation_parameters.csv"))
for (i in seq_len(nrow(transformations))) {
  specification <- transformations[i]
  value <- d[[specification$source]]
  if (specification$kind == "log1p") value <- log1p(value)
  set(d, j = specification$target, value = (value - specification$center) / specification$divisor)
}
d[, `:=`(N_c_sq = N_c^2, rain_c_sq = rain_c^2)]
d[, region := relevel(factor(region, levels = region_levels), "North China Plain")]
d[, country := factor(country)]
d[, irrigation := factor(Irrigation_decision, levels = c("Rainfed", "Irrigation"))]
d[, N_strategy := factor(N_strategy, levels = c("S", "SV", "SVR", "SR", "VR", "V", "R"))]
d[, plant_density_c := plant_density - mean(plant_density)]
d[, ref_site_year_ID := paste(reference_ID, site_year_ID, sep = "__")]
d[, effect_id := data_ID]
for (variable in c("reference_ID", "ref_site_year_ID", "effect_id")) set(d, j = variable, value = factor(d[[variable]]))
env_heat <- unique(d[, .(env_ID, region, EDH30_total)])
stopifnot(nrow(env_heat) == uniqueN(d$env_ID))
env_heat[, heat_log := log1p(EDH30_total)]
env_heat[, heat_region_mean := mean(heat_log), by = region]
env_heat[, heat_wr := heat_log - heat_region_mean]
d[env_heat, on = .(env_ID), `:=`(heat_log = i.heat_log, heat_region_mean = i.heat_region_mean, heat_wr = i.heat_wr)]

# 国家效应嵌套于区域；参考国家按环境数、文献数及国家名称的原版顺序选择。
nested_country_terms <- character()
country_rows <- list()
for (reg in levels(d$region)) {
  in_region <- as.character(d$region) == reg
  country_support <- d[in_region, .(
    effects = .N, references = uniqueN(reference_ID), environments = uniqueN(env_ID)
  ), by = country][order(-environments, -references, as.character(country))]
  reference_country <- as.character(country_support$country[1])
  for (cty in as.character(country_support$country)) {
    country_rows[[length(country_rows) + 1L]] <- data.table(
      region = reg, region_abbr = unname(region_abbr[reg]), country = cty,
      reference_country = cty == reference_country,
      effects = country_support[as.character(country) == cty, effects],
      references = country_support[as.character(country) == cty, references],
      environments = country_support[as.character(country) == cty, environments]
    )
    if (cty != reference_country) {
      term <- make.names(paste("country_within", unname(region_abbr[reg]), cty, sep = "_"))
      term <- make.unique(c(nested_country_terms, term))[length(nested_country_terms) + 1L]
      set(d, j = term, value = as.integer(in_region & as.character(d$country) == cty))
      nested_country_terms <- c(nested_country_terms, term)
    }
  }
}
tab(rbindlist(country_rows), "nested_country_mapping")
support <- d[, .(
  effects = .N,
  references = uniqueN(reference_ID),
  environments = uniqueN(env_ID),
  countries = uniqueN(country),
  N_min = min(N_rate),
  N_q25 = as.numeric(quantile(unique(.SD[, .(env_ID, N_rate)])$N_rate, .25)),
  N_q50 = as.numeric(quantile(unique(.SD[, .(env_ID, N_rate)])$N_rate, .50)),
  N_q75 = as.numeric(quantile(unique(.SD[, .(env_ID, N_rate)])$N_rate, .75)),
  N_max = max(N_rate),
  EDH30_min = min(EDH30_total),
  EDH30_max = max(EDH30_total)
), by = region]
support[, region_abbr := unname(region_abbr[as.character(region)])]
tab(support, "regional_support")
rhs <- c(
  "region * heat_wr * (N_c + N_c_sq)", nested_country_terms,
  "plant_density_c", "N_strategy", "year_c", "GDD_c", "pH_c",
  "SOC_log_c", "TN_log_c", "rain_c", "rain_c_sq", "irrigation"
)
form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
matrix_for <- function(formula, data) {
  X <- model.matrix(formula, data)
  colnames(X)[colnames(X) == "(Intercept)"] <- "intrcpt"
  X
}
X <- matrix_for(form, d)
stopifnot(!anyNA(d), ncol(X) == 108L, qr(X)$rank == ncol(X))
tab(d, "model_input")

# 4. 仅拟合主图使用的国家调整模型，三个随机截距及优化器保持不变。
message("Fitting Figure3b: 4541 comparisons; 108 fixed-effect coefficients; REML.")
fit <- rma.mv(
  yi = yi,
  V = V,
  mods = form,
  random = list(~ 1 | reference_ID, ~ 1 | ref_site_year_ID, ~ 1 | effect_id),
  data = d,
  method = "REML",
  sparse = TRUE,
  control = list(optimizer = "nlminb")
)
stopifnot(fit$k == nrow(d), all(is.finite(coef(fit))))
saveRDS(list(fit = fit, formula = form), file.path(out, "model_REML.rds"))
tab(data.table(
  term = names(coef(fit)), estimate = as.numeric(coef(fit)), se = fit$se,
  ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, p = fit$pval
), "coefficients_REML")
covariance <- as.data.table(as.table(vcov(fit)))
setnames(covariance, c("term_1", "term_2", "covariance"))
tab(covariance, "coefficient_covariance_REML")
tab(data.table(level = c("reference", "reference_site_year", "effect"), sigma2 = fit$sigma2), "variance_components_REML")

# 5. 沿用原版 QR Wald 算法及 12 个区域为一组的 BH 校正。
wald_from_L <- function(L) {
  estimate <- as.numeric(L %*% coef(fit))
  variance <- L %*% vcov(fit) %*% t(L)
  degrees <- qr(variance)$rank
  statistic <- as.numeric(crossprod(estimate, qr.solve(variance, estimate)))
  data.table(QM = statistic, df = degrees, p = pchisq(statistic, degrees, lower.tail = FALSE))
}
coefficient_names <- names(coef(fit))
three_way_terms <- which(vapply(strsplit(coefficient_names, ":", fixed = TRUE), function(parts) {
  any(startsWith(parts, "region")) && "heat_wr" %in% parts && any(parts %in% c("N_c", "N_c_sq"))
}, logical(1)))
stopifnot(length(three_way_terms) == 22L)
template <- copy(d[1])
for (term in nested_country_terms) set(template, j = term, value = 0L)
template[, `:=`(
  plant_density_c = 0, year_c = 0, GDD_c = 0, pH_c = 0,
  SOC_log_c = 0, TN_log_c = 0, rain_c = 0, rain_c_sq = 0,
  irrigation = factor("Rainfed", levels = levels(d$irrigation)),
  N_strategy = factor("S", levels = levels(d$N_strategy))
)]
interaction_contrasts <- function(reg) {
  model_row <- function(h, n1, n2) {
    nd <- copy(template)
    nd[, region := factor(reg, levels = levels(d$region))]
    nd[, `:=`(heat_wr = h, N_c = n1, N_c_sq = n2)]
    matrix_for(form, nd)
  }
  rbind(
    model_row(1, 1, 0) - model_row(1, 0, 0) - model_row(0, 1, 0) + model_row(0, 0, 0),
    model_row(1, 0, 1) - model_row(1, 0, 0) - model_row(0, 0, 1) + model_row(0, 0, 0)
  )
}
regional_tests <- rbindlist(lapply(region_levels, function(reg) {
  cbind(data.table(region = reg, region_abbr = unname(region_abbr[reg])), wald_from_L(interaction_contrasts(reg)))
}))
regional_tests[, p_BH := p.adjust(p, method = "BH")]
tab(regional_tests, "regional_heat_N_joint_tests_REML")
tab(rbindlist(list(
  cbind(data.table(test = "Equality of heat-by-N interactions across regions"),
        wald_from_L(diag(length(coefficient_names))[three_way_terms, , drop = FALSE])),
  cbind(data.table(test = "All region-specific heat-by-N interactions equal zero"),
        wald_from_L(do.call(rbind, lapply(region_levels, interaction_contrasts))))
)), "omnibus_Wald_tests_REML")
tab(data.table(
  effects = fit$k,
  references = uniqueN(d$reference_ID),
  environments = uniqueN(d$env_ID),
  regions = uniqueN(d$region),
  control_groups = length(blocks),
  shared_control_groups = length(shared),
  V_nonzero = nnzero(V),
  design_columns = ncol(X),
  design_rank = qr(X)$rank,
  formula = paste(deparse(form, width.cutoff = 500), collapse = " "),
  method = fit$method,
  optimizer = "nlminb",
  region_reference = "North China Plain",
  irrigation_reference = "Rainfed",
  N_strategy_reference = "S",
  plant_density_center = mean(d$plant_density),
  prediction_weighting = "Equal weight per environment within region; shared within-environment rows",
  N_scenarios = "Regional Q25/Q50/Q75 of unique env_ID and N_rate pairs; type 7",
  exposure_centering = "Mean log1p(EDH30) over unique environments within region",
  regional_p_adjustment = "BH across 12 regional two-df joint tests",
  logLik_REML = as.numeric(logLik(fit)),
  AIC_REML = AIC(fit),
  BIC_REML = BIC(fit),
  QE = fit$QE,
  QE_p = fit$QEp,
  R_version = as.character(getRversion()),
  metafor_version = as.character(packageVersion("metafor"))
), "model_parameters")
message("Completed Figure3b fit. Elapsed fitting seconds: ", round(fit$time, 2))
print(regional_tests)

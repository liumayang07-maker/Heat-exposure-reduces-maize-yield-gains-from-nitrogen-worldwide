#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(metafor)
})

# 1. 默认读取归档输入；可传入最新主表路径，重新提取本面板所需字段。
entry <- grep("^--file=", commandArgs(), value = TRUE)
out <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript fit_Figure2b.R [MASTER_CSV]")
input_file <- file.path(out, "comparison_input.csv")
input_columns <- c(
  "data_ID", "env_ID", "reference_ID", "site_year_ID", "Label_of_control_treatment",
  "mean_c", "sd_c", "n_c", "mean_t", "sd_t", "n_t", "lnRR",
  "EDH30_total", "N_rate", "N_rate_cluster", "plant_density", "N_strategy",
  "country", "Exp_year", "GDD_total", "soil_pH", "SOC", "soil_total_N",
  "total_rainfall", "Irrigation_decision"
)
if (length(args)) {
  master_text <- iconv(readLines(args[[1]], warn = FALSE), from = "GB18030", to = "UTF-8")
  if (anyNA(master_text)) stop("Cannot decode master CSV as GB18030.")
  master <- fread(text = paste(master_text, collapse = "\n"))
  stopifnot(all(input_columns %in% names(master)))
  fwrite(master[, ..input_columns], input_file)
}
d <- fread(input_file, encoding = "UTF-8")
stopifnot(
  nrow(d) == 4613L,
  !anyDuplicated(d$data_ID),
  !anyNA(d[, ..input_columns]),
  all(d$mean_t > 0 & d$mean_c > 0),
  all(d$sd_t >= 0 & d$sd_c >= 0),
  all(d$n_t > 0 & d$n_c > 0),
  max(abs(d$lnRR - log(d$mean_t / d$mean_c))) < 1e-8
)
write_table <- function(x, name) fwrite(x, file.path(out, paste0(name, ".csv")))

# 2. 用当前 SD 重建共享对照协方差；保留原模型含 env_ID 的对照定义。
d[, yi := lnRR]
d[, treatment_var := sd_t^2 / (n_t * mean_t^2)]
d[, control_var := sd_c^2 / (n_c * mean_c^2)]
d[, vi := treatment_var + control_var]
stopifnot(all(is.finite(d$vi) & d$vi > 0))
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
pair_blocks <- blocks[lengths(blocks) > 1L]
pairs <- do.call(cbind, lapply(pair_blocks, combn, m = 2))
stopifnot(all(vapply(blocks, function(index) {
  max(abs(d$control_var[index] - d$control_var[index[1]])) < 1e-12
}, logical(1))))
V <- sparseMatrix(
  i = c(d$model_row, pairs[1, ], pairs[2, ]),
  j = c(d$model_row, pairs[2, ], pairs[1, ]),
  x = c(d$vi, d$control_var[pairs[1, ]], d$control_var[pairs[1, ]]),
  dims = c(nrow(d), nrow(d))
)
stopifnot(isSymmetric(V), max(abs(diag(V) - d$vi)) < 1e-12)

# 3. 保留原版变换参数、参考水平及氮量分组，刷新所有原始协变量。
transformations <- fread(file.path(out, "transformation_parameters.csv"))
for (index in seq_len(nrow(transformations))) {
  specification <- transformations[index]
  value <- d[[specification$source]]
  if (specification$kind == "log1p") value <- log1p(value)
  set(d, j = specification$target, value = (value - specification$center) / specification$divisor)
}
d[, N_c_sq := N_c^2]
d[, rain_c_sq := rain_c^2]
d[, plant_density_c := plant_density - mean(plant_density)]
d[, country := relevel(factor(country), "China")]
d[, irrigation := factor(Irrigation_decision, levels = c("Rainfed", "Irrigation"))]
d[, N_strategy := factor(N_strategy, levels = c("S", "SV", "SVR", "SR", "VR", "V", "R"))]
d[, ref_site_year_ID := paste(reference_ID, site_year_ID, sep = "__")]
d[, effect_id := data_ID]
for (variable in c("reference_ID", "ref_site_year_ID", "effect_id")) {
  set(d, j = variable, value = factor(d[[variable]]))
}
groups <- c("0_100", "100_200", "200_300", "300_400", "400_600")
d[, N_group := cut(
  N_rate,
  breaks = c(0, 100, 200, 300, 400, 600),
  right = FALSE,
  labels = groups
)]
stopifnot(!anyNA(d), all(as.character(d$N_group) == d$N_rate_cluster))
model_formula <- ~ heat_c * N_group + N_c + N_c_sq +
  plant_density_c + N_strategy +
  country + year_c + GDD_c + pH_c + SOC_log_c + TN_log_c +
  rain_c + rain_c_sq + irrigation
X <- model.matrix(model_formula, d)
stopifnot(qr(X)$rank == ncol(X))
write_table(d, "model_input")

# 4. 沿用原版三层随机效应 REML 模型，不调用任何已拟合的历史模型。
message(format(Sys.time(), "%F %T"), " | Fitting Figure2b REML; K = ", nrow(d))
fit <- rma.mv(
  yi = yi,
  V = V,
  mods = model_formula,
  random = list(
    ~ 1 | reference_ID,
    ~ 1 | ref_site_year_ID,
    ~ 1 | effect_id
  ),
  data = d,
  method = "REML",
  sparse = TRUE,
  control = list(optimizer = "nlminb")
)
stopifnot(fit$k == nrow(d), all(is.finite(coef(fit))))
saveRDS(list(fit = fit, formula = model_formula), file.path(out, "model_REML.rds"))
write_table(data.table(
  term = names(coef(fit)),
  estimate = as.numeric(coef(fit)),
  se = fit$se,
  ci_lb = fit$ci.lb,
  ci_ub = fit$ci.ub,
  p = fit$pval
), "coefficients_REML")
covariance <- as.data.table(as.table(vcov(fit)))
setnames(covariance, c("term_1", "term_2", "covariance"))
write_table(covariance, "coefficient_covariance_REML")
write_table(data.table(
  level = c("reference", "reference_site_year", "effect"),
  sigma2 = fit$sigma2
), "variance_components_REML")
interaction_terms <- grep("heat_c:N_group", names(coef(fit)), fixed = TRUE)
stopifnot(length(interaction_terms) == 4L)
joint <- anova(fit, btt = interaction_terms)
write_table(data.table(
  test = "Equal EDH30 slopes across N-rate groups",
  QM = joint$QM,
  df = length(interaction_terms),
  p = joint$QMp
), "interaction_joint_Wald")
write_table(data.table(
  effects = fit$k,
  references = uniqueN(d$reference_ID),
  environments = uniqueN(d$env_ID),
  control_groups = length(blocks),
  shared_control_groups = length(pair_blocks),
  V_nonzero = nnzero(V),
  design_columns = ncol(X),
  design_rank = qr(X)$rank,
  method = fit$method,
  optimizer = "nlminb",
  formula = paste(deparse(model_formula), collapse = " "),
  random_effects = "reference_ID | ref_site_year_ID | effect_id",
  confidence_level = 0.95,
  plant_density_center = mean(d$plant_density),
  country_reference = "China",
  irrigation_reference = "Rainfed",
  N_strategy_reference = "S",
  N_group_reference = "0_100",
  N_group_breaks = "0|100|200|300|400|600; left-closed, right-open",
  logLik_REML = as.numeric(logLik(fit)),
  AIC_REML = AIC(fit),
  BIC_REML = BIC(fit),
  QE = fit$QE,
  QE_p = fit$QEp,
  R_version = as.character(getRversion()),
  metafor_version = as.character(packageVersion("metafor"))
), "model_parameters")
print(data.table(QM = joint$QM, df = 4L, p = joint$QMp))
message(format(Sys.time(), "%F %T"), " | Figure2b fit completed")

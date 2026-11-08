#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(metafor)
})

invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))

message_time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
}

to_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)))))
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "missing"
  trimws(x)
}

# 1. 定位输入与输出；默认只使用本目录文件。
entry <- grep("^--file=", commandArgs(), value = TRUE)
module_dir <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript fit_overall_REML.R [OUTPUT_DIRECTORY]")
table_dir <- if (length(args)) args[[1]] else module_dir
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
data_path <- file.path(module_dir, "comparison_input.csv")
if (!file.exists(data_path)) {
  stop("Input file not found: ", data_path)
}

message_time("Reading input data")
dt <- fread(data_path, encoding = "UTF-8", na.strings = c("", "NA"))
setnames(dt, trimws(names(dt)))
dt[, row_id := .I]

required_cols <- c(
  "reference_ID", "site_year_ID", "Label_of_control_treatment",
  "mean_c", "sd_c", "n_c", "mean_t", "sd_t", "n_t", "lnRR"
)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

numeric_cols <- c("mean_c", "sd_c", "n_c", "mean_t", "sd_t", "n_t", "lnRR")
for (col in numeric_cols) {
  dt[, paste0(col, "_num") := to_num(get(col))]
}

dt[, lnRR_recalculated := log(mean_t_num / mean_c_num)]
dt[, lnRR_abs_diff := abs(lnRR_num - lnRR_recalculated)]
dt[, treatment_var_lnRR := sd_t_num^2 / (n_t_num * mean_t_num^2)]
dt[, control_var_lnRR := sd_c_num^2 / (n_c_num * mean_c_num^2)]
dt[, vi_lnRR := treatment_var_lnRR + control_var_lnRR]
dt[, sei_lnRR := sqrt(vi_lnRR)]

dt[, usable_lnRR := is.finite(lnRR_num) &
     is.finite(treatment_var_lnRR) &
     is.finite(control_var_lnRR) &
     is.finite(vi_lnRR) &
     is.finite(sei_lnRR) &
     vi_lnRR > 0 &
     treatment_var_lnRR >= 0 &
     control_var_lnRR >= 0 &
     mean_c_num > 0 &
     mean_t_num > 0 &
     sd_c_num >= 0 &
     sd_t_num >= 0 &
     n_c_num > 0 &
     n_t_num > 0]

dd <- copy(dt[usable_lnRR == TRUE])
if (nrow(dd) == 0) {
  stop("No rows with finite lnRR and positive vi_lnRR were found.")
}

dd[, reference_ID_chr := safe_chr(reference_ID)]
dd[, site_year_ID_chr := safe_chr(site_year_ID)]
dd[, ref_site_year_ID := paste(reference_ID_chr, site_year_ID_chr, sep = "__")]
dd[, control_id := paste(
  reference_ID_chr,
  site_year_ID_chr,
  safe_chr(Label_of_control_treatment),
  sprintf("%.10g", mean_c_num),
  sprintf("%.10g", sd_c_num),
  sprintf("%.10g", n_c_num),
  sep = "||"
)]

setorder(dd, control_id, row_id)
dd[, model_row := .I]
dd[, control_block := .GRP, by = control_id]
dd[, effect_id := factor(model_row)]
dd[, reference_ID := factor(reference_ID_chr)]
dd[, site_year_ID := factor(site_year_ID_chr)]
dd[, ref_site_year_ID := factor(ref_site_year_ID)]

# 2. 以共享对照构建采样协方差矩阵，保留原分组与排序。
message_time("Building sparse shared-control V matrix")
block_summary <- dd[, .(
  treatment_arms = .N,
  reference_ID = reference_ID_chr[1],
  site_year_ID = site_year_ID_chr[1],
  mean_c = mean_c_num[1],
  sd_c = sd_c_num[1],
  n_c = n_c_num[1],
  control_var_lnRR = control_var_lnRR[1],
  min_vi_lnRR = min(vi_lnRR),
  median_vi_lnRR = stats::median(vi_lnRR),
  max_vi_lnRR = max(vi_lnRR)
), by = .(control_block, control_id)]

shared_ids <- block_summary[treatment_arms > 1, control_id]
pair_list <- dd[control_id %in% shared_ids, {
  if (.N < 2) {
    NULL
  } else {
    pairs <- utils::combn(model_row, 2)
    data.table(
      i = c(pairs[1, ], pairs[2, ]),
      j = c(pairs[2, ], pairs[1, ]),
      x = rep(control_var_lnRR[1], 2 * ncol(pairs)),
      control_id = control_id[1],
      term = "shared_control_covariance"
    )
  }
}, by = control_id]

V_triplet <- rbindlist(
  list(
    dd[, .(i = model_row, j = model_row, x = vi_lnRR, control_id = control_id, term = "variance")],
    pair_list[, .(i, j, x, control_id, term)]
  ),
  use.names = TRUE,
  fill = TRUE
)

V_sparse <- sparseMatrix(
  i = V_triplet$i,
  j = V_triplet$j,
  x = V_triplet$x,
  dims = c(nrow(dd), nrow(dd)),
  symmetric = FALSE
)

# 3. 按原算法拟合总体模型；不执行无关的漏斗图和 Egger 分析。
message_time("Fitting overall hierarchical metafor::rma.mv with shared-control V and REML")
fit <- rma.mv(
  yi = lnRR_num,
  V = V_sparse,
  mods = ~ 1,
  random = list(
    ~ 1 | reference_ID,
    ~ 1 | ref_site_year_ID,
    ~ 1 | effect_id
  ),
  data = dd,
  method = "REML",
  sparse = TRUE,
  control = list(optimizer = "nlminb")
)

est <- as.numeric(coef(fit)[1])
sigma2_sum <- sum(fit$sigma2)
sigma2_prop <- if (sigma2_sum > 0) fit$sigma2 / sigma2_sum else rep(NA_real_, length(fit$sigma2))

overall_summary <- data.table(
  model = "overall_lnRR_hierarchical_shared_control_V",
  package = "metafor",
  function_call = "rma.mv",
  method = fit$method,
  mods = "~ 1",
  random_structure = "~ 1 | reference_ID; ~ 1 | ref_site_year_ID; ~ 1 | effect_id",
  sampling_variance_structure = "sparse shared-control V",
  estimate_lnRR = est,
  se = fit$se[1],
  zval = fit$zval[1],
  pval = fit$pval[1],
  ci_lb_lnRR = fit$ci.lb[1],
  ci_ub_lnRR = fit$ci.ub[1],
  response_ratio = exp(est),
  response_change_percent = (exp(est) - 1) * 100,
  ci_lb_response_change_percent = (exp(fit$ci.lb[1]) - 1) * 100,
  ci_ub_response_change_percent = (exp(fit$ci.ub[1]) - 1) * 100,
  sigma2_reference_ID = fit$sigma2[1],
  sigma2_ref_site_year_ID = fit$sigma2[2],
  sigma2_effect_id = fit$sigma2[3],
  sigma2_reference_ID_percent_of_random_effects = 100 * sigma2_prop[1],
  sigma2_ref_site_year_ID_percent_of_random_effects = 100 * sigma2_prop[2],
  sigma2_effect_id_percent_of_random_effects = 100 * sigma2_prop[3],
  QE = fit$QE,
  QE_df = fit$k - fit$p,
  QE_pval = fit$QEp,
  k = fit$k,
  unique_references = uniqueN(dd$reference_ID),
  unique_site_years = uniqueN(dd$site_year_ID),
  unique_ref_site_years = uniqueN(dd$ref_site_year_ID),
  control_blocks = nrow(block_summary),
  shared_control_blocks = block_summary[treatment_arms > 1, .N],
  logLik = as.numeric(logLik(fit)),
  AIC = AIC(fit),
  BIC = BIC(fit)
)

# 4. 只保存总体拟合结果及截距协方差。
fwrite(overall_summary, file.path(table_dir, "overall_REML_results.csv"))
fwrite(
  data.table(term_i = "intercept", term_j = "intercept", covariance = as.numeric(vcov(fit))),
  file.path(table_dir, "coefficient_covariance.csv")
)
message_time("Saved overall REML results; k = ", fit$k)

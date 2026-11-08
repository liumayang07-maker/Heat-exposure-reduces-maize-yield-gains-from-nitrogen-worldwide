#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))

# 1. 读取当前 4,613 对数据的必要字段快照；不修改主数据。
entry <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript reproduce_SD_imputation.R [OUTPUT_DIRECTORY] [LOOKUP_CSV]")
input_csv <- file.path(root, "current_SD_input.csv")
output_dir <- if (length(args)) args[[1]] else root
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. 沿用原观测 CV 计算和训练记录去重口径。
df <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(df) == 4613L, !anyDuplicated(df$data_ID))
names(df) <- sub("^\xEF\xBB\xBF", "", names(df), useBytes = TRUE)

required_columns <- c(
  "data_ID", "reference_ID", "site_year_ID", "env_ID", "continent",
  "N_rate", "Label of N addition treatment",
  "mean_t", "sd_t", "sd_t_lab", "n_t",
  "mean_c", "sd_c", "sd_c_lab", "n_c",
  "Label_of_control_treatment"
)
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

to_numeric <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
for (cc in c("N_rate", "mean_t", "sd_t", "n_t", "mean_c", "sd_c", "n_c")) {
  df[[cc]] <- to_numeric(df[[cc]])
}

valid_sd_label <- function(x) {
  tolower(trimws(as.character(x))) %in% c("origin", "original", "converted", "strict_converted")
}

has_positive <- function(x) is.finite(x) & x > 0

yield_breaks <- function(x) {
  qs <- as.numeric(stats::quantile(x, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE, type = 7))
  if (!all(is.finite(qs)) || qs[1] >= qs[2]) {
    rng <- range(x, na.rm = TRUE)
    qs <- c(rng[1] + diff(rng) / 3, rng[1] + 2 * diff(rng) / 3)
  }
  qs
}

assign_yield_group <- function(x, breaks) {
  cut(x, breaks = c(-Inf, breaks[1], breaks[2], Inf),
      labels = c("Low", "Medium", "High"), include.lowest = TRUE, right = TRUE)
}

make_training_data <- function(df, status) {
  if (identical(status, "N-fertilized")) {
    keep <- has_positive(df$mean_t) & has_positive(df$sd_t) & is.finite(df$n_t) & valid_sd_label(df$sd_t_lab)
    dat <- df[keep, c("data_ID", "reference_ID", "site_year_ID", "env_ID", "continent",
                      "N_rate", "Label of N addition treatment", "mean_t", "sd_t", "n_t")]
    names(dat)[names(dat) == "mean_t"] <- "Mean"
    names(dat)[names(dat) == "sd_t"] <- "SD"
    names(dat)[names(dat) == "n_t"] <- "n"
    names(dat)[names(dat) == "Label of N addition treatment"] <- "treatment_label"
    dedup_keys <- c("reference_ID", "site_year_ID", "env_ID", "N_rate", "treatment_label", "Mean", "SD", "n")
  } else {
    keep <- has_positive(df$mean_c) & has_positive(df$sd_c) & is.finite(df$n_c) & valid_sd_label(df$sd_c_lab)
    dat <- df[keep, c("data_ID", "reference_ID", "site_year_ID", "env_ID", "continent",
                      "mean_c", "sd_c", "n_c", "Label_of_control_treatment")]
    names(dat)[names(dat) == "mean_c"] <- "Mean"
    names(dat)[names(dat) == "sd_c"] <- "SD"
    names(dat)[names(dat) == "n_c"] <- "n"
    dedup_keys <- c("reference_ID", "site_year_ID", "env_ID", "Label_of_control_treatment", "Mean", "SD", "n")
  }
  dat <- dat[has_positive(dat$Mean) & has_positive(dat$SD), ]
  dat <- dat[!duplicated(dat[, dedup_keys, drop = FALSE]), ]
  dat$Treatment_status <- status
  dat$CV <- dat$SD / dat$Mean
  dat
}

# 3. 计算三层 CV 参数；产量分组阈值按处理状态统一确定。
min_cell_n <- 20L
summarise_cv <- function(observed, groups, level, breaks) {
  key <- interaction(observed[, groups, drop = FALSE], drop = TRUE, sep = "\r")
  do.call(rbind, lapply(split(observed, key), function(z) {
    data.frame(
      level = level,
      Treatment_status = z$Treatment_status[1],
      continent = if ("continent" %in% groups) z$continent[1] else "",
      yield_group = if ("yield_group" %in% groups) z$yield_group[1] else "",
      n_observed_CV = nrow(z),
      median_CV = median(z$CV),
      low_yield_cut = breaks[1],
      high_yield_cut = breaks[2],
      row.names = NULL
    )
  }))
}

lookups <- list()
for (status in c("N-fertilized", "zero-N")) {
  observed <- make_training_data(df, status)
  breaks <- yield_breaks(observed$Mean)
  observed$yield_group <- as.character(assign_yield_group(observed$Mean, breaks))
  lookups[[status]] <- rbind(
    summarise_cv(observed, c("Treatment_status", "continent", "yield_group"),
                 "continent_yield_group", breaks),
    summarise_cv(observed, c("Treatment_status", "continent"), "continent", breaks),
    summarise_cv(observed, "Treatment_status", "global", breaks)
  )
  cat(status, ": observed CV n =", nrow(observed),
      "; yield tertiles =", sprintf("%.12g", breaks), "\n")
}
lookup <- do.call(rbind, lookups)

# 4. 保留观测 SD；其余 SD 按分组、洲、全局 CV 依次重算。
for (arm in c("t", "c")) {
  status <- if (arm == "t") "N-fertilized" else "zero-N"
  sd_col <- paste0("sd_", arm)
  mean_col <- paste0("mean_", arm)
  label_col <- paste0(sd_col, "_lab")
  observed_origin <- valid_sd_label(df[[label_col]]) & has_positive(df[[sd_col]])
  pars <- lookups[[status]]
  yield_group <- as.character(assign_yield_group(
    df[[mean_col]], c(pars$low_yield_cut[1], pars$high_yield_cut[1])
  ))
  source <- rep("observed_origin", nrow(df))
  cv_used <- rep(NA_real_, nrow(df))
  lookup_n <- rep(NA_integer_, nrow(df))
  for (i in which(!observed_origin)) {
    stopifnot(is.finite(df[[mean_col]][i]), df[[mean_col]][i] > 0)
    selected <- pars[pars$level == "continent_yield_group" &
                     pars$continent == df$continent[i] &
                     pars$yield_group == yield_group[i] &
                     pars$n_observed_CV >= min_cell_n, ]
    if (nrow(selected) == 0L) {
      selected <- pars[pars$level == "continent" &
                       pars$continent == df$continent[i] &
                       pars$n_observed_CV >= min_cell_n, ]
    }
    if (nrow(selected) == 0L) selected <- pars[pars$level == "global", ]
    stopifnot(nrow(selected) == 1L, is.finite(selected$median_CV))
    cv_used[i] <- selected$median_CV
    lookup_n[i] <- selected$n_observed_CV
    source[i] <- paste0(selected$level, "_median_CV")
    df[[sd_col]][i] <- df[[mean_col]][i] * cv_used[i]
    df[[label_col]][i] <- "filled"
  }
  df[[paste0(sd_col, "_fill_source")]] <- source
  df[[paste0(sd_col, "_yield_group")]] <- yield_group
  df[[paste0(sd_col, "_median_CV_used")]] <- cv_used
  df[[paste0(sd_col, "_lookup_n_observed_CV")]] <- lookup_n
  stopifnot(all(is.finite(df[[sd_col]])), all(df[[sd_col]] > 0))
  cat(status, ": observed retained =", sum(observed_origin),
      "; SD recomputed =", sum(!observed_origin), "\n")
}

# 5. 仅输出连接主表所需的 SD、来源标签及回退来源。
# 补充表脚本需要 CV 参数时，显式指定临时 CSV 路径。
if (length(args) == 2L) {
  write.csv(lookup, args[[2]], row.names = FALSE, fileEncoding = "UTF-8")
}
output_columns <- c(
  "data_ID",
  "sd_t", "sd_t_lab", "sd_t_fill_source",
  "sd_c", "sd_c_lab", "sd_c_fill_source"
)
write.csv(
  df[, output_columns, drop = FALSE],
  file.path(output_dir, "meta_data_v2_SD_reimputed_4613.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
cat("Saved current 4613-record minimal SD output.\n")

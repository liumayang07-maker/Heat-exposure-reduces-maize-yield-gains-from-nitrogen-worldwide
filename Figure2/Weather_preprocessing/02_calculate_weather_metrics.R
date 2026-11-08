#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
})
setDTthreads(2L)

# 1. 只读取本目录的最简输入；小时气象文件位于外部 --raw-dir。
entry <- grep("^--file=", commandArgs(), value = TRUE)
module <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || length(args) %% 2L != 0L) stop("Usage: Rscript 02_calculate_weather_metrics.R --raw-dir PATH [--output-dir PATH]")
options_cli <- setNames(args[seq(2L, length(args), by = 2L)], args[seq(1L, length(args), by = 2L)])
if (!"--raw-dir" %in% names(options_cli) || any(!names(options_cli) %in% c("--raw-dir", "--output-dir"))) {
  stop("Supply --raw-dir PATH; optional --output-dir PATH.")
}
raw_dir <- options_cli[["--raw-dir"]]
out <- if ("--output-dir" %in% names(options_cli)) options_cli[["--output-dir"]] else module
if (!dir.exists(raw_dir)) stop("Missing raw directory: ", raw_dir)
tasks <- fread(file.path(module, "download_tasks.csv"))
environments <- fread(file.path(module, "environment_input.csv"))
comparisons <- fread(file.path(module, "R1_input.csv"), na.strings = c("", "NA"))
for (column in c("window_start", "window_end", "request_start_utc", "request_end_utc")) {
  tasks[, (column) := as.IDate(get(column))]
}
environments[, `:=`(sowing_date = as.IDate(sowing_date), harvest_date = as.IDate(harvest_date))]
comparisons[, R1_date := as.IDate(R1_date)]
stopifnot(
  nrow(tasks) == 775L,
  !anyDuplicated(tasks$task_id),
  !anyDuplicated(tasks$output_filename),
  nrow(environments) == 856L,
  !anyDuplicated(environments$env_ID),
  setequal(environments$task_id, tasks$task_id),
  nrow(comparisons) == 4613L,
  !anyDuplicated(comparisons$data_ID),
  setequal(comparisons$env_ID, environments$env_ID),
  !anyNA(environments[, .(env_ID, task_id, sowing_date, harvest_date, timezone_name)]),
  all(environments$harvest_date >= environments$sowing_date)
)

# 2. 先核验 UTC 连续小时，再按冻结的当地时区截取含首尾日的物理小时窗口。
physical_hours <- function(start, end, tz) {
  as.integer(round(as.numeric(difftime(
    as.POSIXct(paste(end + 1L, "00:00:00"), tz = tz),
    as.POSIXct(paste(start, "00:00:00"), tz = tz),
    units = "hours"
  ))))
}
period_metrics <- function(raw, local_date, start, end, tz) {
  part <- raw[local_date >= start & local_date <= end]
  expected <- physical_hours(start, end, tz)
  if (nrow(part) != expected || expected <= 0L) stop("Incomplete local-time window: ", start, " / ", end)
  # 原版对小于 -0.001 mm 的小时降水报错；仅将浮点级负值裁为零。
  if (any(part$tp_mm < -0.001)) stop("Material negative hourly precipitation.")
  list(
    part = part,
    GDD = sum(pmin(pmax(part$t2m_c, 10), 30) - 10) / 24,
    EDH30 = sum(pmax(part$t2m_c - 30, 0)),
    precipitation = sum(pmax(part$tp_mm, 0))
  )
}
annual_results <- list()
season_results <- list()
profile_results <- list()
r1_results <- list()
index <- 0L
for (task_number in seq_len(nrow(tasks))) {
  task <- tasks[task_number]
  path <- file.path(raw_dir, task$output_filename)
  if (!file.exists(path)) stop("Missing raw file: ", path)
  raw <- fread(path, select = c("valid_time", "latitude", "longitude", "t2m", "tp"),
               colClasses = list(character = "valid_time"))
  raw[, time_utc := suppressWarnings(ymd_hms(valid_time, tz = "UTC", quiet = TRUE))]
  if (anyNA(raw$time_utc) || anyDuplicated(raw$time_utc)) stop("Invalid UTC timestamps: ", path)
  setorder(raw, time_utc)
  expected_start <- as.POSIXct(task$request_start_utc, tz = "UTC")
  expected_end <- as.POSIXct(task$request_end_utc + 1L, tz = "UTC") - 3600
  expected_n <- as.integer(difftime(expected_end, expected_start, units = "hours")) + 1L
  if (nrow(raw) != expected_n || raw$time_utc[1] != expected_start || tail(raw$time_utc, 1) != expected_end ||
      any(diff(as.numeric(raw$time_utc)) != 3600)) stop("Incomplete UTC hours: ", path)
  if (any(!is.finite(as.matrix(raw[, .(latitude, longitude, t2m, tp)])))) stop("Non-finite weather or coordinates: ", path)
  normalized_lon <- ((raw$longitude + 180) %% 360) - 180
  if (max(abs(raw$latitude - task$download_gp_lat)) > 1e-6 ||
      max(abs(normalized_lon - task$download_gp_lon)) > 1e-6) stop("Incorrect downloaded grid: ", path)
  raw[, `:=`(t2m_c = t2m - 273.15, tp_mm = tp * 1000)]
  rows <- environments[task_id == task$task_id]
  for (i in seq_len(nrow(rows))) {
    row <- rows[i]
    index <- index + 1L
    local_date <- as.IDate(with_tz(raw$time_utc, row$timezone_name), tz = row$timezone_name)
    annual <- period_metrics(raw, local_date, task$window_start, task$window_end, row$timezone_name)
    season <- period_metrics(raw, local_date, row$sowing_date, row$harvest_date, row$timezone_name)
    annual_results[[index]] <- data.table(
      env_ID = row$env_ID, task_id = task$task_id,
      window_start = task$window_start, window_end = task$window_end,
      annual_DD10_30_C_d = annual$GDD,
      annual_EDH30_C_h = annual$EDH30,
      annual_precip_mm = annual$precipitation
    )
    season_results[[index]] <- data.table(
      env_ID = row$env_ID, task_id = task$task_id,
      sowing_date = row$sowing_date, harvest_date = row$harvest_date,
      growing_season_GDD10_30_C_d = season$GDD,
      growing_season_EDH30_C_h = season$EDH30,
      growing_season_precip_mm = season$precipitation
    )

    # 3. 沿用原版逐小时累计 GDD 分箱；每 bin 强度为 EDH30 总量 / 小时数。
    part <- season$part
    part[, GDD_h := (pmin(pmax(t2m_c, 10), 30) - 10) / 24]
    total_gdd_hourly <- sum(part$GDD_h)
    if (!is.finite(total_gdd_hourly) || total_gdd_hourly <= 0) stop("Non-positive growing-season GDD: ", row$env_ID)
    part[, bin := pmin(100L, pmax(1L, as.integer(floor(cumsum(GDD_h) / total_gdd_hourly * 100) + 1L)))]
    part[, EDH30_h := pmax(t2m_c - 30, 0)]
    binned <- part[, .(bin_hours = .N, EDH30_sum = sum(EDH30_h)), by = bin][
      data.table(bin = 1:100), on = "bin"
    ]
    if (anyNA(binned$bin_hours) || any(binned$bin_hours <= 0L)) stop("Empty thermal-progress bin: ", row$env_ID)
    binned[, `:=`(env_ID = row$env_ID, EDH30_per_hour = EDH30_sum / bin_hours)]
    if (abs(sum(binned$EDH30_sum) - season$EDH30) > 1e-7) stop("EDH30 bin totals do not match the growing season.")
    profile_results[[index]] <- binned[, .(env_ID, bin, bin_hours, EDH30_per_hour)]

    # 4. 对每个实际 R1 日期单独计算播种至 R1 当日（含）的积温，不提前合并日期。
    cases <- unique(comparisons[env_ID == row$env_ID, .(env_ID, R1_date)])
    season_dates <- local_date[local_date >= row$sowing_date & local_date <= row$harvest_date]
    cases[, `:=`(
      authoritative_growing_season_GDD10_30_C_d = season$GDD,
      pre_silking_total_GDD = NA_real_,
      post_silking_total_GDD = NA_real_,
      status = "R1_missing"
    )]
    for (j in seq_len(nrow(cases))) {
      r1_date <- cases$R1_date[j]
      if (is.na(r1_date)) next
      if (r1_date < row$sowing_date || r1_date > row$harvest_date) {
        cases[j, status := "R1_outside_growing_period"]
      } else {
        pre <- sum(part$GDD_h[season_dates <= r1_date])
        post <- sum(part$GDD_h[season_dates > r1_date])
        if (abs(pre + post - season$GDD) > 1e-8) stop("R1 GDD decomposition mismatch.")
        cases[j, `:=`(pre_silking_total_GDD = pre, post_silking_total_GDD = post, status = "complete")]
      }
    }
    r1_results[[index]] <- cases
  }
  if (task_number %% 50L == 0L || task_number == nrow(tasks)) message("Processed tasks: ", task_number, "/", nrow(tasks))
}

# 5. 全部窗口通过检查后，才写入四份必要结果；不改主表或下游冻结输入。
annual <- rbindlist(annual_results)
season <- rbindlist(season_results)
profiles <- rbindlist(profile_results)
r1_cases <- rbindlist(r1_results)
r1 <- merge(comparisons, r1_cases, by = c("env_ID", "R1_date"), all.x = TRUE, sort = FALSE)
setorder(annual, env_ID)
setorder(season, env_ID)
setorder(profiles, env_ID, bin)
setorder(r1, data_ID)
stopifnot(
  nrow(annual) == 856L,
  nrow(season) == 856L,
  nrow(profiles) == 85600L,
  nrow(r1) == 4613L,
  !anyNA(r1$status),
  !anyDuplicated(r1$data_ID),
  all(profiles$EDH30_per_hour >= 0)
)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(annual, file.path(out, "annual_metrics.csv"))
fwrite(season, file.path(out, "growing_season_metrics.csv"))
fwrite(profiles, file.path(out, "profiles_input.csv"))
fwrite(r1, file.path(out, "pre_silking_GDD_by_data_ID.csv"))
message("Completed: 775 weather tasks; 856 environments; 85,600 bins; 4,613 comparisons.")
print(r1[, .N, by = status])

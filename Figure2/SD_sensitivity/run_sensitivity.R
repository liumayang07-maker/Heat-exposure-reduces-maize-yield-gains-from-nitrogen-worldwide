#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(metafor)
})

# 1. 原生 R 入口：复用原分析表达式，所有写入均重定向到独立情景目录。
# 共享对照身份沿用原 SD 定义；数值方差和协方差仍由情景 SD 完整重算。
compact_model <- function(x) {
  # 模型公式只引用数据列；移除局部公式环境，避免序列化整套运行状态。
  if (is.list(x)) {
    metadata <- attributes(x)
    x <- lapply(unclass(x), compact_model)
    attributes(x) <- metadata
  }
  if (!is.null(attr(x, ".Environment", exact = TRUE))) attr(x, ".Environment") <- baseenv()
  x
}

freeze_control_identity <- function(expression) {
  if (!is.call(expression)) return(expression)
  if (identical(expression[[1]], as.name("sprintf")) &&
      length(expression) == 3L && identical(expression[[2]], "%.10g") && is.symbol(expression[[3]]) &&
      as.character(expression[[3]]) %in% c("sd_c", "sd_c_num")) {
    expression[[3]] <- as.name("sd_c_identity")
    return(expression)
  }
  for (i in seq_along(expression)[-1]) {
    if (is.call(expression[[i]])) expression[[i]] <- freeze_control_identity(expression[[i]])
  }
  expression
}

evaluate_original <- function(path, environment, stop_at = NULL, freeze = FALSE) {
  expressions <- parse(path, encoding = "UTF-8")
  if (!is.null(stop_at)) {
    hits <- which(vapply(expressions, function(x) {
      grepl(stop_at, paste(deparse(x), collapse = " "), fixed = TRUE)
    }, logical(1)))
    stopifnot(length(hits) == 1L)
    expressions <- expressions[seq_len(hits)]
  }
  for (expression in expressions) {
    if (freeze) expression <- freeze_control_identity(expression)
    eval(expression, envir = environment)
  }
}

run_sensitivity_work <- function(figure, root, output, selected = NULL, resume = FALSE) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  model_names <- if (figure == "Figure2") paste0("Figure2", letters[1:4]) else "Figure3b"
  scenarios <- c(Baseline = NA_real_, CV05 = 0.05, CV10 = 0.10, CV20 = 0.20)
  if (!is.null(selected)) scenarios <- scenarios[names(scenarios) %in% selected]
  stopifnot(length(scenarios) > 0L)
  sd_input <- fread(file.path(root, "Figure2/Supplementary_Table_SD_imputation/current_SD_input.csv"))
  # 归档插补结果提供原 SD 与来源标签，不再依赖 manuscripts 外的主表。
  master <- fread(file.path(root, "Figure2/Supplementary_Table_SD_imputation/meta_data_v2_SD_reimputed_4613.csv"))
  stopifnot(nrow(master) == 4613L, !anyDuplicated(master$data_ID))
  statuses <- list()
  validations <- list()
  if (resume) {
    status_path <- file.path(output, "run_status.csv")
    validation_path <- file.path(output, "validation.csv")
    if (file.exists(status_path)) statuses <- list(fread(status_path))
    if (file.exists(validation_path)) validations <- list(fread(validation_path))
  }
  record <- function(model, scenario, check, difference = 0, passed = TRUE) {
    validations[[length(validations) + 1L]] <<- data.table(
      model = model, scenario = scenario, check = check,
      max_abs_difference = difference, passed = passed
    )
    fwrite(rbindlist(validations), file.path(output, "validation.csv"))
    if (!passed) stop("Validation failed: ", check, "; difference = ", difference)
  }
  prediction_files <- c(
    Figure2b = "N_group_EDH30_percentile_curves_REML.csv",
    Figure2c = "adjusted_N_response_curves_REML.csv",
    Figure2d = "fixed_N_percentile_curves_REML.csv",
    Figure3b = "regional_N_EDH30_predictions_REML.csv"
  )
  source_files <- unlist(lapply(model_names, function(model) {
    folder <- file.path(root, figure, model)
    list.files(folder, pattern = "\\.(R|csv|rds|pdf|png|tif|svg)$", full.names = TRUE)
  }))
  source_files <- c(source_files, file.path(
    root,
    "Figure2/Supplementary_Table_SD_imputation",
    c("current_SD_input.csv", "meta_data_v2_SD_reimputed_4613.csv")
  ))
  original_hashes <- tools::md5sum(source_files)
  if (resume) {
    previous_hashes <- fread(file.path(output, "original_file_checksums.csv"))
    stopifnot(identical(previous_hashes$md5, unname(tools::md5sum(previous_hashes$path))))
  }
  fwrite(data.table(path = names(original_hashes), md5 = unname(original_hashes)),
         file.path(output, "original_file_checksums.csv"))
  writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))

  # 2. 各模型严格使用其归档样本及顺序，不从主表重新筛选分析样本。
  for (model in model_names) {
    source_dir <- file.path(root, figure, model)
    original <- fread(file.path(source_dir, "comparison_input.csv"))
    index <- match(original$data_ID, master$data_ID)
    observed_index <- match(original$data_ID, sd_input$data_ID)
    stopifnot(!anyNA(index), !anyNA(observed_index), !anyDuplicated(original$data_ID))
    for (scenario in names(scenarios)) {
      if (resume && length(statuses)) {
        finished <- rbindlist(statuses)
        if (any(finished$model == model & finished$scenario == scenario & finished$status == "success")) next
        matching <- finished$model == model & finished$scenario == scenario
        previous <- finished[matching]
        if (nrow(previous)) {
          history <- file.path(output, "retry_history.csv")
          fwrite(previous, history, append = file.exists(history))
          statuses <- list(finished[!matching])
        }
      }
      started <- Sys.time()
      warnings <- character()
      target <- file.path(output, "results", model, scenario)
      dir.create(target, recursive = TRUE, showWarnings = FALSE)
      if (!resume && file.exists(file.path(target, "model_REML.rds"))) {
        stop("Existing result found; refusing to overwrite: ", target)
      }
      message(format(Sys.time()), " | ", model, " | ", scenario)
      status <- tryCatch(withCallingHandlers({
        d <- copy(original)
        for (arm in c("t", "c")) {
          column <- paste0("sd_", arm)
          mean_column <- paste0("mean_", arm)
          flags <- master[[paste0(column, "_lab")]][index]
          stopifnot(all(flags %in% c("origin", "original", "converted", "strict_converted", "filled")))
          difference <- max(abs(d[[column]] - master[[column]][index]))
          record(model, scenario, paste0(column, "_matches_current_master"), difference, difference < 1e-8)
          missing <- flags == "filled"
          input_flags <- sd_input[[paste0(column, "_lab")]][observed_index]
          input_sd <- sd_input[[column]][observed_index]
          observed <- input_flags %in% c("origin", "original", "converted", "strict_converted") &
                      is.finite(input_sd) & input_sd > 0
          record(model, scenario, paste0(column, "_provenance_mask_verified"), passed = identical(missing, !observed))
          if (scenario != "Baseline") set(d, which(missing), column, d[[mean_column]][missing] * scenarios[[scenario]])
          record(model, scenario, paste0(column, "_observed_preserved"),
                 max(abs(d[[column]][!missing] - original[[column]][!missing])),
                 identical(d[[column]][!missing], original[[column]][!missing]))
        }
        untouched <- setdiff(names(original), c("sd_t", "sd_c"))
        record(model, scenario, "non_SD_inputs_identical", passed = identical(d[, ..untouched], original[, ..untouched]))
        d[, sd_c_identity := original$sd_c]
        fwrite(d, file.path(target, "comparison_input.csv"))
        parameter_file <- file.path(source_dir, "transformation_parameters.csv")
        if (file.exists(parameter_file)) file.copy(parameter_file, target, overwrite = FALSE)
        environment <- new.env(parent = globalenv())
        command_environment <- new.env(parent = baseenv())
        command_environment$input_path <- file.path(target, "comparison_input.csv")
        environment$commandArgs <- eval(quote(function(trailingOnly = FALSE) {
          if (trailingOnly) character() else paste0("--file=", input_path)
        }), command_environment)
        environment$saveRDS <- function(object, file, ...) {
          base::saveRDS(compact_model(object), file, ...)
        }
        script <- if (model == "Figure2a") "fit_overall_REML.R" else paste0("fit_", model, ".R")
        evaluate_original(file.path(source_dir, script), environment, freeze = TRUE)
        fit <- environment$fit
        fit_data <- if (model == "Figure2a") environment$dd else environment$d
        V <- if (model == "Figure2a") environment$V_sparse else environment$V
        record(model, scenario, "sample_size_preserved", passed = fit$k == nrow(original))
        expected_vi <- fit_data$sd_t^2 / (fit_data$n_t * fit_data$mean_t^2) +
                       fit_data$sd_c^2 / (fit_data$n_c * fit_data$mean_c^2)
        record(model, scenario, "V_recomputed_from_SD", max(abs(diag(V) - expected_vi)),
               max(abs(diag(V) - expected_vi)) < 1e-12)
        if (model == "Figure2a") saveRDS(compact_model(list(fit = fit, formula = ~1)), file.path(target, "model_REML.rds"))
        # 保存稀疏 V 及逐行身份，便于独立核验共享对照结构与更新后的协方差。
        saveRDS(V, file.path(target, "sampling_covariance.rds"))
        fwrite(fit_data[, .(data_ID, control_id)], file.path(target, "control_membership.csv"))
        reference_fit_file <- file.path(source_dir, "model_REML.rds")
        if (file.exists(reference_fit_file)) {
          reference <- readRDS(reference_fit_file)$fit
          record(model, scenario, "design_matrix_identical", max(abs(fit$X - reference$X)),
                 identical(dimnames(fit$X), dimnames(reference$X)) && max(abs(fit$X - reference$X)) < 1e-12)
          record(model, scenario, "response_identical", max(abs(fit$yi - reference$yi)),
                 max(abs(fit$yi - reference$yi)) < 1e-12)
          if (scenario == "Baseline") {
            difference <- max(abs(coef(fit) - coef(reference)))
            record(model, scenario, "baseline_coefficients_reproduced", difference, difference < 1e-6)
            difference <- max(abs(vcov(fit) - vcov(reference)))
            record(model, scenario, "baseline_coefficient_covariance_reproduced", difference, difference < 1e-6)
          }
        } else if (model == "Figure2a" && scenario == "Baseline") {
          reference <- fread(file.path(source_dir, "overall_REML_results.csv"))
          difference <- max(abs(c(coef(fit)[1], fit$ci.lb[1], fit$ci.ub[1]) -
                               c(reference$estimate_lnRR, reference$ci_lb_lnRR, reference$ci_ub_lnRR)))
          record(model, scenario, "baseline_overall_reproduced", difference, difference < 1e-6)
        }
        baseline_dir <- file.path(output, "results", model, "Baseline")
        if (scenario != "Baseline") {
          membership <- fread(file.path(baseline_dir, "control_membership.csv"))
          record(model, scenario, "control_membership_identical", passed =
                   identical(membership, fread(file.path(target, "control_membership.csv"))))
          base_fit <- readRDS(file.path(baseline_dir, "model_REML.rds"))$fit
          record(model, scenario, "baseline_design_identical", max(abs(fit$X - base_fit$X)),
                 max(abs(fit$X - base_fit$X)) < 1e-12)
        }
        if (model == "Figure3b") {
          rows <- list()
          for (reg in environment$region_levels) {
            nd <- copy(environment$template)
            nd[, region := factor(reg, levels = levels(environment$d$region))]
            nd[, `:=`(heat_wr = 0, N_c = 0, N_c_sq = 0)]
            x0 <- environment$matrix_for(environment$form, nd)
            nd[, heat_wr := 1]
            x1 <- environment$matrix_for(environment$form, nd)
            L <- rbind(x1 - x0, environment$interaction_contrasts(reg))
            estimate <- drop(L %*% coef(fit))
            se <- sqrt(diag(L %*% vcov(fit) %*% t(L)))
            rows[[reg]] <- data.table(
              region = reg, region_abbr = unname(environment$region_abbr[reg]),
              term = c("heat_wr at centered N", "heat_wr:N_c", "heat_wr:N_c_sq"),
              estimate = estimate, se = se,
              ci_lb = estimate - qnorm(.975) * se,
              ci_ub = estimate + qnorm(.975) * se,
              p = 2 * pnorm(-abs(estimate / se))
            )
          }
          fwrite(rbindlist(rows), file.path(target, "regional_heat_coefficients.csv"))
        }
        # 3. 执行原绘图脚本的预测与检验部分，止于绘图之前，原网格不作抽样。
        if (model != "Figure2a") {
          stop_at <- if (model == "Figure2d") "tab(fixed_N_tests," else prediction_files[[model]]
          evaluate_original(file.path(source_dir, paste0("plot_", model, ".R")), environment, stop_at)
          predicted <- fread(file.path(target, prediction_files[[model]]))
          reference <- fread(file.path(source_dir, prediction_files[[model]]))
          grid_columns <- intersect(c("N_group", "Env_Type", "region", "region_abbr", "N_scenario",
                                      "N_rate", "N_control_kg_ha", "percentile", "log1p_EDH30", "EDH30_total"), names(predicted))
          record(model, scenario, "prediction_grid_identical", passed =
                   isTRUE(all.equal(predicted[, ..grid_columns], reference[, ..grid_columns], tolerance = 1e-12)))
          if (scenario == "Baseline") {
            numeric_columns <- names(predicted)[vapply(predicted, is.numeric, logical(1))]
            difference <- max(abs(as.matrix(predicted[, ..numeric_columns]) - as.matrix(reference[, ..numeric_columns])))
            record(model, scenario, "baseline_predictions_reproduced", difference, difference < 1e-5)
          }
        }
        "success"
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        # 记录但不隐藏数值或收敛告警，不自动换优化器或剔除样本。
      }), error = function(e) paste0("ERROR: ", conditionMessage(e)))
      statuses[[length(statuses) + 1L]] <- data.table(
        model = model, scenario = scenario,
        status = if (status == "success") "success" else "failed",
        error = if (status == "success") "" else status,
        warnings = paste(unique(warnings), collapse = " | "),
        elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
      )
      fwrite(rbindlist(statuses), file.path(output, "run_status.csv"))
      message(model, " / ", scenario, ": ", status)
      gc()
    }
  }
  record(figure, "All", "original_files_unchanged", passed = identical(original_hashes, tools::md5sum(source_files)))
  if (any(rbindlist(statuses)$status != "success")) stop("Some scenarios failed; see run_status.csv")
}

run_sensitivity <- function(figure, root, selected = NULL, resume = FALSE) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  stopifnot(figure %in% c("Figure2", "Figure3"))
  output <- file.path(root, figure, "SD_sensitivity")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  state_file <- file.path(output, "workspace_state.csv")
  temporary_root <- normalizePath(Sys.getenv("TEMP", unset = tempdir()), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(output, "model_status.csv"))) {
    stop("Completed core CSVs already exist; refusing to overwrite: ", output)
  }

  # 4. 临时结果沿用原目录结构；失败保留工作目录，成功后只发布核心 CSV。
  if (resume) {
    if (!file.exists(state_file)) stop("No temporary workspace to resume: ", state_file)
    state <- fread(state_file)
    stopifnot(nrow(state) == 1L, state$figure == figure, state$status != "success")
    work <- normalizePath(state$workspace, winslash = "/", mustWork = TRUE)
    stopifnot(identical(dirname(work), temporary_root), grepl("^SD_sensitivity_", basename(work)))
  } else {
    if (file.exists(state_file) || file.exists(file.path(output, "run_status.csv"))) {
      stop("Existing run found; refusing to overwrite. Resume requires workspace_state.csv.")
    }
    work <- tempfile(pattern = "SD_sensitivity_", tmpdir = temporary_root)
    dir.create(work)
  }
  state <- data.table(figure = figure, workspace = work, status = "running")
  fwrite(state, state_file)
  completed <- FALSE
  on.exit({
    if (!completed) {
      state[, status := "incomplete"]
      fwrite(state, state_file)
      message("Temporary results retained for --resume: ", work)
    }
  }, add = TRUE)
  run_sensitivity_work(figure, root, work, selected = selected, resume = resume)

  # 汇总拟合结果；仓库仅保留分析代码及输入输出。
  retained <- file.path(work, "retained")
  dir.create(retained, showWarnings = FALSE)
  python <- Sys.which("py")
  python_args <- "-3"
  if (!nzchar(python)) {
    python <- Sys.which("python")
    python_args <- character()
  }
  if (!nzchar(python)) stop("Python is required for the core CSV summarizer.")
  result <- system2(
    command = python,
    args = c(
      python_args,
      shQuote(file.path(root, "Figure2/SD_sensitivity/summarize_sensitivity.py")),
      "--figure", figure,
      "--input-root", shQuote(work),
      "--output-dir", shQuote(retained)
    )
  )
  if (result != 0L) stop("Core CSV summarizer failed; temporary results retained.")
  files <- list.files(retained, full.names = TRUE)
  stopifnot(length(files) > 0L, all(grepl("\\.csv$", files)),
            "model_status.csv" %in% basename(files))
  destinations <- file.path(output, basename(files))
  existing <- file.exists(destinations)
  # 允许继续一次中断的发布，但不覆盖内容不同的已有结果。
  if (any(existing)) {
    stopifnot(identical(unname(tools::md5sum(files[existing])),
                        unname(tools::md5sum(destinations[existing]))))
  }
  # 完成标记最后发布，避免中途复制失败被误判为已完成。
  publish_order <- order(basename(files) == "model_status.csv")
  pending <- publish_order[!existing[publish_order]]
  for (i in pending) {
    stopifnot(file.copy(files[i], destinations[i], overwrite = FALSE),
              identical(unname(tools::md5sum(files[i])), unname(tools::md5sum(destinations[i]))))
  }
  stopifnot(identical(unname(tools::md5sum(files)), unname(tools::md5sum(destinations))))
  # work 仅能是本函数在系统临时目录下创建并验证的专用目录。
  stopifnot(identical(dirname(normalizePath(work, winslash = "/")), temporary_root),
            grepl("^SD_sensitivity_", basename(work)))
  stopifnot(unlink(work, recursive = TRUE) == 0L)
  stopifnot(unlink(state_file) == 0L)
  completed <- TRUE
  message(figure, ": core CSVs retained in ", output)
}

if (sys.nframe() == 0L) {
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  root <- normalizePath(file.path(dirname(entry), "../.."), winslash = "/")
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(all(args %in% "--resume"))
  run_sensitivity("Figure2", root, resume = "--resume" %in% args)
}

#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({
  library(data.table)
  library(mclust)
  library(cluster)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(openxlsx)
})

# 1. 仅使用本目录的核心输入；K=3 为正式结果，其他 K 仅用于诊断。
entry <- grep("^--file=", commandArgs(), value = TRUE)
out <- dirname(normalizePath(sub("^--file=", "", entry), winslash = "/"))
bin_csv <- file.path(out, "profiles_input.csv")
cv_repeats <- 20L
stability_repeats <- 100L
seed <- 20260910L
k_values <- 2:8
cv_folds <- 5L
full_nstart <- 200L
resample_nstart <- 25L
stability_fraction <- 0.80

required_cols <- c("env_ID", "bin", "bin_hours", "EDH30_per_hour")
bin_dt <- fread(bin_csv)
missing_cols <- setdiff(required_cols, names(bin_dt))
if (length(missing_cols) > 0L) stop("Binned input is missing: ", paste(missing_cols, collapse = ", "))
bin_dt[, bin := as.integer(bin)]
bin_dt[, bin_hours := as.numeric(bin_hours)]
bin_dt[, value_raw := as.numeric(EDH30_per_hour)]
if (anyNA(bin_dt$value_raw) || any(!is.finite(bin_dt$value_raw)) || any(bin_dt$value_raw < 0)) {
  stop("EDH30_per_hour must be finite and non-negative.")
}
if (anyNA(bin_dt$bin_hours) || any(!is.finite(bin_dt$bin_hours)) || any(bin_dt$bin_hours <= 0)) {
  stop("bin_hours must be finite and positive.")
}
if (anyDuplicated(bin_dt[, .(env_ID, bin)])) stop("Duplicate env_ID-bin rows detected.")
env_qc <- bin_dt[, .(
  n_rows = .N,
  n_bins = uniqueN(bin),
  min_bin = min(bin),
  max_bin = max(bin),
  nonzero_bins = sum(value_raw > 0),
  zero_fraction = mean(value_raw == 0),
  cumulative_EDH30 = sum(value_raw * bin_hours),
  profile_sum = sum(value_raw)
), by = env_ID]
if (any(env_qc$n_rows != 100L | env_qc$n_bins != 100L |
        env_qc$min_bin != 1L | env_qc$max_bin != 100L)) {
  stop("Every environment must have exactly one row for bins 1 through 100.")
}
stopifnot(nrow(env_qc) == 856L)
message(format(Sys.time(), "%F %T"), " | Fitting positive-exposure Gaussian mixture")
positive <- env_qc[cumulative_EDH30 > 0]
zero_envs <- env_qc[cumulative_EDH30 == 0, env_ID]
if (length(zero_envs) == 0L || nrow(positive) < max(k_values)) stop("Invalid zero/positive split.")
intensity_x <- log1p(positive$cumulative_EDH30)
set.seed(seed)
intensity_fit <- Mclust(
  intensity_x, G = 1:5, modelNames = c("E", "V"), verbose = FALSE
)
if (intensity_fit$G < 2L) stop("The selected intensity model has no separable Low component.")

component_order <- order(drop(intensity_fit$parameters$mean))
component_rank <- setNames(seq_along(component_order), component_order)
component_names <- if (intensity_fit$G == 2L) {
  c("Low", "High")
} else if (intensity_fit$G == 3L) {
  c("Low", "Moderate", "High")
} else {
  paste("Intensity", seq_len(intensity_fit$G))
}
map_rank <- as.integer(component_rank[as.character(intensity_fit$classification)])
map_label <- component_names[map_rank]
posterior_ordered <- intensity_fit$z[, component_order, drop = FALSE]
colnames(posterior_ordered) <- paste0("posterior_", tolower(component_names))
# 2. 保留原版 MAP 分型、逐环境标准化、K-means 和重采样算法。
positive_intensity <- cbind(
  positive,
  data.table(
    log1p_cumulative_EDH30 = intensity_x,
    MAP_component = intensity_fit$classification,
    MAP_intensity_rank = map_rank,
    MAP_intensity_class = map_label,
    MAP_probability = apply(intensity_fit$z, 1, max),
    classification_uncertainty = intensity_fit$uncertainty
  ),
  as.data.table(posterior_ordered)
)
if (anyDuplicated(names(positive_intensity))) {
  duplicate_names <- unique(names(positive_intensity)[duplicated(names(positive_intensity))])
  stop("Duplicate columns in positive_intensity: ", paste(duplicate_names, collapse = ", "))
}
zero_intensity <- copy(env_qc[env_ID %in% zero_envs])
zero_intensity[, `:=`(
  log1p_cumulative_EDH30 = 0,
  MAP_component = 0L,
  MAP_intensity_rank = 0L,
  MAP_intensity_class = "No exposure",
  MAP_probability = 1,
  classification_uncertainty = 0,
  posterior_low = NA_real_
)]
for (nm in colnames(posterior_ordered)) zero_intensity[, (nm) := NA_real_]
intensity_assignment <- rbindlist(list(zero_intensity, positive_intensity), fill = TRUE)
setorder(intensity_assignment, MAP_intensity_rank, cumulative_EDH30, env_ID)

bic_dt <- as.data.table(as.table(as.matrix(intensity_fit$BIC)))
setnames(bic_dt, c("G", "model", "BIC"))
bic_dt[, G := as.integer(as.character(G))]
bic_dt[, selected := G == intensity_fit$G & model == intensity_fit$modelName]
component_variance <- intensity_fit$parameters$variance$sigmasq
if (length(component_variance) == 1L) component_variance <- rep(component_variance, intensity_fit$G)
component_dt <- data.table(
  component = component_order,
  intensity_rank = seq_along(component_order),
  intensity_class = component_names,
  mixing_probability = intensity_fit$parameters$pro[component_order],
  mean_log1p = intensity_fit$parameters$mean[component_order],
  sd_log1p = sqrt(component_variance[component_order]),
  mean_location_backtransformed = expm1(intensity_fit$parameters$mean[component_order]),
  MAP_n = as.integer(table(factor(map_rank, levels = seq_len(intensity_fit$G))))
)

map_nonlow_envs <- positive_intensity[MAP_intensity_class %in% c("Moderate", "High"), env_ID]
wide <- dcast(bin_dt, env_ID ~ bin, value.var = "value_raw")
bin_cols <- as.character(1:100)
setcolorder(wide, c("env_ID", bin_cols))
x_all <- as.matrix(wide[, ..bin_cols])
rownames(x_all) <- wide$env_ID
x_raw <- x_all[map_nonlow_envs, , drop = FALSE]
row_sd <- apply(x_raw, 1, sd)
if (any(!is.finite(row_sd) | row_sd == 0)) stop("A MAP non-Low profile cannot be shape-standardized.")
x_z <- t(scale(t(x_raw), center = TRUE, scale = TRUE))
n_env <- nrow(x_z)
if (n_env != length(map_nonlow_envs) || any(!is.finite(x_z))) stop("Invalid MAP non-Low shape matrix.")

cor_mat <- cor(t(x_raw), method = "pearson")
if (any(!is.finite(cor_mat))) stop("Correlation matrix contains non-finite values.")
diag(cor_mat) <- 1
dist_mat <- 1 - cor_mat
dist_mat[dist_mat < 0] <- 0
cor_dist <- as.dist(dist_mat)
profile_total <- rowSums(x_raw)
profile_timing <- rowSums(sweep(x_raw, 2, 1:100, "*")) / profile_total

comb2 <- function(v) v * (v - 1) / 2
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  nn <- sum(tab)
  if (nn < 2L) return(NA_real_)
  index <- sum(comb2(tab))
  row_pairs <- sum(comb2(rowSums(tab)))
  col_pairs <- sum(comb2(colSums(tab)))
  expected <- row_pairs * col_pairs / comb2(nn)
  denominator <- 0.5 * (row_pairs + col_pairs) - expected
  if (denominator == 0) return(NA_real_)
  (index - expected) / denominator
}
relabel_by_timing <- function(cl, timing) {
  med <- tapply(timing, cl, median, na.rm = TRUE)
  map <- setNames(seq_along(med), names(sort(med)))
  as.integer(map[as.character(cl)])
}
cluster_centers <- function(x, cl) {
  centers <- rowsum(x, group = cl, reorder = TRUE) / as.numeric(table(cl))
  as.matrix(centers)
}
nearest_center_mse <- function(test_x, centers) {
  d2 <- outer(rowSums(test_x^2), rowSums(centers^2), "+") - 2 * tcrossprod(test_x, centers)
  mean(pmax(0, apply(d2, 1, min))) / ncol(test_x)
}
calc_metrics <- function(cl) {
  k <- uniqueN(cl)
  total_ss <- sum(sweep(x_z, 2, colMeans(x_z), "-")^2)
  within_ss <- sum(vapply(split(seq_len(n_env), cl), function(ii) {
    xx <- x_z[ii, , drop = FALSE]
    sum(sweep(xx, 2, colMeans(xx), "-")^2)
  }, numeric(1)))
  ch <- ((total_ss - within_ss) / (k - 1)) / (within_ss / (n_env - k))
  groups <- split(seq_len(n_env), cl)
  max_intra <- max(vapply(groups, function(ii) {
    if (length(ii) < 2L) return(0)
    max(dist_mat[ii, ii, drop = FALSE])
  }, numeric(1)))
  gids <- seq_along(groups)
  min_inter <- min(unlist(lapply(gids[-length(gids)], function(i) {
    vapply(gids[gids > i], function(j) min(dist_mat[groups[[i]], groups[[j]], drop = FALSE]), numeric(1))
  })))
  sil_width <- cluster::silhouette(cl, cor_dist)[, "sil_width"]
  data.table(
    silhouette = mean(sil_width),
    silhouette_negative_n = sum(sil_width < 0),
    silhouette_negative_pct = 100 * mean(sil_width < 0),
    calinski_harabasz = ch,
    dunn = min_inter / max_intra,
    within_ss = within_ss,
    min_cluster_n = min(table(cl)),
    max_cluster_n = max(table(cl)),
    cluster_sizes = paste(as.integer(table(cl)), collapse = ",")
  )
}

full_assignments <- list()
metric_parts <- list()
for (k in k_values) {
  set.seed(seed + k)
  km <- kmeans(x_z, centers = k, nstart = full_nstart, iter.max = 1000, algorithm = "Hartigan-Wong")
  cl <- relabel_by_timing(km$cluster, profile_timing)
  full_assignments[[as.character(k)]] <- cl
  metric_parts[[length(metric_parts) + 1L]] <- cbind(data.table(k = k), calc_metrics(cl))
}
metrics <- rbindlist(metric_parts)
stability_parts <- vector("list", stability_repeats * length(k_values))
pos <- 0L
stability_n <- floor(stability_fraction * n_env)
for (b in seq_len(stability_repeats)) {
  if (b %% 20L == 0L) message("Stability subsample ", b, "/", stability_repeats)
  set.seed(seed + 1000L + b)
  ii <- sort(sample.int(n_env, stability_n))
  for (k in k_values) {
    set.seed(seed + 100000L + b * 100L + k)
    sub_cl <- kmeans(
      x_z[ii, , drop = FALSE], centers = k, nstart = resample_nstart,
      iter.max = 1000, algorithm = "Hartigan-Wong"
    )$cluster
    pos <- pos + 1L
    stability_parts[[pos]] <- data.table(
      resample_id = b, k = k,
      ari = adjusted_rand(full_assignments[[as.character(k)]][ii], sub_cl)
    )
  }
}
stability_detail <- rbindlist(stability_parts)
stability_summary <- stability_detail[, .(
  stability_mean_ari = mean(ari),
  stability_sd_ari = sd(ari),
  stability_p05_ari = as.numeric(quantile(ari, 0.05)),
  stability_p95_ari = as.numeric(quantile(ari, 0.95))
), by = k]

cv_parts <- vector("list", cv_repeats * cv_folds * length(k_values))
pos <- 0L
for (r in seq_len(cv_repeats)) {
  if (r %% 5L == 0L) message("CV repeat ", r, "/", cv_repeats)
  set.seed(seed + 2000L + r)
  fold_id <- sample(rep(seq_len(cv_folds), length.out = n_env))
  for (fold in seq_len(cv_folds)) {
    test_idx <- which(fold_id == fold)
    train_idx <- which(fold_id != fold)
    train_x <- x_z[train_idx, , drop = FALSE]
    test_x <- x_z[test_idx, , drop = FALSE]
    for (k in k_values) {
      set.seed(seed + 200000L + r * 1000L + fold * 100L + k)
      fit <- kmeans(
        train_x, centers = k, nstart = resample_nstart,
        iter.max = 1000, algorithm = "Hartigan-Wong"
      )
      pos <- pos + 1L
      cv_parts[[pos]] <- data.table(
        cv_repeat = r, fold = fold, k = k,
        n_train = length(train_idx), n_test = length(test_idx),
        cv_mse = nearest_center_mse(test_x, fit$centers)
      )
    }
  }
}
cv_detail <- rbindlist(cv_parts)
cv_repeat_summary <- cv_detail[, .(cv_mse = weighted.mean(cv_mse, n_test)), by = .(cv_repeat, k)]
cv_summary <- cv_repeat_summary[, .(
  cv_mse_mean = mean(cv_mse),
  cv_mse_sd = sd(cv_mse),
  cv_mse_se = sd(cv_mse) / sqrt(.N)
), by = k]
metrics <- merge(metrics, stability_summary, by = "k")
metrics <- merge(metrics, cv_summary, by = "k")
setorder(metrics, k)

k1_wss <- sum(sweep(x_z, 2, colMeans(x_z), "-")^2)
wss_dt <- rbind(data.table(k = 1L, within_ss = k1_wss), metrics[, .(k, within_ss)])
wss_dt[, `:=`(
  wss_drop_from_previous = shift(within_ss) - within_ss,
  wss_drop_from_previous_pct = 100 * (shift(within_ss) - within_ss) / shift(within_ss)
)]
metrics <- merge(
  metrics,
  wss_dt[k >= 2L, .(k, wss_drop_from_previous, wss_drop_from_previous_pct)],
  by = "k"
)
setorder(metrics, k)
shape_assignment <- data.table(
  env_ID = rownames(x_raw),
  cumulative_EDH30 = positive_intensity[match(rownames(x_raw), env_ID), cumulative_EDH30],
  MAP_intensity_class = positive_intensity[match(rownames(x_raw), env_ID), MAP_intensity_class],
  MAP_probability = positive_intensity[match(rownames(x_raw), env_ID), MAP_probability],
  posterior_low = positive_intensity[match(rownames(x_raw), env_ID), posterior_low],
  exposure_timing_bin = profile_timing
)
for (k in k_values) shape_assignment[, (sprintf("kmeans_k%d", k)) := full_assignments[[as.character(k)]]]

final_assignment <- merge(
  intensity_assignment,
  shape_assignment[, .(env_ID, exposure_timing_bin, kmeans_k3)],
  by = "env_ID", all.x = TRUE, sort = FALSE
)
final_assignment[, analysis_group := fcase(
  MAP_intensity_class == "No exposure", "No exposure",
  MAP_intensity_class == "Low", "MAP Low (excluded from shape clustering)",
  kmeans_k3 == 1L, "Early heat pattern",
  kmeans_k3 == 2L, "Middle heat pattern",
  kmeans_k3 == 3L, "Late heat pattern",
  default = NA_character_
)]
if (anyNA(final_assignment$analysis_group)) stop("Incomplete all-environment MAP classification.")
if (any(final_assignment[MAP_intensity_class %in% c("No exposure", "Low"), !is.na(kmeans_k3)])) {
  stop("An excluded environment received a K-means class.")
}
# 3. 仅输出 K=3 环境分类；零暴露与 MAP Low 不进入时序聚类。
classification <- final_assignment[, .(
  env_ID, cumulative_EDH30, log1p_cumulative_EDH30,
  MAP_component, MAP_intensity_rank, MAP_intensity_class,
  MAP_probability, classification_uncertainty,
  posterior_low, posterior_moderate, posterior_high,
  exposure_timing_bin, kmeans_k3
)]
classification[, Env_Type := fcase(
  MAP_intensity_class %in% c("No exposure", "Low"), "Low heat exposure",
  kmeans_k3 == 1L, "Early heat exposure",
  kmeans_k3 == 2L, "Middle heat exposure",
  kmeans_k3 == 3L, "Late heat exposure"
)]
stopifnot(
  nrow(classification) == 856L,
  n_env == 610L,
  !anyNA(classification$Env_Type),
  identical(as.integer(table(factor(classification$Env_Type, levels = c(
    "Low heat exposure", "Early heat exposure", "Middle heat exposure", "Late heat exposure"
  )))), c(246L, 270L, 201L, 139L))
)
fwrite(classification, file.path(out, "environment_classification.csv"))
# 4. 统一紧凑版式；仅保留坐标、图例和面板编号，不添加标题或页脚。
theme_pub <- function() {
  theme_classic(base_size = 8, base_family = "Arial") + theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1.2, "mm"),
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(size = 8),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    plot.tag = element_text(face = "bold", size = 10),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    plot.margin = margin(4, 5, 3, 4),
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width = unit(5, "mm"),
    legend.text = element_text(size = 7),
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0)
  )
}

density_x <- seq(min(intensity_x), max(intensity_x), length.out = 1000)
density_parts <- rbindlist(lapply(seq_len(intensity_fit$G), function(j) {
  data.table(
    x = density_x,
    density = intensity_fit$parameters$pro[j] * dnorm(
      density_x,
      mean = intensity_fit$parameters$mean[j],
      sd = sqrt(component_variance[j])
    ),
    component = component_names[component_rank[as.character(j)]]
  )
}))
density_total <- density_parts[, .(density = sum(density)), by = x]
intensity_colors <- setNames(
  c("#0072B2", "#E69F00", "#D55E00", "#009E73")[seq_len(intensity_fit$G)],
  component_names
)
p_density <- ggplot(positive_intensity, aes(x = log1p_cumulative_EDH30)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = "grey88", color = "white", linewidth = 0.2) +
  geom_line(data = density_parts, aes(x = x, y = density, color = component), linewidth = 0.65) +
  geom_line(data = density_total, aes(x = x, y = density), linewidth = 0.8, color = "black") +
  scale_color_manual(values = intensity_colors, breaks = component_names) +
  labs(tag = "a", x = "ln(1+EDH30)", y = "Density", color = NULL) +
  theme_pub()
p_bic <- ggplot(bic_dt, aes(x = G, y = BIC, color = model)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.6) +
  geom_point(data = bic_dt[selected == TRUE], shape = 21, size = 3, stroke = 0.8, fill = "white") +
  scale_x_continuous(breaks = 1:5) +
  scale_color_manual(values = c(E = "#3366CC", V = "#D55E00"),
                     labels = c(E = "Equal variance", V = "Unequal variance")) +
  labs(tag = "b", x = "Number of components", y = "BIC", color = NULL) +
  theme_pub() +
  theme(legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0))
p_uncertainty <- ggplot(
  positive_intensity,
  aes(x = cumulative_EDH30, y = classification_uncertainty, color = MAP_intensity_class)
) +
  geom_point(size = 1, alpha = 0.55) +
  scale_x_log10(labels = label_number()) +
  scale_color_manual(values = intensity_colors, breaks = component_names) +
  labs(
    tag = "c",
    x = "EDH30 (\u00b0C\u00b7h)",
    y = "Classification uncertainty", color = NULL
  ) +
  theme_pub()
p_status <- ggplot(
  positive_intensity,
  aes(
    x = factor(MAP_intensity_class, levels = component_names),
    y = cumulative_EDH30,
    fill = MAP_intensity_class
  )
) +
  geom_boxplot(width = 0.65, linewidth = 0.35, outlier.size = 0.5) +
  scale_y_log10(labels = label_number()) +
  scale_fill_manual(values = intensity_colors, guide = "none") +
  scale_x_discrete(labels = c(
    "Low" = sprintf("Low\n(n=%d)", sum(map_label == "Low")),
    "Moderate" = sprintf("Moderate\n(n=%d)", sum(map_label == "Moderate")),
    "High" = sprintf("High\n(n=%d)", sum(map_label == "High"))
  )) +
  labs(tag = "d", x = NULL, y = "EDH30 (\u00b0C\u00b7h)") +
  theme_pub() + theme(axis.text.x = element_text(size = 7))
p_intensity <- wrap_plots(p_density, p_bic, p_uncertainty, p_status, ncol = 2)

# 仅展示分离程度与重抽样稳定性；其余 K 值诊断保留于 Source Data。
p_mean_silhouette <- ggplot(metrics, aes(x = k, y = silhouette)) +
  geom_vline(xintercept = 3, linetype = "dashed", linewidth = 0.3, color = "#D55E00") +
  geom_line(linewidth = 0.5, color = "#0072B2") +
  geom_point(size = 1.4, color = "#0072B2") +
  scale_x_continuous(breaks = 2:8) +
  labs(tag = "b", x = "K", y = "Mean silhouette") +
  theme_pub()
p_stability <- ggplot(metrics, aes(x = k, y = stability_mean_ari)) +
  geom_vline(xintercept = 3, linetype = "dashed", linewidth = 0.3, color = "#D55E00") +
  geom_line(linewidth = 0.5, color = "#0072B2") +
  geom_point(size = 1.4, color = "#0072B2") +
  scale_x_continuous(breaks = 2:8) +
  labs(tag = "c", x = "K", y = "Subsample stability (ARI)") +
  theme_pub()

# 5. 从 K=3 的 100-bin 聚类输入计算轮廓系数及 PCA 投影。
assignment <- shape_assignment[, .(env_ID, cluster = kmeans_k3)]
cluster_labels <- c("1" = "Early heat pattern", "2" = "Middle heat pattern", "3" = "Late heat pattern")
cluster_colors <- c("1" = "#0072B2", "2" = "#E69F00", "3" = "#009E73")
attr(x_z, "scaled:center") <- NULL
attr(x_z, "scaled:scale") <- NULL
pca_fit <- prcomp(x_z, center = TRUE, scale. = FALSE)
variance_pct <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
pca_scores <- as.data.table(pca_fit$x, keep.rownames = "env_ID")
pca_scores <- merge(pca_scores, assignment, by = "env_ID", sort = FALSE)
pca_scores[, cluster_factor := factor(cluster, levels = 1:3, labels = unname(cluster_labels))]
pca_loadings <- as.data.table(copy(pca_fit$rotation), keep.rownames = "bin")
pca_loadings[, bin := as.integer(bin)]
pca_variance <- data.table(
  component = paste0("PC", seq_along(variance_pct)),
  component_number = seq_along(variance_pct),
  variance_pct = variance_pct,
  cumulative_variance_pct = cumsum(variance_pct)
)

pca_centroids <- pca_scores[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = .(cluster, cluster_factor)]
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = factor(cluster))) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.45, show.legend = FALSE) +
  geom_point(size = 1.15, alpha = 0.55) +
  geom_point(
    data = pca_centroids, aes(x = PC1, y = PC2, color = factor(cluster)),
    shape = 4, stroke = 1.1, size = 4, inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_colors, labels = c("Early", "Middle", "Late"), name = NULL) +
  labs(
    tag = "e",
    x = sprintf("PC1 (%.1f%%)", variance_pct[1]),
    y = sprintf("PC2 (%.1f%%)", variance_pct[2])
  ) +
  coord_equal(
    xlim = range(c(pca_scores$PC1, pca_scores$PC2)) + c(-1, 1),
    ylim = range(c(pca_scores$PC1, pca_scores$PC2)) + c(-1, 1)
  ) +
  theme_pub() +
  theme(legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0))

centers_z <- do.call(rbind, lapply(1:3, function(k) {
  colMeans(x_z[assignment$cluster == k, , drop = FALSE])
}))
projected_centers <- predict(pca_fit, newdata = centers_z)
stopifnot(max(abs(projected_centers[, 1:2] - as.matrix(pca_centroids[order(cluster), .(PC1, PC2)]))) < 1e-8)
observed_wss <- sum((x_z - centers_z[assignment$cluster, ])^2)
stopifnot(abs(observed_wss - metrics[k == 3, within_ss]) < 1e-6)
elbow_dt <- rbind(
  data.table(k = 1L, within_ss = sum(scale(x_z, center = TRUE, scale = FALSE)^2)),
  metrics[, .(k, within_ss)]
)
distance_matrix <- 1 - cor(t(x_raw))
distance_matrix[distance_matrix < 0] <- 0
diag(distance_matrix) <- 0
sil <- cluster::silhouette(assignment$cluster, as.dist(distance_matrix))
sil_dt <- data.table(env_ID = assignment$env_ID, cluster = assignment$cluster,
                     silhouette = sil[, "sil_width"])
stopifnot(abs(mean(sil_dt$silhouette) - metrics[k == 3, silhouette]) < 1e-8)
setorder(sil_dt, cluster, silhouette)
sil_dt[, display_y := .I + 12L * (cluster - 1L)]
sil_labels <- sil_dt[, .(display_y = mean(display_y), n = .N), by = cluster]
diagnostic_theme <- theme_pub() + theme(aspect.ratio = 1)
p_elbow <- ggplot(elbow_dt, aes(k, within_ss)) +
  geom_line(color = "#326C82", linewidth = 0.5) +
  geom_point(color = "#326C82", size = 1.8) +
  geom_point(data = elbow_dt[k == 3], shape = 1, size = 3.5, color = "#B33B3B", stroke = 0.8) +
  scale_x_continuous(breaks = 1:8) +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(tag = "a", x = "K", y = "Within-cluster\nsum of squares") +
  theme_pub()
p_silhouette <- ggplot(sil_dt) +
  geom_segment(aes(x = 0, xend = silhouette, y = display_y, yend = display_y,
                   color = factor(cluster)), linewidth = 0.25, show.legend = FALSE) +
  geom_vline(xintercept = mean(sil_dt$silhouette), linetype = "dashed", color = "#B33B3B", linewidth = 0.4) +
  scale_color_manual(values = cluster_colors) +
  scale_y_reverse(breaks = sil_labels$display_y,
                  labels = sprintf("%s\n(n=%d)", c("Early", "Middle", "Late"), sil_labels$n)) +
  scale_x_continuous(
    limits = c(min(-0.1, min(sil_dt$silhouette) - 0.02), max(0.5, max(sil_dt$silhouette) + 0.02)),
    breaks = c(-0.1, 0, 0.2, 0.4)
  ) +
  labs(tag = "d", x = "Silhouette coefficient", y = NULL) + diagnostic_theme
p_k3_diagnostics <- wrap_plots(
  p_elbow,
  p_mean_silhouette,
  p_stability,
  p_silhouette,
  p_pca,
  design = "AABBCC\nDDDEEE",
  heights = c(0.65, 1)
)


ggsave(file.path(out, "Supplementary_intensity_classification.pdf"),
       p_intensity, width = 180, height = 118, units = "mm", device = cairo_pdf)
ggsave(file.path(out, "Supplementary_intensity_classification.tif"),
       p_intensity, width = 180, height = 118, units = "mm",
       dpi = 300, device = "tiff", type = "cairo", compression = "lzw")
ggsave(file.path(out, "Supplementary_K3_diagnostics.pdf"),
       p_k3_diagnostics, width = 180, height = 138, units = "mm", device = cairo_pdf)

# 6. Source Data 只保留实际输入、作图数值与拟合参数，全部使用英文。
mixture_curves <- rbind(
  density_parts[, .(x, density, component)],
  density_total[, .(x, density, component = "Mixture")]
)
histogram <- as.data.table(ggplot_build(p_density)$data[[1]])[, .(xmin, xmax, count, density)]
boxplot_stats <- as.data.table(ggplot_build(p_status)$data[[1]])[
  , .(intensity_class = component_names[group], ymin, lower, middle, upper, ymax)
]
shape_input <- as.data.table(x_z, keep.rownames = "env_ID")
setnames(shape_input, as.character(1:100), sprintf("bin_%03d", 1:100))
centroids <- as.data.table(cluster_centers(x_z, full_assignments[["3"]]), keep.rownames = "cluster")
setnames(centroids, as.character(1:100), sprintf("bin_%03d", 1:100))
parameters <- data.table(
  selected_shape_K = 3L,
  intensity_model = intensity_fit$modelName,
  intensity_G = intensity_fit$G,
  intensity_loglik = intensity_fit$loglik,
  intensity_BIC = intensity_fit$bic,
  intensity_input = "log1p(cumulative_EDH30), positive environments only",
  eligibility = "MAP Moderate or High; no posterior threshold",
  intensity_environment_n = nrow(positive),
  zero_environment_n = length(zero_envs),
  shape_environment_n = n_env,
  shape_standardization = "Within-environment sample mean and sample SD across 100 bins",
  kmeans_algorithm = "Hartigan-Wong",
  seed_base = seed,
  full_seed = "seed_base + K",
  full_nstart = full_nstart,
  resample_nstart = resample_nstart,
  iter_max = 1000L,
  silhouette_distance = "1 - Pearson correlation",
  stability_fraction = stability_fraction,
  stability_repeats = stability_repeats,
  cv_repeats = cv_repeats,
  cv_folds = cv_folds,
  PCA_center = TRUE,
  PCA_scale = FALSE,
  ellipse = "Descriptive 95% normal contour; not a cluster-validity confidence region",
  R_version = as.character(getRversion()),
  mclust_version = as.character(packageVersion("mclust")),
  cluster_version = as.character(packageVersion("cluster"))
)
sheets <- list(
  Environment_data = classification,
  Shape_input = shape_input,
  Mixture_parameters = component_dt,
  Mixture_curves = mixture_curves,
  Histogram_bins = histogram,
  Boxplot_log10_stats = boxplot_stats,
  Intensity_BIC = bic_dt,
  K_selection = metrics,
  Elbow_curve = elbow_dt,
  K3_silhouette = sil_dt,
  K3_centroids = centroids,
  Stability_data = stability_detail,
  Cross_validation = cv_detail,
  PCA_scores = pca_scores[, .(env_ID, cluster, PC1, PC2)],
  PCA_loadings = pca_loadings[, .(bin, PC1, PC2)],
  PCA_variance = pca_variance[component_number <= 2L],
  PCA_centroids = pca_centroids[, .(cluster, PC1, PC2)],
  Model_parameters = parameters
)
wb <- createWorkbook()
header_style <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold", fgFill = "#E8EDF2")
modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
for (name in names(sheets)) {
  data <- sheets[[name]]
  text_values <- unlist(lapply(data, function(x) if (is.character(x) || is.factor(x)) as.character(x) else character()))
  stopifnot(!any(grepl("[^\\x01-\\x7F]", text_values, perl = TRUE), na.rm = TRUE))
  addWorksheet(wb, name)
  writeData(wb, name, data, headerStyle = header_style, withFilter = TRUE, keepNA = FALSE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_along(data), widths = 18)
}
saveWorkbook(wb, file.path(out, "Supplementary_Classification_Source_Data.xlsx"), overwrite = TRUE)
message(format(Sys.time(), "%F %T"), " | Classification, diagnostics and Source Data completed")
print(classification[, .N, by = Env_Type])
print(metrics[k == 3, .(k, silhouette, stability_mean_ari, within_ss)])

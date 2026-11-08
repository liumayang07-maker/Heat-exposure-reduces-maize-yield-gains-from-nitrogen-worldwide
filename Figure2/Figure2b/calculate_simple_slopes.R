#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages(library(metafor))

# 1. 读取正式模型，不重新拟合；斜率单位为 lnRR / ln(1 + EDH30)。
entry <- grep("^--file=", commandArgs(), value = TRUE)
output_dir <- dirname(normalizePath(sub("^--file=", "", entry)))
fit <- readRDS(file.path(output_dir, "model_REML.rds"))$fit
transformations <- read.csv(file.path(output_dir, "transformation_parameters.csv"))
groups <- c("0_100", "100_200", "200_300", "300_400", "400_600")
confidence_level <- 0.95
stopifnot(fit$test == "z", fit$method == "REML")
stopifnot(transformations$divisor[transformations$target == "heat_c"] == 1)

# 2. 构造各组线性对比，保留主效应与交互项之间的协方差。
beta <- coef(fit)
covariance <- vcov(fit)
contrast <- matrix(0, nrow = length(groups), ncol = length(beta),
                   dimnames = list(groups, names(beta)))
contrast[, "heat_c"] <- 1
for (index in 2:length(groups)) {
  contrast[index, paste0("heat_c:N_group", groups[index])] <- 1
}
slope <- as.numeric(contrast %*% beta)
standard_error <- sqrt(diag(contrast %*% covariance %*% t(contrast)))
critical_value <- qnorm((1 + confidence_level) / 2)
z_value <- slope / standard_error
results <- data.frame(
  N_group = groups,
  simple_slope = slope,
  SE = standard_error,
  CI_lower = slope - critical_value * standard_error,
  CI_upper = slope + critical_value * standard_error,
  z = z_value,
  P_two_sided = 2 * pnorm(-abs(z_value)),
  row.names = NULL
)

# 3. 与 metafor 的单组 Wald 对比核验；P 值不作多重检验校正。
for (index in seq_along(groups)) {
  check <- anova(fit, X = contrast[index, , drop = FALSE])
  stopifnot(isTRUE(all.equal(results$P_two_sided[index], as.numeric(check$QMp),
                            tolerance = 1e-10)))
}
stopifnot(nrow(results) == 5L, all(is.finite(as.matrix(results[, -1]))))
write.csv(results, file.path(output_dir, "simple_slopes_REML.csv"), row.names = FALSE)
print(results, digits = 8, row.names = FALSE)

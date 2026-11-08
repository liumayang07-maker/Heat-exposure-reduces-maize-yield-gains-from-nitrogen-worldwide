#!/usr/bin/env Rscript
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
# 复用临时拟合、完整验证及核心 CSV 汇总流程；成功后只保留核心 CSV。
entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- normalizePath(file.path(dirname(entry), "../.."), winslash = "/")
source(file.path(root, "Figure2/SD_sensitivity/run_sensitivity.R"), encoding = "UTF-8")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(all(args %in% "--resume"))
run_sensitivity("Figure3", root, resume = "--resume" %in% args)

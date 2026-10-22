# Compare the benchmark predictions made at different nsim (07_predict_1km_tiled.R with W2W_SAMPLE).
#   Rscript R/bench_nsim_compare.R model2_baseline_2010_1km      (stem = prediction name)
if (!exists("DIR_DATA")) source("R/00_config.R")
stem <- commandArgs(trailingOnly = TRUE)[1]
fs <- list.files(file.path(DIR_OUT, "bench_1km"), pattern = paste0("^", stem, "_nsim[0-9]+\\.rds$"), full.names = TRUE)
if (length(fs) < 2) stop("need benchmark files for at least two nsim values")
r <- lapply(fs, readRDS); ns <- vapply(r, function(x) x$nsim, 1L); r <- r[order(ns)]; ns <- sort(ns)
ref <- r[[length(r)]]; cat("reference = nsim", max(ns), "\n")
for (i in seq_along(r)) {
  d <- r[[i]]$mean - ref$mean; tot_d <- rowSums(r[[i]]$mean) - rowSums(ref$mean)
  cat(sprintf("nsim %3d | %5.1f ms/cell | total BA: mean |diff| %.3f (ref mean %.2f) | per-species mean |diff| %.4f = %.1f%% of mean se\n",
              ns[i], 1000 * r[[i]]$secs_per_cell, mean(abs(tot_d)), mean(rowSums(ref$mean)),
              mean(abs(d)), 100 * mean(abs(d)) / mean(ref$se)))
}

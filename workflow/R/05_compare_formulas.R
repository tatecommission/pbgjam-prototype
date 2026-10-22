# =============================================================================
# 05_compare_formulas.R -- fit all 8 Model 1 formulas (vegDRI removed, see
# 00_config.R header) and all 16 Model 2 formulas at an exploratory MCMC
# budget (ng = NG_EXPLORE, burnin = BURNIN_EXPLORE), score each with DIC
# and WAIC ONLY.
#
# PARALLEL + RESUMABLE (this revision):
#   - Each formula fits independently, so all 24 are dispatched across
#     cores with parallel::mclapply (fork-based -- xdata_c/ymat are loaded
#     ONCE in this process and shared copy-on-write across forks, so there
#     is no per-core data-reload cost). mclapply is macOS/Linux only; it
#     will not parallelize on Windows (falls back to serial there).
#   - Every finished fit -- success OR failure -- is checkpointed to
#     out/formula_comparison/fits/<name>.rds THE MOMENT it completes, via
#     a write-to-tempfile-then-rename so a checkpoint file is never left
#     half-written if interrupted mid-save.
#   - On every run, already-checkpointed formulas are skipped entirely --
#     you can Ctrl+C, close the laptop, or just stop after N formulas and
#     re-run this script later; it only fits whatever's missing.
#   - Final ranking/report is rebuilt from ALL checkpoint files on disk
#     each time, so partial runs still produce a partial report, and the
#     final run (once everything's cached) produces the complete one.
#
# NO holdout R2, NO climate-transfer RMSE, NO monotonicity fraction -- all
# explicitly dropped per instruction; DIC/WAIC only.
#
# DIC/WAIC extraction is DEFENSIVE (gjam's exact output field was not
# verified against your installed package version) -- see get_dic()/
# get_waic() below; if both come back NA for everything, check
# names(fit)/names(fit$fit) on any one saved fit object and fix the
# extractor.
# =============================================================================

# NUMBER OF CORES
N_CORES <- 4

if (!exists("DIR_DATA")) source("R/00_config.R")
library(gjam)
library(parallel)

xdata_c <- readRDS(file.path(DIR_DATA, "xdata_centered.rds"))
ymat    <- as.matrix(readRDS(file.path(DIR_DATA, "ydata_aligned.rds")))
FORMULAS_M1 <- readRDS(file.path(DIR_DATA, "candidate_formulas_m1.rds"))
FORMULAS_M2 <- readRDS(file.path(DIR_DATA, "candidate_formulas_m2.rds"))

S <- ncol(ymat)

## ---- checkpoint directory ---------------------------------------------------
FITS_DIR <- file.path(DIR_OUT, "formula_comparison", "fits")
dir.create(FITS_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(DIR_OUT, "formula_comparison", "run_log.txt")

log_line <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)),
      file = LOG_FILE, append = TRUE)
}

## ---- how many cores -----------------------------------------------------
## Override by setting PARALLEL_CORES in R/00_config.R if you want a
## specific number (e.g. to leave more headroom for other work). Default:
## all cores minus 1, capped at 6 so a big-core machine doesn't try to run
## 20 gjam fits in memory simultaneously.

is_forking_ok <- .Platform$OS.type == "unix"
if (!is_forking_ok)
  message("NOTE: mclapply forking isn't available on Windows -- this will run serially.")
message(sprintf("using %d core(s) for parallel fitting", if (is_forking_ok) N_CORES else 1))

## ---- DIC / WAIC extraction, both defensive ---------------------------------
get_dic <- function(fit) {
  candidates <- list(
    tryCatch(fit$fit$DIC, error = function(e) NULL),
    tryCatch(fit$DIC, error = function(e) NULL),
    tryCatch(fit$fit$dic, error = function(e) NULL)
  )
  val <- Find(Negate(is.null), candidates)
  if (is.null(val)) NA_real_ else as.numeric(val)
}

get_waic <- function(fit) {
  candidates <- list(
    tryCatch(fit$fit$WAIC, error = function(e) NULL),
    tryCatch(fit$WAIC, error = function(e) NULL),
    tryCatch(fit$fit$waic, error = function(e) NULL)
  )
  val <- Find(Negate(is.null), candidates)
  if (is.null(val)) NA_real_ else as.numeric(val)
}

## ---- atomic checkpoint write ------------------------------------------------
save_checkpoint <- function(path, obj) {
  tmp <- paste0(path, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, path)   # rename is atomic on the same filesystem
}

## ---- fit + score + checkpoint one formula -----------------------------------
fit_one <- function(name, form, model_label, save_full_fit = TRUE) {
  ckpt <- file.path(FITS_DIR, paste0(name, ".rds"))
  if (file.exists(ckpt)) return(invisible(NULL))   # already done, skip
  
  log_line("[%s] starting (ng=%d, burnin=%d)", name, NG_EXPLORE, BURNIN_EXPLORE)
  t0 <- Sys.time()
  
  ml <- list(ng = NG_EXPLORE, burnin = BURNIN_EXPLORE,
             typeNames = rep("CA", S),
             reductList = list(r = REDUCT_R, N = REDUCT_N),
             PREDICTX = FALSE)
  
  fit <- tryCatch(
    gjam(form, xdata = xdata_c, ydata = ymat, modelList = ml),
    error = function(e) e)
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  
  if (inherits(fit, "error")) {
    log_line("[%s] FAILED after %.1f min: %s", name, elapsed, conditionMessage(fit))
    result <- list(name = name, model = model_label, success = FALSE,
                   dic = NA_real_, waic = NA_real_,
                   error = conditionMessage(fit), elapsed_min = elapsed)
  } else {
    dic  <- get_dic(fit)
    waic <- get_waic(fit)
    log_line("[%s] done in %.1f min (DIC=%.2f, WAIC=%.2f)", name, elapsed, dic, waic)
    result <- list(name = name, model = model_label, success = TRUE,
                   dic = dic, waic = waic, elapsed_min = elapsed,
                   fit = if (save_full_fit) fit else NULL)
  }
  
  save_checkpoint(ckpt, result)
  invisible(NULL)
}

## ---- build the combined task queue (both models, so cores are shared -----
## across the full 24 rather than doing 8 then 16 sequentially) ---------------
tasks <- c(
  Map(function(n, f) list(name = n, form = f, model = "M1"),
      names(FORMULAS_M1), FORMULAS_M1),
  Map(function(n, f) list(name = n, form = f, model = "M2"),
      names(FORMULAS_M2), FORMULAS_M2)
)

pending <- Filter(function(tk) !file.exists(file.path(FITS_DIR, paste0(tk$name, ".rds"))), tasks)
message(sprintf("%d of %d formulas already checkpointed -- %d remaining to fit",
                length(tasks) - length(pending), length(tasks), length(pending)))

if (length(pending) > 0) {
  if (is_forking_ok) {
    invisible(mclapply(pending, function(tk) {
      fit_one(tk$name, tk$form, tk$model)
    }, mc.cores = N_CORES, mc.preschedule = FALSE))
  } else {
    for (tk in pending) fit_one(tk$name, tk$form, tk$model)
  }
} else {
  message("nothing left to fit -- all checkpoints already present")
}

## =============================================================================
## Aggregate report from EVERY checkpoint currently on disk (not just this
## run's) -- so a partial run still gives you a partial table, and the
## final run gives you the complete one.
## =============================================================================
ckpt_files <- list.files(FITS_DIR, pattern = "\\.rds$", full.names = TRUE)
all_results <- lapply(ckpt_files, readRDS)

tab <- do.call(rbind, lapply(all_results, function(r) data.frame(
  model   = r$model,
  formula = r$name,
  success = r$success,
  dic     = r$dic,
  waic    = r$waic,
  elapsed_min = r$elapsed_min,
  stringsAsFactors = FALSE
)))
tab <- tab[order(tab$model, tab$formula), ]

n_total <- length(tasks)
n_done  <- nrow(tab)
if (n_done < n_total)
  message(sprintf("\nNOTE: %d of %d formulas fitted so far -- re-run this script later ",
                  n_done, n_total),
          "to pick up the rest (already-done ones will be skipped).")

failed <- tab[!tab$success, ]
if (nrow(failed) > 0) {
  message("\nThe following formulas FAILED and were NOT included in ranking ",
          "(see run_log.txt / the checkpoint's $error field for details):")
  print(failed[, c("model", "formula")])
}

tab_ok <- tab[tab$success, ]

if (nrow(tab_ok) > 0 && all(is.na(tab_ok$dic)))
  warning("DIC is NA for every successful fit -- get_dic() likely doesn't match ",
          "your installed gjam version's output structure. Check names(fit)/",
          "names(fit$fit) on any saved fits/<name>.rds$fit and fix get_dic() above.")
if (nrow(tab_ok) > 0 && all(is.na(tab_ok$waic)))
  warning("WAIC is NA for every successful fit -- either get_waic() doesn't match ",
          "your installed gjam version, or this version doesn't expose WAIC at ",
          "all. DIC alone is still usable.")

## ---- rank WITHIN each model (comparing M1 to M2 by DIC/WAIC isn't ---------
## meaningful -- they're fit on different predictor sets for different
## purposes, not competing to explain the same data the same way)
tab_m1_sorted <- tab_ok[tab_ok$model == "M1", ]
tab_m1_sorted <- tab_m1_sorted[order(tab_m1_sorted$dic), ]
tab_m2_sorted <- tab_ok[tab_ok$model == "M2", ]
tab_m2_sorted <- tab_m2_sorted[order(tab_m2_sorted$dic), ]

cat("\n--- Model 1, ranked by DIC (lower = better) ---\n")
print(tab_m1_sorted[, c("formula", "dic", "waic", "elapsed_min")])
cat("\n--- Model 2, ranked by DIC (lower = better) ---\n")
print(tab_m2_sorted[, c("formula", "dic", "waic", "elapsed_min")])

if (n_done == n_total) {
  cat("\nAll formulas fitted. Check whether DIC and WAIC agree on the same winner\n",
      "within each model -- if they disagree, that's worth reading closely rather\n",
      "than picking whichever you prefer. At 8 (Model 1) and 16 (Model 2)\n",
      "formulas, read the full table, not just the top row.\n", sep = "")
}

write.csv(tab, file.path(DIR_OUT, "formula_comparison", "scores.csv"), row.names = FALSE)
message("\nwrote out/formula_comparison/scores.csv (", n_done, "/", n_total, " formulas)")
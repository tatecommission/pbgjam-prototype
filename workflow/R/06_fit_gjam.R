# =============================================================================
# 06_fit_gjam.R -- full-budget production fit of BOTH models' winning
# formulas. Runs TWICE (once per model) -- this is a real, deliberate
# increase in compute cost over the pre-rewrite single-model version,
# inherent to the two-model design, not something to shortcut.
#
# Set FORMULA_NAME_M1 / FORMULA_NAME_M2 after inspecting
# out/formula_comparison/scores.csv (05_compare_formulas.R). Override only
# with a stated reason if DIC and WAIC disagree on the winner.
# =============================================================================

source("R/00_config.R")
Rcpp::sourceCpp(file.path(GJAM_SRC, "src", "cppFns.cpp"))
for (f in list.files(file.path(GJAM_SRC, "R"), full.names = TRUE, pattern = "\\.[Rr]$"))
  source(f)

FORMULA_NAME_M1 <- "M1_01_sand30_precipGPM_LST"               # placeholder -- replace with 05's actual M1 winner
FORMULA_NAME_M2 <- "M2_01_sand30_meanTemp_annualPrec"        # placeholder -- replace with 05's actual M2 winner

FORMULAS_M1 <- readRDS(file.path(DIR_DATA, "candidate_formulas_m1.rds"))
FORMULAS_M2 <- readRDS(file.path(DIR_DATA, "candidate_formulas_m2.rds"))
CENTERS     <- readRDS(file.path(DIR_DATA, "predictor_centers.rds"))
xdata_c     <- readRDS(file.path(DIR_DATA, "xdata_centered.rds"))
ymat        <- as.matrix(readRDS(file.path(DIR_DATA, "ydata_aligned.rds")))

form_m1 <- FORMULAS_M1[[FORMULA_NAME_M1]]
form_m2 <- FORMULAS_M2[[FORMULA_NAME_M2]]
stopifnot("FORMULA_NAME_M1 not found -- check out/formula_comparison/scores.csv" = !is.null(form_m1))
stopifnot("FORMULA_NAME_M2 not found -- check out/formula_comparison/scores.csv" = !is.null(form_m2))

if (!DEM_COLS_READY) stop("DEM_COLS_READY is FALSE -- both formulas reference north_component/east_component.")
if (!MODEL1_COLS_READY) stop("MODEL1_COLS_READY is FALSE -- Model 1's formula references columns that are still NA.")
if (!MODEL2_COLS_READY) stop("MODEL2_COLS_READY is FALSE -- Model 2's formula references columns that are still NA.")

S <- ncol(ymat)

fit_one <- function(form, label) {
  message(sprintf("\n=== fitting %s (%s): %d plots x %d species, ng=%d ===",
                  label, deparse(form), nrow(ymat), S, NG))
  ml <- list(ng = NG, burnin = BURNIN, typeNames = rep("CA", S),
            reductList = list(r = REDUCT_R, N = REDUCT_N), PREDICTX = FALSE)
  t0 <- Sys.time()
  fit <- gjam(formula = form, xdata = xdata_c, ydata = ymat, modelList = ml)
  message(sprintf("elapsed: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  fit
}

fit_m1 <- fit_one(form_m1, "Model 1 (remote sensing)")
fit_m2 <- fit_one(form_m2, "Model 2 (climate-only)")

## ---- store each fit with its centering constants + formula name ----------
# CENTERS is the SAME object for both models (built once in 04, covers all
# quadratic-eligible variables across both formula sets), so it's fine to
# bundle the whole thing with each fit rather than trying to subset it.
saveRDS(list(fit = fit_m1, formula_name = FORMULA_NAME_M1, formula = form_m1, centers = CENTERS),
       file.path(DIR_OUT, "fits", "model1_gjam_out.rds"))
saveRDS(list(fit = fit_m2, formula_name = FORMULA_NAME_M2, formula = form_m2, centers = CENTERS),
       file.path(DIR_OUT, "fits", "model2_gjam_out.rds"))

message("\nwrote out/fits/model1_gjam_out.rds and out/fits/model2_gjam_out.rds")

## ---- quick diagnostics for both --------------------------------------------
cat("\n--- Model 1 summary ---\n"); summary(fit_m1)
cat("\n--- Model 2 summary ---\n"); summary(fit_m2)

gjamPlot(fit_m1, plotPars = list(PLOTALLY = TRUE, SAVEPLOTS = TRUE,
                                 outFolder = file.path(DIR_OUT, "fits", "gjamPlots_model1")))
gjamPlot(fit_m2, plotPars = list(PLOTALLY = TRUE, SAVEPLOTS = TRUE,
                                 outFolder = file.path(DIR_OUT, "fits", "gjamPlots_model2")))
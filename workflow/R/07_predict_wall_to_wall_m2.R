# =============================================================================
# 07_predict_wall_to_wall_m2.R -- Model 2's CMIP6 wall-to-wall prediction, both
# stages: historical baseline (scenario "historical", years <= 2014) and future
# scenarios (ssp245 / ssp585, 2015-2100). SAME formula and coefficients as the
# fit -- no reduced model.
#
# CHANGES in this version (raster pipeline) -- same as 07_predict_wall_to_wall_m1.R:
#   GeoTIFF inputs joined by cell number, whole-US land cells (see LAND_FRAC_MIN /
#   FOREST_FRAC_MIN), chunked prediction
#   with progress/ETA, finished (scenario, year) outputs are skipped on re-run,
#   library(gjam), gjam-safe names, "Ctr" suffix, fitted centres.
#
# Prediction years are DECADES (PRED_YEARS_5YR = seq(1950, 2100, by = 10) in 00_config.R):
# 7 historical baseline years + 9 future years x 2 scenarios = 25 predictions.
#
# CLIMATE WINDOWS (this version): each prediction year is a DECADE LABEL; its climate is the
# mean of the yearly rasters in the window W2W_WINDOW_BEFORE .. W2W_WINDOW_AFTER around it
# (default -5/+4: "2010" = 2005-2014, "2020" = 2015-2024, "2100" = 2095-2100), clipped to the
# scenario's years (historical 1950-2014, ssp 2015-2100). The fit used 2001-2020 mean
# climate, so the model is fed mean climate, not single-year weather.
#
# Predictors read from the climate rasters: winTmin and annualPrec (the locked-in
# formula M2_06 does not use meanTemp). Ensemble label comes from
# GCM_ENSEMBLE_LABEL in 00_config.R ("13MME").
# Output file names are unchanged: out/pred_tables/model2_<tag>_{mean,se}.csv with
# tag = baseline_<year> or <scenario>_<year>.
# =============================================================================
if (!exists("DIR_DATA")) source("R/00_config.R")
source("R/w2w_helpers.R")
library(gjam)
stopifnot("gjam did not load as an installed package -- reinstall before running" =
            environmentName(environment(gjam)) == "gjam")

# Model 2 window = 20 years (the fit's 2001-2020 length); winTmin = min over it, annualPrec = mean
if (!exists("W2W_M2_BEFORE")) W2W_M2_BEFORE <- 10
if (!exists("W2W_M2_AFTER"))  W2W_M2_AFTER  <- 9
fit_obj <- readRDS(file.path(DIR_OUT, "fits", "model2_gjam_out.rds"))
fit     <- fit_obj$fit
CENTERS <- fit_obj$centers
xdata_c <- readRDS(file.path(DIR_DATA, "xdata_centered.rds"))

st     <- w2w_load_static()
static <- st$tab
meta   <- st$meta
message(sprintf("static table: %s cells (grid %s%s)", format(nrow(static), big.mark = ","),
                W2W_TAG, W2W_SUFFIX))

# Optional command-line argument: a task number 1..25 (SLURM job array) selects ONE prediction;
# with no argument every prediction runs in turn. Order = baseline years, then each scenario's years.
ARGS     <- commandArgs(trailingOnly = TRUE)
ALL_TAGS <- c(sprintf("baseline_%d", CMIP6_HISTORICAL_YEARS),
              unlist(lapply(SCENARIOS, function(s) sprintf("%s_%d", s, PRED_YEARS_FUTURE))))
ONLY <- NULL
if (length(ARGS) == 1 && grepl("^[0-9]+$", ARGS)) {
  ONLY <- ALL_TAGS[as.integer(ARGS)]
  if (is.na(ONLY)) stop("task number ", ARGS, " is out of range 1..", length(ALL_TAGS))
  message("array task ", ARGS, " of ", length(ALL_TAGS), " -> ", ONLY)
}

# Optional command-line argument: a task number 1..25 (SLURM job array) selects ONE prediction;
# with no argument every prediction runs in turn. Order = baseline years, then each scenario's years.
ARGS     <- commandArgs(trailingOnly = TRUE)
ALL_TAGS <- c(sprintf("baseline_%d", CMIP6_HISTORICAL_YEARS),
              unlist(lapply(SCENARIOS, function(s) sprintf("%s_%d", s, PRED_YEARS_FUTURE))))
ONLY <- NULL
if (length(ARGS) == 1 && grepl("^[0-9]+$", ARGS)) {
  ONLY <- ALL_TAGS[as.integer(ARGS)]
  if (is.na(ONLY)) stop("task number ", ARGS, " is out of range 1..", length(ALL_TAGS))
  message("array task ", ARGS, " of ", length(ALL_TAGS), " -> ", ONLY)
}

predict_one <- function(scenario, yr, tag) {
  if (!is.null(ONLY) && tag != ONLY) return(invisible(NULL))
  if (!is.null(ONLY) && tag != ONLY) return(invisible(NULL))
  stem_out <- sprintf("model2_%s%s", tag, W2W_SUFFIX)
  if (w2w_pred_done(stem_out)) { message("[", tag, "] already predicted -- skipping"); return(invisible(NULL)) }

  # climate = mean of the yearly rasters in this label year's window, clipped to the
  # scenario's own range so a "historical" window never borrows ssp years (or vice versa)
  rng   <- w2w_scenario_range(scenario)
  yrs   <- w2w_window_years(yr, rng[1], rng[2], W2W_M2_BEFORE, W2W_M2_AFTER)
  stems <- vapply(yrs, function(y) w2w_stem_m2(scenario, y), "")
  message(sprintf("[%s] climate = mean of %d-%d (%d years)", tag, min(yrs), max(yrs), length(yrs)))
  clim <- w2w_climate_window(stems, M2_BANDS, static$geeID, meta, min_bands = "winTmin")
  if (is.null(clim)) {
    message("[", tag, "] skipping -- export the missing years (python gee/wall_to_wall_climate_m2.py --scenarios ",
            scenario, " --years ", paste(yrs, collapse = " "), ") and download them to ", W2W_DIR)
    return(invisible(NULL))
  }

  grid <- cbind(static, clim)
  need <- c("north_component", "east_component", "cec30", "sand30", "winTmin", "annualPrec")
  ok <- complete.cases(as.data.frame(grid)[, need])
  if (!all(ok)) message(sprintf("[%s] dropping %d of %d cells with a missing predictor",
                                tag, sum(!ok), length(ok)))
  grid <- grid[ok]
  message(sprintf("[%s] predicting %s cells", tag, format(nrow(grid), big.mark = ",")))

  newx <- data.frame(
    northComponent = grid$north_component,
    eastComponent  = grid$east_component,
    cec30          = grid$cec30,
    sand30         = grid$sand30,
    clay30         = grid$clay30,
    winTmin        = grid$winTmin,
    annualPrec     = grid$annualPrec
  )
  w2w_range_check(newx, xdata_c, sprintf("[%s]", tag))
  for (v in names(CENTERS)) if (v %in% names(newx)) newx[[paste0(v, "Ctr")]] <- newx[[v]] - CENTERS[[v]]

  pred <- w2w_predict_chunked(fit, newx)
  w2w_write_pred(grid$geeID, pred, stem_out)
  message(sprintf("[%s] wrote out/pred_tables/%s_{mean,se}.csv", tag, stem_out))
}

## ---- historical baseline (CMIP6 historical scenario, capped at 2014) ---------
for (yr in CMIP6_HISTORICAL_YEARS) predict_one("historical", yr, sprintf("baseline_%d", yr))

## ---- future scenarios -----------------------------------------------------------
for (scen in SCENARIOS) for (yr in PRED_YEARS_FUTURE) predict_one(scen, yr, sprintf("%s_%d", scen, yr))

message("\nModel 2 wall-to-wall prediction finished for the available year/scenario combinations.")

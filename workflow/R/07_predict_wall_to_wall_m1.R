# =============================================================================
# 07_predict_wall_to_wall_m1.R -- Model 1's historical-OBSERVED wall-to-wall
# prediction. For each year in MODEL1_WALLTOWALL_YEARS: read that year's
# remote-sensing GeoTIFF, join it to the forest-cell static table by cell
# number, predict with gjamPredict, write mean/SE tables.
#
# CLIMATE WINDOWS (this version): each prediction year is a DECADE LABEL, and its
#   climate is the MEAN of the yearly rasters in a window around it
#   (W2W_WINDOW_BEFORE / W2W_WINDOW_AFTER, default -5/+4: "2010" = 2005-2014; the last
#   label is clipped to the observed record, so "2020" = 2015-2020). The model was fit on
#   2001-2020 mean climate, so it is fed mean climate here, not single-year weather.
#   Export EVERY year in the windows (see w2w_missing()).
#
# CHANGES in this version (raster pipeline):
#   - Inputs are GeoTIFFs on one fixed grid (gee/wall_to_wall_*.py), joined by
#     cell number, not CSVs merged on a string geeID.
#   - Every cell in the static table is predicted: whole-US land by default
#     (LAND_FRAC_MIN); set FOREST_FRAC_MIN > 0 in 00_config.R to restrict to forest.
#   - diffLSTairT is computed here as LST - airT_TC (Surface minus Air, the
#     locked-in convention) instead of being exported.
#   - Prediction runs in chunks of PRED_CHUNK cells with a progress/ETA line, and
#     a year whose output already exists is skipped, so an interrupted run
#     resumes where it stopped (set W2W_OVERWRITE <- TRUE to redo).
#   - Earlier fixes kept: library(gjam) (not the sourceCpp() hack), gjam-safe
#     column names (northComponent, eastComponent, diffLSTairT), "Ctr" suffix
#     for centred terms, centres taken from the FIT, never recomputed here.
# Output file names are unchanged: out/pred_tables/model1_historical_<year>_{mean,se}.csv
# =============================================================================
if (!exists("DIR_DATA")) source("R/00_config.R")
source("R/w2w_helpers.R")
library(gjam)
stopifnot("gjam did not load as an installed package -- reinstall before running" =
            environmentName(environment(gjam)) == "gjam")

fit_obj <- readRDS(file.path(DIR_OUT, "fits", "model1_gjam_out.rds"))
fit     <- fit_obj$fit
CENTERS <- fit_obj$centers
xdata_c <- readRDS(file.path(DIR_DATA, "xdata_centered.rds"))   # fitted ranges, for the extrapolation check

st     <- w2w_load_static()
static <- st$tab
meta   <- st$meta
message(sprintf("static table: %s cells (grid %s%s)", format(nrow(static), big.mark = ","),
                W2W_TAG, W2W_SUFFIX))

YEARS <- MODEL1_WALLTOWALL_YEARS   # decades since 2001 up to 2014 (= 2010); _EXT adds 2020

for (yr in YEARS) {
  stem_out <- sprintf("model1_historical_%d%s", yr, W2W_SUFFIX)
  if (w2w_pred_done(stem_out)) { message("[", yr, "] already predicted -- skipping"); next }

  # climate = mean of the yearly rasters in this label year's window (e.g. 2010 -> 2005-2014)
  yrs   <- w2w_window_years(yr, min(w2w_m1_range()), max(w2w_m1_range()))
  stems <- vapply(yrs, w2w_stem_m1, "")
  message(sprintf("[%d] climate = mean of %d-%d (%d years)", yr, min(yrs), max(yrs), length(yrs)))
  clim <- w2w_climate_window(stems, M1_BANDS, static$geeID, meta)
  if (is.null(clim)) {
    message("[", yr, "] skipping -- export the missing years (python gee/wall_to_wall_climate_m1.py --years ",
            paste(yrs, collapse = " "), ") and download them to ", W2W_DIR)
    next
  }

  grid <- cbind(static, clim)
  grid[, diffLST_airT := LST - airT_TC]            # Surface minus Air
  need <- c("north_component", "east_component", "cec30", "clay30",
            "precipGPM", "LST", "diffLST_airT")
  ok <- complete.cases(as.data.frame(grid)[, need])
  if (!all(ok)) message(sprintf("[%d] dropping %d of %d cells with a missing predictor",
                                yr, sum(!ok), length(ok)))
  grid <- grid[ok]
  message(sprintf("[%d] predicting %s cells", yr, format(nrow(grid), big.mark = ",")))

  newx <- data.frame(
    northComponent = grid$north_component,
    eastComponent  = grid$east_component,
    cec30          = grid$cec30,
    sand30         = grid$sand30,
    clay30         = grid$clay30,
    precipGPM      = grid$precipGPM,
    LST            = grid$LST,
    diffLSTairT    = grid$diffLST_airT
  )
  w2w_range_check(newx, xdata_c, sprintf("[%d]", yr))

  # centred terms use the SAME centres as the fit -- never recompute on this grid
  for (v in names(CENTERS)) if (v %in% names(newx)) newx[[paste0(v, "Ctr")]] <- newx[[v]] - CENTERS[[v]]

  pred <- w2w_predict_chunked(fit, newx)
  w2w_write_pred(grid$geeID, pred, stem_out)
  message(sprintf("[%d] wrote out/pred_tables/%s_{mean,se}.csv", yr, stem_out))
}

message("\nModel 1 wall-to-wall prediction finished for the available years.")

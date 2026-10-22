# =============================================================================
# 07_predict_1km_tiled.R -- 1 km wall-to-wall prediction, ONE TILE of ONE prediction per run.
#
#   Rscript R/07_predict_1km_tiled.R list            # print the 26 predictions
#   Rscript R/07_predict_1km_tiled.R <pred 1..26> <tile 1..NTILES>
#
# Same model, formula, centres, climate windows and predictors as 07_predict_wall_to_wall_m1/m2.R
# (4 km); only the grid is 1 km and the cells are split into NTILES blocks (env W2W_NTILES,
# default 30) so each SLURM array task does ~1/30 of the US. A finished tile is saved as
# out/pred_tiles_1km/<stem>_tNNN_ofNNN.rds and skipped on a re-run (resume after a timeout).
# Does NOT read or write any 4 km output.
#
# Environment variables:
#   W2W_NTILES  number of tiles (default 30; must be the same for every task and for 07c)
#   W2W_NSIM    gjamPredict nsim (default NSIM_PRED from 00_config.R)
#   W2W_SAMPLE  if > 0: BENCHMARK MODE -- predict this many random cells (whole US) and save to
#               out/bench_1km/<stem>_nsim<N>.rds instead of a tile; used to choose nsim.
# =============================================================================
if (!exists("DIR_DATA")) source("R/00_config.R")
GRID_RES_M <- 1000L                      # this script only; the shared config file is never edited
source("R/w2w_helpers.R")
stopifnot("W2W_TAG is not 1km -- GRID_RES_M override failed" = W2W_TAG == "1km")
source("R/w2w_1km_specs.R")
SPECS <- w2w_specs()

ARGS <- commandArgs(trailingOnly = TRUE)
if (length(ARGS) == 1 && ARGS == "list") { print(SPECS[, c("idx", "model", "scen", "yr", "stem")]); quit(save = "no") }
if (length(ARGS) != 2 || !all(grepl("^[0-9]+$", ARGS)))
  stop("usage: Rscript R/07_predict_1km_tiled.R <pred 1..", nrow(SPECS), "> <tile>   (or: list)")
PI   <- as.integer(ARGS[1]); TILE <- as.integer(ARGS[2])
NT   <- as.integer(Sys.getenv("W2W_NTILES", "30"))
NSIM <- as.integer(Sys.getenv("W2W_NSIM", as.character(NSIM_PRED)))
SAMPLE_N <- as.integer(Sys.getenv("W2W_SAMPLE", "0"))
if (is.na(PI) || PI < 1 || PI > nrow(SPECS)) stop("pred must be 1..", nrow(SPECS))
if (TILE < 1 || TILE > NT) stop("tile must be 1..", NT, " (W2W_NTILES)")
sp <- SPECS[PI, ]
message(sprintf("1 km | pred %d = %s (model %d, %s %d) | tile %d of %d | nsim %d%s",
                PI, sp$stem, sp$model, sp$scen, sp$yr, TILE, NT, NSIM,
                if (SAMPLE_N > 0) sprintf(" | BENCHMARK sample %d", SAMPLE_N) else ""))

dir.create(file.path(DIR_OUT, "pred_tiles_1km"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_OUT, "bench_1km"),      showWarnings = FALSE, recursive = TRUE)
f_tile  <- w2w_tile_file(sp$stem, TILE, NT)
f_bench <- file.path(DIR_OUT, "bench_1km", sprintf("%s_nsim%d.rds", sp$stem, NSIM))
if (SAMPLE_N == 0 && file.exists(f_tile) && !isTRUE(W2W_OVERWRITE)) { message("tile already done -- skipping"); quit(save = "no") }

library(gjam)
stopifnot("gjam did not load as an installed package" =
            environmentName(environment(gjam)) == "gjam")
fit_obj <- readRDS(file.path(DIR_OUT, "fits", sprintf("model%d_gjam_out.rds", sp$model)))
fit     <- fit_obj$fit
CENTERS <- fit_obj$centers
xdata_c <- readRDS(file.path(DIR_DATA, "xdata_centered.rds"))

st     <- w2w_load_static()
static <- st$tab
meta   <- st$meta
n <- nrow(static)
message(sprintf("static table: %s land cells (grid %s)", format(n, big.mark = ","), W2W_TAG))
if (SAMPLE_N > 0) {                                    # benchmark: same random cells every run
  set.seed(1); rows <- sort(sample(n, min(SAMPLE_N, n)))
} else {                                               # tile = contiguous block of the cell-number-ordered table
  b <- floor(seq(0, n, length.out = NT + 1)); rows <- (b[TILE] + 1):b[TILE + 1]
}
static <- static[rows]
message(sprintf("this task: %s cells", format(nrow(static), big.mark = ",")))

if (sp$model == 1L) {
  yrs   <- w2w_window_years(sp$yr, min(w2w_m1_range()), max(w2w_m1_range()))
  stems <- vapply(yrs, w2w_stem_m1, "")
  clim  <- w2w_climate_window(stems, M1_BANDS, static$geeID, meta)
} else {
  rng   <- w2w_scenario_range(sp$scen)
  yrs   <- w2w_window_years(sp$yr, rng[1], rng[2], 10, 9)       # 20-year window, as at 4 km
  stems <- vapply(yrs, function(y) w2w_stem_m2(sp$scen, y), "")
  clim  <- w2w_climate_window(stems, M2_BANDS, static$geeID, meta, min_bands = "winTmin")
}
message(sprintf("climate window %d-%d (%d years)", min(yrs), max(yrs), length(yrs)))
if (is.null(clim)) stop("climate window incomplete -- download the missing 1 km files first (see w2w_missing())")

grid <- cbind(static, clim)
if (sp$model == 1L) {
  grid[, diffLST_airT := LST - airT_TC]                       # Surface minus Air (locked-in)
  need <- c("north_component", "east_component", "cec30", "clay30", "precipGPM", "LST", "diffLST_airT")
} else {
  need <- c("north_component", "east_component", "cec30", "sand30", "winTmin", "annualPrec")
}
ok <- complete.cases(as.data.frame(grid)[, need])
if (!all(ok)) message(sprintf("dropping %d of %d cells with a missing predictor", sum(!ok), length(ok)))
grid <- grid[ok]

if (sp$model == 1L) {
  newx <- data.frame(northComponent = grid$north_component, eastComponent = grid$east_component,
                     cec30 = grid$cec30, sand30 = grid$sand30, clay30 = grid$clay30,
                     precipGPM = grid$precipGPM, LST = grid$LST, diffLSTairT = grid$diffLST_airT)
} else {
  newx <- data.frame(northComponent = grid$north_component, eastComponent = grid$east_component,
                     cec30 = grid$cec30, sand30 = grid$sand30, clay30 = grid$clay30,
                     winTmin = grid$winTmin, annualPrec = grid$annualPrec)
}
w2w_range_check(newx, xdata_c, "[1km]")
for (v in names(CENTERS)) if (v %in% names(newx)) newx[[paste0(v, "Ctr")]] <- newx[[v]] - CENTERS[[v]]

t0   <- Sys.time()
pred <- w2w_predict_chunked(fit, newx, nsim = NSIM)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(sprintf("predicted %s cells in %.1f min (%.1f ms/cell)", format(nrow(newx), big.mark = ","),
                secs / 60, 1000 * secs / nrow(newx)))

res <- list(geeID = grid$geeID, mean = pred$mean, se = pred$se, nsim = NSIM,
            stem = sp$stem, secs_per_cell = secs / nrow(newx))
f_out <- if (SAMPLE_N > 0) f_bench else f_tile
saveRDS(res, paste0(f_out, ".tmp")); file.rename(paste0(f_out, ".tmp"), f_out)
message("wrote ", f_out)

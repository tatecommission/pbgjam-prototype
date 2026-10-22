# =============================================================================
# w2w_helpers.R -- shared helpers for the raster wall-to-wall pipeline.
# Source this after 00_config.R. Mirrors gee/wall_to_wall_common.py.
#
# THE IDEA: every GEE export is a GeoTIFF on ONE fixed EPSG:5070 grid, so a
# cell's number (terra's row-major cell index) is the same in every file.
# That cell number is the `geeID` column the downstream scripts already use.
# =============================================================================
if (!exists("DIR_DATA")) source("R/00_config.R")
suppressPackageStartupMessages({ library(terra); library(data.table) })

W2W_NODATA   <- -9999                               # sentinel written by the GEE scripts
W2W_TAG      <- sprintf("%dkm", as.integer(GRID_RES_M / 1000))
W2W_DIR      <- file.path(DIR_DATA, "gee_exports")
W2W_SUFFIX   <- if (exists("W2W_TEST") && isTRUE(W2W_TEST)) "_test" else ""
# Domain: whole-US land by default. FOREST_FRAC_MIN = 0 means "no forest filter";
# raise it (e.g. 0.20) to restrict to forest cells again. LAND_FRAC_MIN keeps cells
# that are at least this much NLCD land (0.25 keeps thin coastal strips such as
# the Outer Banks while dropping open water and anything outside the US).
if (!exists("FOREST_FRAC_MIN")) FOREST_FRAC_MIN <- 0
if (!exists("LAND_FRAC_MIN"))   LAND_FRAC_MIN   <- 0.25
if (!exists("PRED_CHUNK"))      PRED_CHUNK      <- 20000L
if (!exists("W2W_OVERWRITE"))   W2W_OVERWRITE   <- FALSE
# Climate window for each prediction (label) year: the label year's window is
# (label - W2W_WINDOW_BEFORE) .. (label + W2W_WINDOW_AFTER), clipped to the years the
# data / scenario cover. 5 and 4 give "2010" = 2005-2014 and "2020" = 2015-2024, so no
# window straddles the CMIP6 historical (<= 2014) / scenario (>= 2015) boundary.
# 0 and 0 reproduces single-year snapshots.
if (!exists("W2W_WINDOW_BEFORE")) W2W_WINDOW_BEFORE <- 5
if (!exists("W2W_WINDOW_AFTER"))  W2W_WINDOW_AFTER  <- 4

STATIC_BANDS <- c("north_component", "east_component", "cec30", "clay30", "sand30",
                  "forest_frac", "land_frac")
M1_BANDS     <- c("precipGPM", "LST", "airT_TC")
M2_BANDS     <- c("winTmin", "annualPrec")

## ---- file-name stems (must match gee/wall_to_wall_common.py) ----------------
w2w_stem_static <- function(suffix = W2W_SUFFIX) sprintf("gee_static_conus_%s%s", W2W_TAG, suffix)
w2w_stem_m1     <- function(year, suffix = W2W_SUFFIX) sprintf("gee_climate_m1_%s_%d%s", W2W_TAG, year, suffix)
w2w_stem_m2     <- function(scenario, year, suffix = W2W_SUFFIX)
  sprintf("gee_climate_m2_%s_%s_%s_%d%s", GCM_ENSEMBLE_LABEL, W2W_TAG, scenario, year, suffix)

w2w_tif_files <- function(stem) {
  # a large export can be split by GEE into name-0000000000-0000000000.tif tiles
  list.files(W2W_DIR, pattern = paste0("^", stem, "(-[0-9]+-[0-9]+)?\\.tif$"), full.names = TRUE)
}

## ---- read one export as a SpatRaster with NODATA -> NA and known band names --
w2w_read <- function(stem, expected) {
  f <- w2w_tif_files(stem)
  if (length(f) == 0) return(NULL)
  r <- if (length(f) == 1) rast(f) else vrt(sort(f))
  nm <- names(r)
  if (all(expected %in% nm)) {
    r <- r[[expected]]
  } else {
    if (nlyr(r) != length(expected))
      stop(stem, ".tif has ", nlyr(r), " bands but ", length(expected), " were expected (",
           paste(expected, collapse = ", "), ") -- was it exported with different --vars?")
    names(r) <- expected          # GEE writes bands in the order the script lists them
  }
  subst(r, W2W_NODATA, NA)
}

w2w_meta <- function(r) list(nrow = nrow(r), ncol = ncol(r), ext = as.vector(ext(r)))

## ---- static GeoTIFF -> table of land cells (same columns as the old CSV, plus forest_frac, land_frac) ---
build_static_table <- function(suffix = W2W_SUFFIX) {
  stem <- w2w_stem_static(suffix)
  r <- w2w_read(stem, STATIC_BANDS)
  if (is.null(r)) stop("not found: ", file.path(W2W_DIR, paste0(stem, ".tif")),
                       " -- run gee/wall_to_wall_static.py and download its output first")
  v <- values(r, mat = TRUE)
  # unit guard: clay30 and sand30 are percentages, so values above 100 mean this raster
  # holds SoilGrids' native units (10x too large) from before the unit conversion.
  if (any(v[, c("clay30", "sand30")] > 100, na.rm = TRUE))
    stop(stem, ".tif has clay30/sand30 above 100, which a percentage cannot be -- it holds ",
         "SoilGrids' native units (10x too large). Re-export with the current ",
         "gee/wall_to_wall_static.py (add --force), delete the old file from Drive, and re-download.")
  is_land <- !is.na(v[, "land_frac"]) & v[, "land_frac"] >= LAND_FRAC_MIN
  is_dom  <- is_land & !is.na(v[, "forest_frac"]) & v[, "forest_frac"] >= FOREST_FRAC_MIN
  keep <- which(is_dom & complete.cases(v))
  message(sprintf(paste0("grid %d x %d = %s cells | %s have land_frac >= %.2f | %s also have forest_frac >= %.2f | ",
                         "%s of those have complete topography + soil"),
                  ncol(r), nrow(r), format(ncell(r), big.mark = ","),
                  format(sum(is_land), big.mark = ","), LAND_FRAC_MIN,
                  format(sum(is_dom), big.mark = ","), FOREST_FRAC_MIN,
                  format(length(keep), big.mark = ",")))
  if (length(keep) == 0) stop("no cells survived the land/forest filter -- check the static raster")
  xy  <- xyFromCell(r, keep)
  ll  <- crds(project(vect(xy, type = "points", crs = crs(r)), "EPSG:4326"))
  tab <- cbind(data.table(geeID = keep, lon = ll[, 1], lat = ll[, 2]),
               as.data.table(v[keep, , drop = FALSE]))
  fwrite(tab, file.path(W2W_DIR, paste0(stem, ".csv")))
  saveRDS(w2w_meta(r), file.path(W2W_DIR, paste0(stem, "_meta.rds")))
  invisible(tab)
}

w2w_load_static <- function(suffix = W2W_SUFFIX) {
  stem <- w2w_stem_static(suffix)
  csv  <- file.path(W2W_DIR, paste0(stem, ".csv"))
  mf   <- file.path(W2W_DIR, paste0(stem, "_meta.rds"))
  if (!file.exists(csv) || !file.exists(mf))
    stop("static table not built -- run R/07a_build_static_table.R first (needs ", stem, ".tif)")
  list(tab = fread(csv), meta = readRDS(mf))
}

## ---- climate raster(s) -> values at the kept cells, coastal gaps filled -------
# Land-only climate products (TerraClimate, GCM output) have no data a few km
# offshore, so coastal cells can come back NA. Each pass copies the mean of the
# 3x3 neighbourhood into NA cells (one cell = res_m of reach).
w2w_check_grid <- function(r, meta, stem) {
  same <- nrow(r) == meta$nrow && ncol(r) == meta$ncol &&
    isTRUE(all.equal(as.vector(ext(r)), meta$ext))
  if (!same) stop(stem, ".tif is not on the same grid as the static raster -- ",
                  "re-export both with the same --res (and both with or without --test)")
  invisible(TRUE)
}

w2w_fill_cells <- function(r, ids, label, max_fill_iter = 10) {
  v <- values(r, mat = TRUE)[ids, , drop = FALSE]
  n0 <- sum(!complete.cases(v)); it <- 0
  while (any(!complete.cases(v)) && it < max_fill_iter) {
    r  <- focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
    v  <- values(r, mat = TRUE)[ids, , drop = FALSE]
    it <- it + 1
  }
  if (n0 > 0)
    message(sprintf("  %s: %d of %d cells had no climate value; %d still missing after %d fill passes",
                    label, n0, length(ids), sum(!complete.cases(v)), it))
  as.data.table(v)
}

# one year (kept for single-year snapshots; the 07_* scripts now use the window reader)
w2w_climate_cells <- function(stem, expected, ids, meta, max_fill_iter = 10) {
  r <- w2w_read(stem, expected)
  if (is.null(r)) return(NULL)
  w2w_check_grid(r, meta, stem)
  w2w_fill_cells(r, ids, stem, max_fill_iter)
}

## ---- climate WINDOWS: mean of several yearly rasters, then gap fill --------------
# Years a prediction label covers, clipped to [lo, hi] (the data / scenario range).
w2w_window_years <- function(label, lo, hi, before = W2W_WINDOW_BEFORE, after = W2W_WINDOW_AFTER) {
  a <- max(lo, label - before); b <- min(hi, label + after)
  if (a > b) integer(0) else seq(a, b)
}
w2w_scenario_range <- function(scenario) if (scenario == "historical") c(1950, 2014) else c(2015, 2100)
w2w_m1_range       <- function() range(FIT_YEARS)       # observed data used for the fit: 2001-2020

# "1950-1960, 1962" style summary of a set of years
w2w_runs <- function(y) {
  if (!length(y)) return("none")
  y  <- sort(unique(y)); br <- c(0, which(diff(y) != 1), length(y))
  paste(vapply(seq_len(length(br) - 1), function(i) {
    a <- y[br[i] + 1]; b <- y[br[i + 1]]
    if (a == b) as.character(a) else paste0(a, "-", b)
  }, ""), collapse = ", ")
}

# Mean over the yearly rasters in `stems` at each cell (ignoring NA years), then
# coastal-gap fill. Returns NULL (and says which files are missing) if any year
# of the window has not been downloaded, so a window is never silently shorter.
w2w_climate_window <- function(stems, expected, ids, meta, max_fill_iter = 10, min_bands = character()) {
  have <- vapply(stems, function(s) length(w2w_tif_files(s)) > 0, TRUE)
  if (!all(have)) {
    message("  window incomplete -- missing: ", paste(stems[!have], collapse = ", "))
    return(NULL)
  }
  sums <- NULL; cnt <- NULL; template <- NULL; mins <- NULL
  for (s in stems) {
    r <- w2w_read(s, expected)
    w2w_check_grid(r, meta, s)
    v  <- values(r, mat = TRUE)
    ok <- !is.na(v)
    vm <- v; vm[!ok] <- Inf                       # for min bands: NA never wins
    mins <- if (is.null(mins)) vm else pmin(mins, vm)
    v[!ok] <- 0
    if (is.null(sums)) { sums <- v; cnt <- ok + 0; template <- r } else { sums <- sums + v; cnt <- cnt + ok }
  }
  avg <- sums / cnt; avg[cnt == 0] <- NA
  # bands in min_bands (winTmin) take the MINIMUM across the window years, not the mean:
  # the fit's winTmin is the coldest monthly value over 2001-2020, an extreme, not an average
  mins[is.infinite(mins)] <- NA
  mb <- intersect(min_bands, colnames(avg))
  if (length(mb)) avg[, mb] <- mins[, mb]
  values(template) <- avg
  w2w_fill_cells(template, ids, sprintf("%d-year window", length(stems)), max_fill_iter)
}

## ---- diagnostic: how much of the grid is outside the fitted predictor range ---
# Catches unit mismatches (e.g. SoilGrids mapped units vs the fit-side soil
# columns) and shows how much of each map is extrapolation.
w2w_range_check <- function(newx, ref, label = "") {
  cols <- intersect(names(newx), names(ref))
  frac <- sapply(cols, function(cn) {
    rng <- range(ref[[cn]], na.rm = TRUE)
    mean(newx[[cn]] < rng[1] | newx[[cn]] > rng[2], na.rm = TRUE)
  })
  message(label, " cells outside the fitted range -- ",
          paste(sprintf("%s %.1f%%", cols, 100 * frac), collapse = "; "))
  invisible(frac)
}

## ---- gjamPredict workaround for models with squared terms ------------------------
# gjam's .getStandX looks up every name in fit$inputs$standRows as a COLUMN of the new
# data, including design terms like "I(precipGPMCtr^2)", which only exist inside the model
# matrix. Without them gjamPredict stops with "undefined columns selected". Add them
# (value = square of the raw column; gjam rebuilds the real standardised square itself).
# Also refuse factor columns: with any factor present .getStandX writes the standardised
# values into the wrong columns (it indexes by position) and predictions are nonsense.
# Checked against the fit's stored design matrix: difference 0.
w2w_add_poly_dummies <- function(fit, newx) {
  if (any(vapply(newx, is.factor, TRUE)))
    stop("newx has a factor column; gjamPredict would mis-standardise (drop it)")
  for (nm in setdiff(names(fit$inputs$standRows), names(newx))) {
    m <- regmatches(nm, regexec("^I\\((.+)\\^2\\)$", nm))[[1]]
    if (length(m) != 2 || !(m[2] %in% names(newx)))
      stop("the fit expects column '", nm, "' and it cannot be built from newx")
    newx[[nm]] <- newx[[m[2]]]^2
  }
  newx
}

## ---- gjamPredict in chunks with progress/ETA ---------------------------------
w2w_predict_chunked <- function(fit, newx, nsim = NSIM_PRED, chunk = PRED_CHUNK) {
  newx <- w2w_add_poly_dummies(fit, newx)
  n <- nrow(newx)
  starts <- seq(1, n, by = chunk)
  means <- vector("list", length(starts)); ses <- vector("list", length(starts))
  t0 <- Sys.time()
  for (i in seq_along(starts)) {
    idx  <- starts[i]:min(starts[i] + chunk - 1, n)
    pred <- gjamPredict(fit, newdata = list(xdata = newx[idx, , drop = FALSE], nsim = nsim))
    if (!all(c("yMu", "yPe") %in% names(pred$sdList)))
      stop("unexpected gjamPredict output; sdList has: ", paste(names(pred$sdList), collapse = ", "))
    means[[i]] <- pred$sdList$yMu
    if (i == 1) {                       # sanity: the fitted data average ~20 total basal area
      tot <- mean(rowSums(pred$sdList$yMu))
      message(sprintf("  first chunk: mean total basal area %.2f (fit data in-sample ~21)", tot))
      if (!is.finite(tot) || tot > 500) stop("implausible predictions (mean total ", tot, ") -- stopping")
    }
    ses[[i]]   <- pred$sdList$yPe
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  chunk %d/%d (%d cells) | %.1f min elapsed | ~%.1f min remaining",
                    i, length(starts), length(idx), el, el / i * (length(starts) - i)))
  }
  list(mean = do.call(rbind, means), se = do.call(rbind, ses))
}

## ---- write mean/se tables (same file layout as before; atomic rename) ----------
w2w_write_pred <- function(ids, pred, stem_out) {
  out <- data.frame(geeID = ids, pred$mean)
  names(out)[-1] <- paste0(names(out)[-1], "_mean")
  se  <- data.frame(geeID = ids, pred$se)
  names(se)[-1]  <- paste0(names(se)[-1], "_se")
  dir.create(file.path(DIR_OUT, "pred_tables"), showWarnings = FALSE, recursive = TRUE)
  f_se   <- file.path(DIR_OUT, "pred_tables", paste0(stem_out, "_se.csv"))
  f_mean <- file.path(DIR_OUT, "pred_tables", paste0(stem_out, "_mean.csv"))
  fwrite(se, paste0(f_se, ".tmp"));     file.rename(paste0(f_se, ".tmp"), f_se)
  fwrite(out, paste0(f_mean, ".tmp"));  file.rename(paste0(f_mean, ".tmp"), f_mean)   # mean last = "done" marker
}
w2w_pred_done <- function(stem_out)
  file.exists(file.path(DIR_OUT, "pred_tables", paste0(stem_out, "_mean.csv"))) && !isTRUE(W2W_OVERWRITE)

## ---- which yearly GeoTIFFs are still missing from data/gee_exports? -----------------
# Every year inside every prediction window is needed. Prints, per series, what is
# missing and the command that exports it.
w2w_missing <- function(suffix = W2W_SUFFIX) {
  has <- function(stem) length(w2w_tif_files(stem)) > 0
  need <- function(labels, lo, hi) sort(unique(unlist(lapply(labels, w2w_window_years, lo = lo, hi = hi))))
  series <- list(
    list(name = "static", stems = w2w_stem_static(suffix), years = NA, cmd = "python gee/wall_to_wall_static.py"))
  y1 <- need(MODEL1_WALLTOWALL_YEARS, min(FIT_YEARS), max(FIT_YEARS))
  series[[2]] <- list(name = "M1 observed", years = y1, stems = vapply(y1, w2w_stem_m1, "", suffix = suffix),
                      cmd = "python gee/wall_to_wall_climate_m1.py --years")
  yh <- need(CMIP6_HISTORICAL_YEARS, 1950, 2014)
  series[[3]] <- list(name = "M2 historical", years = yh,
                      stems = vapply(yh, function(y) w2w_stem_m2("historical", y, suffix), ""),
                      cmd = "python gee/wall_to_wall_climate_m2.py --scenarios historical --years")
  yf <- need(PRED_YEARS_FUTURE, 2015, 2100)
  for (sc in SCENARIOS)
    series[[length(series) + 1]] <- list(name = paste("M2", sc), years = yf,
        stems = vapply(yf, function(y) w2w_stem_m2(sc, y, suffix), ""),
        cmd = paste("python gee/wall_to_wall_climate_m2.py --scenarios", sc, "--years"))
  total_missing <- 0
  for (sr in series) {
    miss <- !vapply(sr$stems, has, TRUE)
    total_missing <- total_missing + sum(miss)
    has_years <- !is.na(sr$years[1])
    message(sprintf("%-14s %3d of %3d present%s", sr$name, sum(!miss), length(miss),
                    if (any(miss) && has_years) paste0("  | missing years: ", w2w_runs(sr$years[miss])) else ""))
    if (any(miss))
      message("    ", sr$cmd, if (has_years) paste0(" ", paste(sr$years[miss], collapse = " ")) else "")
  }
  invisible(total_missing)
}

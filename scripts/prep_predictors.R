# prep_predictors.R -- run ONCE (RStudio: Source). Makes the predictor maps shown in the app.
# Uses the 4 km files already in v4_gjam_workflow/data/gee_exports.
#   Climate (5 periods x 2 variables): label L uses years (L - WIN_BEFORE) .. (L + WIN_AFTER)
#       winTmin    = coldest value over the window (as in the model code)
#       annualPrec = mean over the window
#   Static (6 files): soil CEC, soil sand, and the two slope-aspect components, masked to land.
# Output: maps_v4/predictors/
library(terra)

# !! UNVERIFIED: must equal W2W_WINDOW_BEFORE / W2W_WINDOW_AFTER in the 00_config.R the 1 km run used on the DCC.
#    Check on the DCC:  grep -n W2W_WINDOW /work/tvc11/PBGJAM/v4_gjam_workflow/R/00_config.R
#    5 / 4 = 10-year windows.  10 / 9 = 20-year windows.  Keep app.R's BA_W_BEFORE / BA_W_AFTER the same.
WIN_BEFORE <- 5
WIN_AFTER  <- 4

gx  <- path.expand("~/Documents/PBGJAM-data-explorer/v4_gjam_workflow/data/gee_exports")
out <- path.expand("~/Documents/PBGJAM-data-explorer/maps_v4/predictors")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

sets <- list(baseline_2000 = list(scen = "historical", label = 2000),
             ssp245_2060   = list(scen = "ssp245",     label = 2060),
             ssp245_2090   = list(scen = "ssp245",     label = 2090),
             ssp585_2060   = list(scen = "ssp585",     label = 2060),
             ssp585_2090   = list(scen = "ssp585",     label = 2090))
for (tag in names(sets)) {
  st    <- sets[[tag]]
  yrs   <- (st$label - WIN_BEFORE):(st$label + WIN_AFTER)
  files <- file.path(gx, sprintf("gee_climate_m2_13MME_4km_%s_%d.tif", st$scen, yrs))
  miss  <- files[!file.exists(files)]
  if (length(miss)) stop(tag, ": missing ", paste(basename(miss), collapse = ", "))
  rs <- lapply(files, rast)
  wt <- app(rast(lapply(rs, function(r) r[["winTmin"]])),    "min",  na.rm = TRUE)
  ap <- app(rast(lapply(rs, function(r) r[["annualPrec"]])), "mean", na.rm = TRUE)
  names(wt) <- "winTmin"; names(ap) <- "annualPrec"
  writeRaster(wt, file.path(out, paste0("winTmin_",    tag, ".tif")), overwrite = TRUE)
  writeRaster(ap, file.path(out, paste0("annualPrec_", tag, ".tif")), overwrite = TRUE)
  message(tag, ": years ", min(yrs), "-", max(yrs), " done")
}

# ---- static: soil + topography ------------------------------------------------------------
STATIC_BANDS <- c("north_component", "east_component", "cec30", "clay30", "sand30", "forest_frac", "land_frac")
s <- rast(file.path(gx, "gee_static_conus_4km.tif"))
if (!all(STATIC_BANDS %in% names(s))) names(s) <- STATIC_BANDS
land <- s[["land_frac"]] >= 0.25                       # same land rule as the predictions
for (b in c("cec30", "sand30", "north_component", "east_component")) {
  x <- mask(s[[b]], land, maskvalues = FALSE)
  x[x < -9000] <- NA                                   # -9999 = 'no soil data' code, not a real value
  writeRaster(x, file.path(out, paste0("static_", b, ".tif")), overwrite = TRUE)
  g <- global(x, c("min", "max"), na.rm = TRUE)
  message("static ", b, ": range ", signif(g[1, 1], 4), " to ", signif(g[1, 2], 4))
}
message("Done. Check the ranges above: sand should be about 0-100 (%), CEC roughly 0-60 (cmol(c)/kg),",
        " slope indices between -1 and 1. Files in ", out)

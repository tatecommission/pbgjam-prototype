# =============================================================================
# 07a_build_static_table.R -- turn the static GeoTIFF (gee/wall_to_wall_static.py)
# into the table of cells that the 07_* prediction scripts and
# 08_maps_per_species.R read:
#   data/gee_exports/gee_static_conus_<res>km.csv   (geeID, lon, lat, north_component,
#                                                    east_component, cec30, clay30, sand30, forest_frac)
#   data/gee_exports/gee_static_conus_<res>km_meta.rds   (grid geometry, used to check
#                                                    that every climate raster matches)
#
# geeID is the terra cell number on the fixed grid -- identical in every raster.
# Cells are kept if land_frac >= LAND_FRAC_MIN (whole-US land), forest_frac >=
# FOREST_FRAC_MIN (0 = no forest filter) and every static layer is present.
# forest_frac stays in the table so maps can be masked to forest afterwards.
#
# Set W2W_TEST <- TRUE in 00_config.R to build from the *_test.tif North Carolina
# window instead (do this first).
# =============================================================================
if (!exists("DIR_DATA")) source("R/00_config.R")
source("R/w2w_helpers.R")

tab <- build_static_table()

cat("\n--- static layers at the forest cells ---\n")
print(summary(tab[, .(north_component, east_component, cec30, clay30, sand30, forest_frac, land_frac)]))

## ---- compare against what the model was FITTED on ---------------------------
# SoilGrids on Earth Engine is in "mapped units" (e.g. clay/sand in g/kg), and the
# fit-side columns came from the lab's own extraction. If the medians below are
# different by a constant factor, the soil terms are in different units and
# every prediction would be wrong -- fix that before running 07.
xdata_c <- readRDS(file.path(DIR_DATA, "xdata_centered.rds"))
chk <- data.frame(northComponent = tab$north_component, eastComponent = tab$east_component,
                  cec30 = tab$cec30, clay30 = tab$clay30, sand30 = tab$sand30)
cmp <- data.frame(
  variable    = names(chk),
  grid_median = sapply(chk, median, na.rm = TRUE),
  fit_median  = sapply(names(chk), function(n) median(xdata_c[[n]], na.rm = TRUE)),
  grid_range  = sapply(chk, function(x) sprintf("%.3g to %.3g", min(x), max(x))),
  fit_range   = sapply(names(chk), function(n) sprintf("%.3g to %.3g", min(xdata_c[[n]], na.rm = TRUE),
                                                        max(xdata_c[[n]], na.rm = TRUE))),
  row.names = NULL)
cat("\n--- grid vs fitted data (medians and ranges should be the same order of magnitude) ---\n")
print(cmp)
w2w_range_check(chk, xdata_c, "static:")

## ---- quick look: is the domain really CONUS land (no Canada/Mexico/ocean)? ------
# Open out/static_quicklook.png. Left = land_frac on the whole grid, right = the
# cells that were kept. Check the northern and southern borders and the coasts.
r <- w2w_read(w2w_stem_static(), STATIC_BANDS)
kept <- r[["land_frac"]]; values(kept) <- NA; kept[tab$geeID] <- 1
png(file.path(DIR_OUT, paste0("static_quicklook", W2W_SUFFIX, ".png")), width = 1800, height = 800, res = 120)
par(mfrow = c(1, 2)); plot(r[["land_frac"]], main = "land_frac"); plot(kept, main = "cells kept", legend = FALSE)
dev.off()
message("wrote ", file.path(DIR_OUT, paste0("static_quicklook", W2W_SUFFIX, ".png")))

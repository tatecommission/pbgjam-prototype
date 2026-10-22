# =============================================================================
# 08_maps_per_species.R -- one GeoTIFF per species, per prediction table
# (i.e. per model/year/scenario). REVISED per Tate's instruction: NO
# combined/weighted total-community-BA layer anymore -- previously this
# script (then named 08_maps_and_weighted_ba.R) also summed all species
# into one total_community_BA.tif per tag; that step is removed entirely,
# not just made optional. Output is purely per-species maps now.
#
# Still newly written for this rewrite -- there is no earlier version of
# this script among what was uploaded at the start of this project, so
# review it more carefully than the rest of this workflow before trusting
# its output.
#
# Uses the static grid's lon/lat centroid columns (see
# gee/wall_to_wall_static.py) to place each geeID cell in space and
# rasterize with terra::rast(type = "xyz").
# =============================================================================

source("R/00_config.R")
library(data.table)

static_grid <- fread(file.path(DIR_DATA, "gee_exports", "gee_static_conus_1km.csv"))
stopifnot(all(c("geeID", "lon", "lat") %in% names(static_grid)))

species <- readRDS(file.path(DIR_DATA, "species.rds"))

pred_files <- list.files(file.path(DIR_OUT, "pred_tables"), pattern = "_mean\\.csv$", full.names = TRUE)
if (!length(pred_files)) stop("No prediction tables found in out/pred_tables/ -- run 07_predict_wall_to_wall_m1.R / _m2.R first.")

map_one <- function(pred_path) {
  tag  <- sub("_mean\\.csv$", "", basename(pred_path))
  pred <- fread(pred_path)

  coords <- merge(static_grid[, .(geeID, lon, lat)], pred, by = "geeID")
  if (!nrow(coords)) { message("[", tag, "] no matching geeIDs after merge -- skipping"); return(invisible(NULL)) }

  sp_cols <- intersect(paste0(species, "_mean"), names(coords))
  if (!length(sp_cols)) { message("[", tag, "] no species columns found -- skipping"); return(invisible(NULL)) }

  out_dir <- file.path(DIR_OUT, "pred_maps", tag)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ## ---- one GeoTIFF per species, this tag (model/year/scenario) only -----
  for (sp in sp_cols) {
    r <- terra::rast(coords[, c("lon", "lat", sp), with = FALSE], type = "xyz", crs = CRS_WGS)
    sp_name <- gsub("_mean$", "", sp)
    terra::writeRaster(r, file.path(out_dir, paste0(sp_name, ".tif")), overwrite = TRUE)
  }

  message(sprintf("[%s] wrote %d per-species rasters to %s", tag, length(sp_cols), out_dir))
}

invisible(lapply(pred_files, map_one))

message("\nAll available prediction tables rasterized -- one GeoTIFF per species per tag ",
        "(model/year/scenario). No combined/total-BA layer is produced, per instruction. ",
        "Check out/pred_maps/<tag>/<species>.tif.")

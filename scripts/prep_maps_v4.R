# prep_maps_v4.R -- run ONCE (RStudio: open and click Source).
# Splits the 5 big maps_v4 files (100 species each) into one single-band GeoTIFF per
# species, so the Shiny app can load one species in about a second instead of reading
# a 2.4 GB file. Output: maps_v4/by_species/<scenario>_<year>/<species>.tif (native 1 km).
# Safe to re-run: a finished set is skipped. Needs ~12 GB of free disk.
library(terra)

maps_dir <- path.expand("~/Documents/PBGJAM-data-explorer/maps_v4")
out_dir  <- file.path(maps_dir, "by_species")
sets <- c(baseline_2000 = "model2_baseline_2000_1km_mean.tif",
          ssp245_2060   = "model2_ssp245_2060_1km_mean.tif",
          ssp245_2090   = "model2_ssp245_2090_1km_mean.tif",
          ssp585_2060   = "model2_ssp585_2060_1km_mean.tif",
          ssp585_2090   = "model2_ssp585_2090_1km_mean.tif")
strip_rows <- 128L   # rows read at a time (about 0.5 GB of RAM)
group_size <- 50L    # species files kept open at once

for (tag in names(sets)) {
  dest <- file.path(out_dir, tag)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  src <- file.path(maps_dir, sets[[tag]])
  if (!file.exists(src)) stop("missing ", src)
  r   <- rast(src)
  nms <- names(r)
  if (nlyr(r) != 100L) stop(tag, ": expected 100 species layers, found ", nlyr(r))
  files <- file.path(dest, paste0(nms, ".tif"))
  if (all(file.exists(files))) { message(tag, ": already done, skipping"); next }
  message(tag, ": splitting ", length(nms), " species ...")
  t0 <- Sys.time()

  for (g in split(seq_along(nms), ceiling(seq_along(nms) / group_size))) {
    outs <- lapply(g, function(i) { o <- rast(r, nlyrs = 1); names(o) <- nms[i]; o })
    for (k in seq_along(g))
      writeStart(outs[[k]], files[g[k]], overwrite = TRUE, datatype = "FLT4S",
                 gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    readStart(r)
    for (start in seq(1L, nrow(r), by = strip_rows)) {
      n <- min(strip_rows, nrow(r) - start + 1L)
      v <- readValues(r, row = start, nrows = n, mat = TRUE)
      for (k in seq_along(g)) writeValues(outs[[k]], v[, g[k]], start, n)
    }
    readStop(r)
    for (k in seq_along(g)) writeStop(outs[[k]])
  }

  # check: first species, split file vs original layer
  a <- global(r[[1]], "mean", na.rm = TRUE)[1, 1]
  b <- global(rast(files[1]), "mean", na.rm = TRUE)[1, 1]
  if (!isTRUE(all.equal(a, b, tolerance = 1e-6))) stop(tag, ": check FAILED (", a, " vs ", b, ")")
  message(tag, ": done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min; check OK")
}
message("All sets written to ", out_dir)

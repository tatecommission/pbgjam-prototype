library(terra)

# ── Paths ─────────────────────────────────────────────────────────────────────
gjam_dir <- path.expand(
  "~/Documents/PBGJAM-data-explorer/data4Tate/gjamCreateMap/Trees/PredRasScale/rcp45/2040_2069"
)
out_dir <- "/Users/tatecommission/Documents/PBGJAM-data-explorer/pred_model_tool_prototype/prototype/dummy_rasters"

dir.create(file.path(out_dir, "timber"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "wildlife"), showWarnings = FALSE, recursive = TRUE)

# ── Build land-mask template at ~0.1° resolution ─────────────────────────────
# Use one GJAM raster for the CONUS land mask; resample to coarser target grid
template_path <- list.files(gjam_dir, pattern = "\\.tif$", full.names = TRUE)[1]
tmpl_fine <- terra::rast(template_path)   # 597 × 1388 at 0.0417°

# Target: ~0.1° resolution (≈2.4× coarser than GJAM, clearly visible at zoom 4)
tmpl <- terra::rast(
  nrows  = 249,
  ncols  = 579,
  extent = terra::ext(tmpl_fine),
  crs    = terra::crs(tmpl_fine)
)
# Resample to get land mask (NA = ocean/outside CONUS)
tmpl_mask <- terra::resample(tmpl_fine, tmpl, method = "near")
tmpl_mask <- terra::classify(tmpl_mask, cbind(NA, NA), others = 1)  # 1 = land

cat("Target grid:", nrow(tmpl_mask), "rows x", ncol(tmpl_mask), "cols\n")
cat("Land cells:", sum(!is.na(terra::values(tmpl_mask))), "\n\n")

# ── Species list from GJAM filenames ─────────────────────────────────────────
tif_files    <- list.files(gjam_dir, pattern = "\\.tif$", full.names = FALSE)
species_keys <- sub("^mean_", "", tif_files)
species_keys <- sub("_2040_2069_rcp45\\.tif$", "", species_keys)
cat("Generating rasters for", length(species_keys), "species × 2 layers\n\n")

# ── Spatial pattern generator ─────────────────────────────────────────────────
# Smooth random field: bilinear-resample a coarse random seed grid onto the
# target resolution, then mask to CONUS land. Produces spatially autocorrelated
# patterns (realistic gradients) that compress very well with DEFLATE.
make_smooth_raster <- function(seed_rows = 22, seed_cols = 50, template, mask) {
  seed <- terra::rast(
    nrows  = seed_rows,
    ncols  = seed_cols,
    extent = terra::ext(template),
    crs    = terra::crs(template)
  )
  terra::values(seed) <- runif(terra::ncell(seed))
  r <- terra::resample(seed, template, method = "bilinear")
  terra::mask(r, mask, maskvalues = NA)
}

# Write options: DEFLATE + floating-point predictor for best compression on [0,1] floats
write_opts <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")

# ── Generate all rasters ──────────────────────────────────────────────────────
set.seed(2025)

t_start <- proc.time()
for (i in seq_along(species_keys)) {
  sp <- species_keys[i]

  timber_r   <- make_smooth_raster(template = tmpl_mask, mask = tmpl_mask)
  wildlife_r <- make_smooth_raster(template = tmpl_mask, mask = tmpl_mask)

  terra::writeRaster(
    timber_r,
    file.path(out_dir, "timber",   paste0(sp, "_timber.tif")),
    overwrite = TRUE,
    gdal      = write_opts,
    datatype  = "FLT4S"
  )
  terra::writeRaster(
    wildlife_r,
    file.path(out_dir, "wildlife", paste0(sp, "_wildlife.tif")),
    overwrite = TRUE,
    gdal      = write_opts,
    datatype  = "FLT4S"
  )

  if (i %% 20 == 0 || i == length(species_keys)) {
    elapsed <- round((proc.time() - t_start)[3], 1)
    cat(sprintf("  %3d / %d  [%.1f s]\n", i, length(species_keys), elapsed))
  }
}

elapsed_total <- round((proc.time() - t_start)[3], 1)
t_files <- list.files(file.path(out_dir, "timber"),   pattern = "\\.tif$")
w_files <- list.files(file.path(out_dir, "wildlife"), pattern = "\\.tif$")
t_sizes <- file.size(file.path(out_dir, "timber",   t_files))
w_sizes <- file.size(file.path(out_dir, "wildlife", w_files))

cat("\n── Done ──────────────────────────────────────────\n")
cat(sprintf("Timber:   %d files, avg %.0f KB\n", length(t_files), mean(t_sizes) / 1024))
cat(sprintf("Wildlife: %d files, avg %.0f KB\n", length(w_files), mean(w_sizes) / 1024))
cat(sprintf("Total disk: %.1f MB\n", sum(t_sizes, w_sizes) / 1024^2))
cat(sprintf("Wall time:  %.1f s\n", elapsed_total))

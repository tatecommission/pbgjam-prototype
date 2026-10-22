# =============================================================================
# 07c_assemble_1km.R -- stitch the tiles of ONE 1 km prediction into GeoTIFFs and sanity-check it.
#   Rscript R/07c_assemble_1km.R <pred 1..26>
# Writes out/final_1km/<stem>_mean.tif and _se.tif (one layer per species, Float32, compressed),
# <stem>_total.tif (summed basal area), and compares the aggregated 1 km total with the 4 km run.
# Stops (loudly) if any tile is missing.
# =============================================================================
if (!exists("DIR_DATA")) source("R/00_config.R")
GRID_RES_M <- 1000L
source("R/w2w_helpers.R")
stopifnot(W2W_TAG == "1km")
source("R/w2w_1km_specs.R")
SPECS <- w2w_specs()
ARGS <- commandArgs(trailingOnly = TRUE)
if (length(ARGS) != 1 || !grepl("^[0-9]+$", ARGS)) stop("usage: Rscript R/07c_assemble_1km.R <pred 1..", nrow(SPECS), ">")
sp <- SPECS[as.integer(ARGS), ]
NT <- as.integer(Sys.getenv("W2W_NTILES", "30"))
outdir <- file.path(DIR_OUT, "final_1km"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

files <- vapply(seq_len(NT), function(t) w2w_tile_file(sp$stem, t, NT), "")
miss  <- !file.exists(files)
if (any(miss)) stop(sum(miss), " tile(s) missing for ", sp$stem, ": tiles ", paste(which(miss), collapse = ","),
                    " -- rerun those array tasks (finished tiles are skipped), then assemble again")
message("reading ", NT, " tiles for ", sp$stem)
ids <- vector("list", NT); mu <- vector("list", NT); se <- vector("list", NT)
for (t in seq_len(NT)) { x <- readRDS(files[t]); ids[[t]] <- x$geeID; mu[[t]] <- x$mean; se[[t]] <- x$se }
ids <- unlist(ids); M <- do.call(rbind, mu); rm(mu); S <- do.call(rbind, se); rm(se); invisible(gc())
stopifnot(!anyDuplicated(ids), nrow(M) == length(ids), nrow(S) == length(ids))
message(sprintf("%s cells x %d species", format(length(ids), big.mark = ","), ncol(M)))

tmpl <- w2w_read(w2w_stem_static(), STATIC_BANDS)[[1]]
write_stack <- function(X, out, nm) {
  parts <- character()
  for (b in split(seq_len(ncol(X)), ceiling(seq_len(ncol(X)) / 10))) {
    r <- rast(tmpl, nlyrs = length(b))
    v <- matrix(NA_real_, ncell(r), length(b)); v[ids, ] <- X[, b, drop = FALSE]
    values(r) <- v; names(r) <- nm[b]
    f <- file.path(outdir, sprintf("tmp_%s_%02d.tif", basename(out), length(parts) + 1))
    writeRaster(r, f, datatype = "FLT4S", overwrite = TRUE); parts <- c(parts, f)
  }
  writeRaster(rast(parts), out, datatype = "FLT4S", overwrite = TRUE,
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=YES"))
  unlink(parts)
}
sp_names <- sub("_mean$", "", colnames(M))
write_stack(M, file.path(outdir, paste0(sp$stem, "_mean.tif")), sp_names)
write_stack(S, file.path(outdir, paste0(sp$stem, "_se.tif")),   sp_names)
tot <- rowSums(M)
write_stack(matrix(tot, ncol = 1), file.path(outdir, paste0(sp$stem, "_total.tif")), "total_basal_area")
message(sprintf("total basal area per cell: mean %.2f | median %.2f | 99th pct %.2f | max %.2f",
                mean(tot), median(tot), quantile(tot, .99), max(tot)))
if (any(!is.finite(M)) || max(tot) > 500) warning("non-finite or implausible values in the assembled prediction")

## ---- compare with the 4 km run of the same prediction (aggregate 1 km -> 4 km) --------
f4 <- file.path(DIR_OUT, "pred_tables", paste0(sp$stem4, "_mean.csv"))
if (file.exists(f4)) {
  t4 <- data.table::fread(f4); id4 <- t4[[1]]; tot4 <- rowSums(as.matrix(t4[, -1]))
  rt <- rast(tmpl); values(rt) <- NA; rt[ids] <- tot
  r4 <- aggregate(rt, fact = 4, fun = "mean", na.rm = TRUE)        # 1214 x 788 grid = the 4 km grid
  a  <- values(r4)[id4, 1]; k <- is.finite(a) & is.finite(tot4)
  message(sprintf("vs 4 km run (%s cells): correlation of total basal area %.3f | mean 1 km/4 km ratio %.3f",
                  format(sum(k), big.mark = ","), cor(a[k], tot4[k]), mean(a[k]) / mean(tot4[k])))
} else message("no 4 km table found for comparison: ", f4)
message("assembled ", sp$stem, " -> ", outdir)

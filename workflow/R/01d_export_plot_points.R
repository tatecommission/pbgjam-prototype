# =============================================================================
# 01d_export_plot_points.R
#
# Stage 1 of the GEE round-trip, shared by BOTH new extraction scripts
# (02a_get_terraclimate_m2.R and 02b_get_remote_sensing_m1.R) -- there is
# only one set of FIA plot points, so it is exported once here rather than
# duplicated in each downstream script.
#
# REUSES the already-verified plot-to-cluster alignment from 01c: covar
# (covariates_4_5_2021.rds) is confirmed to align positionally with
# km$cluster. No new key-recovery work needed here -- identical to how the
# pre-rewrite 01d_get_nasa_climate.R did this.
# =============================================================================

source("R/00_config.R")

covar <- readRDS(file.path(DIR_FIA_SERVER, "covariates_4_5_2021.rds"))
km_path <- file.path(DIR_FIA_SERVER, "kmeans_cluster_4_5_2021.rds")
stopifnot("kmeans_cluster_4_5_2021.rds not found" = file.exists(km_path))
km <- readRDS(km_path)
stopifnot("covar row count does not match km$cluster length -- re-verify alignment (see 01c)" =
            nrow(covar) == length(km$cluster))

plot_points <- data.frame(rowid = seq_len(nrow(covar)), lon = covar$lon, lat = covar$lat)
pts_path <- file.path(DIR_DATA, "gee_exports", "plot_points_for_climate.csv")
write.csv(plot_points, pts_path, row.names = FALSE)

message("wrote ", pts_path, " (", nrow(plot_points), " plots)")
message("\nNEXT (manual, outside R), ONCE:")
message("  1. Upload ", pts_path, " as an Earth Engine table asset")
message("     (Code Editor Assets tab, or `earthengine upload table`).")
message("  2. Run BOTH extraction scripts against that SAME asset:")
message("       python gee/extract_terraclimate_at_plots.py   --asset <asset-id> --out terraclimate_at_plots")
message("       python gee/extract_remote_sensing_at_plots.py --asset <asset-id> --out remote_sensing_at_plots")
message("     Both are far lighter than the old 720-band NEX-GDDP export --")
message("     each is under 25 bands total -- so no special tileScale tuning")
message("     or multi-year test run should be needed, but consider a 2-year")
message("     test on the remote-sensing script first since it touches four")
message("     different collections at once.")
message("  3. Download both results to data/gee_exports/")
message("  4. Run 02a_get_terraclimate_m2.R and 02b_get_remote_sensing_m1.R")
message("     (each is also two-stage; they will pick up from here)")

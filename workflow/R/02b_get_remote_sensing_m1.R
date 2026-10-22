# =============================================================================
# 02b_get_remote_sensing_m1.R -- Stage 2: turn the downloaded remote-sensing
# CSV (gee/extract_remote_sensing_at_plots.py's output) into cluster-level
# precipGPM/LST/diffLST_airT for Model 1.
#
# REVISED: vegDRI REMOVED entirely (Tong's call -- see 00_config.R header).
# Three columns now, not four; no more sub-window inconsistency, since all
# three use the same FIT_YEARS window.
# =============================================================================

source("R/00_config.R")
library(data.table)

covar <- readRDS(file.path(DIR_FIA_SERVER, "covariates_4_5_2021.rds"))
km    <- readRDS(file.path(DIR_FIA_SERVER, "kmeans_cluster_4_5_2021.rds"))
load(XDATA_CLUSTER_RDATA)   # -> xdataCluster (for rownames / cluster count only)

rs_csv <- file.path(DIR_DATA, "gee_exports", "remote_sensing_at_plots.csv")

if (!file.exists(rs_csv)) {
  stop("data/gee_exports/remote_sensing_at_plots.csv not found. Run ",
      "01d_export_plot_points.R, then ",
      "gee/extract_remote_sensing_at_plots.py, then download the result here.")
}

rs <- fread(rs_csv)
stopifnot("remote_sensing_at_plots.csv row count doesn't match covar -- check the GEE export" =
            nrow(rs) == nrow(covar))

rs <- rs[order(rs$rowid), ]
stopifnot(identical(rs$rowid, seq_len(nrow(covar))))
stopifnot(all(c("precipGPM", "LST", "airT_TC", "diffLST_airT") %in% names(rs)))

## ---- verify the diff term server-side computation --------------------------
diff_check <- rs$LST - rs$airT_TC
max_mismatch <- max(abs(diff_check - rs$diffLST_airT), na.rm = TRUE)
if (max_mismatch > 0.01)
  warning(sprintf(
    "diffLST_airT does not equal LST - airT_TC (max mismatch %.3f) -- ",
    max_mismatch),
    "check for a scale-factor or unit error in the GEE extraction before proceeding.")

message(sprintf("plot-level precipGPM range: %.0f to %.0f mm/yr",
                min(rs$precipGPM, na.rm = TRUE), max(rs$precipGPM, na.rm = TRUE)))
message(sprintf("plot-level LST range: %.1f to %.1f C",
                min(rs$LST, na.rm = TRUE), max(rs$LST, na.rm = TRUE)))
message(sprintf("plot-level airT_TC range: %.1f to %.1f C",
                min(rs$airT_TC, na.rm = TRUE), max(rs$airT_TC, na.rm = TRUE)))
message(sprintf("plot-level diffLST_airT range: %.1f to %.1f C",
                min(rs$diffLST_airT, na.rm = TRUE), max(rs$diffLST_airT, na.rm = TRUE)))
message("\nSanity check LST is typically warmer than airT_TC for forested land ",
        "(mostly-positive diffLST_airT expected but not guaranteed) before proceeding.")

## ---- aggregate to cluster level -------------------------------------------
i <- km$cluster

precipGPM_cluster    <- signif(tapply(rs$precipGPM,    i, mean, na.rm = TRUE), 3)
LST_cluster          <- signif(tapply(rs$LST,          i, mean, na.rm = TRUE), 3)
diffLST_airT_cluster <- signif(tapply(rs$diffLST_airT, i, mean, na.rm = TRUE), 3)

stopifnot(length(precipGPM_cluster)    == nrow(xdataCluster),
         length(LST_cluster)          == nrow(xdataCluster),
         length(diffLST_airT_cluster) == nrow(xdataCluster))

model1_climate_cluster <- data.frame(
  precipGPM    = precipGPM_cluster[rownames(xdataCluster)],
  LST          = LST_cluster[rownames(xdataCluster)],
  diffLST_airT = diffLST_airT_cluster[rownames(xdataCluster)],
  row.names    = rownames(xdataCluster)
)

message(sprintf(
  "\ncluster-level NA counts -- precipGPM: %d, LST: %d, diffLST_airT: %d (of %d clusters)",
  sum(is.na(model1_climate_cluster$precipGPM)),
  sum(is.na(model1_climate_cluster$LST)),
  sum(is.na(model1_climate_cluster$diffLST_airT)),
  nrow(xdataCluster)))

saveRDS(model1_climate_cluster, file.path(DIR_DATA, "xdataCluster_model1_climate.rds"))
message("wrote data/xdataCluster_model1_climate.rds")
message("NEXT: set MODEL1_COLS_READY <- TRUE in R/00_config.R, then run 03_load_fit_data.R")

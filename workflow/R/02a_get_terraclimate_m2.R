# =============================================================================
# 02a_get_terraclimate_m2.R -- Stage 2: turn the downloaded TerraClimate CSV
# (gee/extract_terraclimate_at_plots.py's output) into cluster-level
# meanTemp/winTmin/annualPrec for Model 2.
#
# Same plot -> cluster aggregation convention as the pre-rewrite
# 01d_get_nasa_climate.R: plain arithmetic mean via tapply(), since none of
# these three variables are circular (unlike ASPECT).
# =============================================================================

source("R/00_config.R")
library(data.table)

covar <- readRDS(file.path(DIR_FIA_SERVER, "covariates_4_5_2021.rds"))
km    <- readRDS(file.path(DIR_FIA_SERVER, "kmeans_cluster_4_5_2021.rds"))
load(XDATA_CLUSTER_RDATA)   # -> xdataCluster (for rownames / cluster count only)

tc_csv <- file.path(DIR_DATA, "gee_exports", "terraclimate_at_plots.csv")

if (!file.exists(tc_csv)) {
  stop("data/gee_exports/terraclimate_at_plots.csv not found. Run ",
      "01d_export_plot_points.R, then ",
      "gee/extract_terraclimate_at_plots.py, then download the result here.")
}

tc <- fread(tc_csv)
stopifnot("terraclimate_at_plots.csv row count doesn't match covar -- check the GEE export" =
            nrow(tc) == nrow(covar))

tc <- tc[order(tc$rowid), ]
stopifnot(identical(tc$rowid, seq_len(nrow(covar))))
stopifnot(all(c("meanTemp", "winTmin", "annualPrec") %in% names(tc)))

message(sprintf("plot-level meanTemp range: %.1f to %.1f C",
                min(tc$meanTemp, na.rm = TRUE), max(tc$meanTemp, na.rm = TRUE)))
message(sprintf("plot-level winTmin range: %.1f to %.1f C",
                min(tc$winTmin, na.rm = TRUE), max(tc$winTmin, na.rm = TRUE)))
message(sprintf("plot-level annualPrec range: %.0f to %.0f mm",
                min(tc$annualPrec, na.rm = TRUE), max(tc$annualPrec, na.rm = TRUE)))
message("Sanity check these ranges by eye before proceeding -- e.g. meanTemp ",
        "well outside roughly 0-25C for a NC-scoped fit, or annualPrec well ",
        "outside roughly 900-2200mm, would suggest a scale-factor or unit error.")

## ---- aggregate to cluster level -------------------------------------------
i <- km$cluster   # positionally aligned with covar/tc, per 01c's verification

meanTemp_cluster   <- signif(tapply(tc$meanTemp,   i, mean, na.rm = TRUE), 3)
winTmin_cluster    <- signif(tapply(tc$winTmin,    i, mean, na.rm = TRUE), 3)
annualPrec_cluster <- signif(tapply(tc$annualPrec, i, mean, na.rm = TRUE), 3)

stopifnot(length(meanTemp_cluster)   == nrow(xdataCluster),
         length(winTmin_cluster)    == nrow(xdataCluster),
         length(annualPrec_cluster) == nrow(xdataCluster))

model2_climate_cluster <- data.frame(
  meanTemp   = meanTemp_cluster[rownames(xdataCluster)],
  winTmin    = winTmin_cluster[rownames(xdataCluster)],
  annualPrec = annualPrec_cluster[rownames(xdataCluster)],
  row.names  = rownames(xdataCluster)
)

message(sprintf("\ncluster-level NA counts -- meanTemp: %d, winTmin: %d, annualPrec: %d (of %d clusters)",
                sum(is.na(model2_climate_cluster$meanTemp)),
                sum(is.na(model2_climate_cluster$winTmin)),
                sum(is.na(model2_climate_cluster$annualPrec)),
                nrow(xdataCluster)))

saveRDS(model2_climate_cluster, file.path(DIR_DATA, "xdataCluster_model2_climate.rds"))
message("wrote data/xdataCluster_model2_climate.rds")
message("NEXT: set MODEL2_COLS_READY <- TRUE in R/00_config.R, then run 03_load_fit_data.R")

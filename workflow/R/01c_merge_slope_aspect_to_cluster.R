# =============================================================================
# 01c_merge_slope_aspect_to_cluster.R
#
# RESOLVED 2026-08-XX. This script previously stopped for manual confirmation
# of the plot -> cluster mapping. That investigation is DONE. Summary of what
# was found, for anyone re-reading this later:
#
#   - km$cluster (length 193,160) does NOT align positionally with xdata_topo
#     (196,765 rows, from 01b). That was the first hypothesis tried; it failed
#     an empirical lon/lat alignment check.
#   - The real input to kmeansEqualGroups() was covariates_4_5_2021.rds
#     ("covar"): exactly 193,160 rows, confirmed via lon/lat reconstruction
#     landing on the identity line against xdataCluster's own lon/lat.
#   - covar$plot is NOT a small FIA "PLOT" number -- despite looking like one
#     in head() (values like 3, 4, 5...), range() showed values up to 15
#     digits and length(unique(covar$plot)) == nrow(covar). It is PLT_CN
#     itself, stored as integer64, joinable directly against plotTopo$PLT_CN
#     from 01a. No state-abbreviation crosswalk needed.
#   - COND.csv has multiple rows (conditions) per PLT_CN; 01a's
#     mean(SLOPE, na.rm=TRUE) / circular.mean(ASPECT) already handles this
#     correctly by ignoring NA conditions per plot.
#
# Confirmed empirically with: row-count checks, random 20-row sample
# inspection, and a full lon/lat scatter against xdataCluster's own
# coordinates. See conversation history for the verification sequence if
# this ever needs to be re-derived (e.g. for a fresh xdataCluster vintage).
# =============================================================================

source("R/00_config.R")

## ---- 1. load everything relevant -------------------------------------------
plotTopo <- readRDS(file.path(DIR_DATA, "plotTopo_slope_aspect.rds"))      # from 01a
covar    <- readRDS(file.path(DIR_FIA_SERVER, "covariates_4_5_2021.rds"))  # kmeans input, confirmed
load(file.path(DIR_FIA_SERVER, "xdataCluster_4_5_2021.rdata"))             # -> xdataCluster
stopifnot(exists("xdataCluster"))

km_path <- file.path(DIR_FIA_SERVER, "kmeans_cluster_4_5_2021.rds")
stopifnot("kmeans_cluster_4_5_2021.rds not found" = file.exists(km_path))
km <- readRDS(km_path)

stopifnot("covar row count no longer matches km$cluster length -- re-verify alignment" =
            nrow(covar) == length(km$cluster))

## ---- 2. re-source circular.mean() exactly, matching 01a --------------------
circular.mean <- function(x){
  if (sum(!is.na(x)) == 0){
    circmean <- NA
  } else {
    x    <- x / 360 * 2 * pi
    sinr <- sum(sin(x), na.rm = TRUE)
    cosr <- sum(cos(x), na.rm = TRUE)
    circmean <- atan2(sinr, cosr) / (2 * pi) * 360
  }
  as.integer(circmean)
}

## ---- 3. join SLOPE/ASPECT onto covar by PLT_CN (covar$plot IS PLT_CN) -----
covar_topo <- merge(covar, plotTopo, by.x = "plot", by.y = "PLT_CN", all.x = TRUE)
stopifnot(nrow(covar_topo) == nrow(covar))  # no duplication/loss from the merge

# wrap ASPECT to 0-360 (atan2 returns -180..180) -- same fix applied in 01a
covar_topo$ASPECT[!is.na(covar_topo$ASPECT) & covar_topo$ASPECT < 0] <-
  covar_topo$ASPECT[!is.na(covar_topo$ASPECT) & covar_topo$ASPECT < 0] + 360

message(sprintf("SLOPE  matched: %d, NA: %d",
                sum(!is.na(covar_topo$SLOPE)), sum(is.na(covar_topo$SLOPE))))
message(sprintf("ASPECT matched: %d, NA: %d",
                sum(!is.na(covar_topo$ASPECT)), sum(is.na(covar_topo$ASPECT))))

## ---- 4. aggregate to cluster level, mirroring clusterFIA.R's own pattern --
# clusterFIA.R aggregates every other continuous variable with:
#   xk <- signif( tapply(xdata[,cc[k]], i, mean), 3 )
# SLOPE follows that exact pattern. ASPECT must NOT use plain mean() -- it's
# circular -- so it uses circular.mean() via split()/sapply() instead.
i <- km$cluster  # positionally aligned with covar_topo, confirmed above

SLOPE_cluster  <- signif(tapply(covar_topo$SLOPE, i, mean, na.rm = TRUE), 3)
ASPECT_cluster <- sapply(split(covar_topo$ASPECT, i), circular.mean)

stopifnot(length(SLOPE_cluster) == nrow(xdataCluster))
stopifnot(length(ASPECT_cluster) == nrow(xdataCluster))

## ---- 5. attach to xdataCluster, matched by cluster ID / rowname -----------
xdataCluster$SLOPE  <- SLOPE_cluster[rownames(xdataCluster)]
xdataCluster$ASPECT <- ASPECT_cluster[rownames(xdataCluster)]

message(sprintf("xdataCluster$SLOPE  NA: %d of %d",
                sum(is.na(xdataCluster$SLOPE)), nrow(xdataCluster)))
message(sprintf("xdataCluster$ASPECT NA: %d of %d",
                sum(is.na(xdataCluster$ASPECT)), nrow(xdataCluster)))

## ---- 6. write output --------------------------------------------------------
saveRDS(xdataCluster, file.path(DIR_DATA, "xdataCluster_with_slope_aspect.rds"))
message("wrote data/xdataCluster_with_slope_aspect.rds")
message("NEXT: set DEM_COLS_READY <- TRUE in R/00_config.R, then re-run ",
        "R/02_load_fit_data.R")
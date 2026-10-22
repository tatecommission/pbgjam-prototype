# =============================================================================
# 01b_link_slope_aspect_to_xdata_key.R
#
# RESOLVED 2026-08-XX. This script previously stopped for manual confirmation
# of the join key between COND.csv's PLT_CN and xdata.rdata's row key.
# Summary of what was found, for anyone re-reading this later:
#
#   - xdata$plot is an undocumented composite string: state abbreviation
#     (2 letters) + zero-padded PLT_CN (padded to match the longest PLT_CN
#     in the whole dataset). Traced to clusterFIA.R line 36:
#       plot <- paste(tree$STATE, plot, sep = '')
#     with the zero-padding done just above that (lines ~26-33).
#   - Confirmed by: stripping the 2-letter prefix and leading zeros,
#     converting to integer64, and matching against plotTopo$PLT_CN --
#     verified on a random 20-row sample (20/20 matched) and with a full
#     lon/lat scatter (xdata$lon vs. PLOT.csv's LAT for the recovered key)
#     landing on the identity line.
#
# See conversation history for the full verification sequence if this ever
# needs to be re-derived (e.g. a different xdata vintage with a different
# plot-key convention).
# =============================================================================

source("R/00_config.R")
library(data.table)
library(bit64)

## ---- 1. load the two sides of the bridge ----------------------------------
plotTopo <- readRDS(file.path(DIR_DATA, "plotTopo_slope_aspect.rds"))  # PLT_CN, SLOPE, ASPECT

plot_path <- file.path(DIR_FIA_SERVER, "PLOT.csv")
stopifnot("PLOT.csv not found -- check DIR_FIA_SERVER in 00_config.R" =
            file.exists(plot_path))
plot_tab <- fread(plot_path, select = c("CN", "STATECD", "PLOT", "LAT", "LON"))

load(file.path(DIR_FIA_SERVER, "xdata.rdata"))   # introduces `xdata`
stopifnot("expected object `xdata` after load()" = exists("xdata"))

message("xdata rows: ", nrow(xdata))

## ---- 2. recover PLT_CN from xdata$plot's composite key ---------------------
xdata$PLT_CN_recovered <- xdata$plot |>
  substr(3, nchar(xdata$plot)) |>   # drop 2-letter state prefix
  sub("^0+", "", x = _) |>          # drop leading zeros from padding
  as.integer64()

## ---- 3. merge slope/aspect onto xdata by the recovered key -----------------
xdata_with_topo <- merge(xdata, plotTopo,
                         by.x = "PLT_CN_recovered", by.y = "PLT_CN", all.x = TRUE)

## ---- 4. verify -------------------------------------------------------------
# Row count should be unchanged from the original xdata -- all.x = TRUE means
# unmatched plots keep their row with NA slope/aspect, not get dropped.
stopifnot(nrow(xdata_with_topo) == nrow(xdata))

message(sprintf("SLOPE  merged: %d matched, %d NA",
                sum(!is.na(xdata_with_topo$SLOPE)), sum(is.na(xdata_with_topo$SLOPE))))
message(sprintf("ASPECT merged: %d matched, %d NA",
                sum(!is.na(xdata_with_topo$ASPECT)), sum(is.na(xdata_with_topo$ASPECT))))

## ---- 5. write output ---------------------------------------------------------
saveRDS(xdata_with_topo, file.path(DIR_DATA, "xdata_with_slope_aspect.rds"))
message("wrote data/xdata_with_slope_aspect.rds")
message("NEXT: 01c_merge_slope_aspect_to_cluster.R maps these plots to ",
        "xdataCluster's spatial clusters.")
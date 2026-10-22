source("R/00_config.R")
library(data.table)

## ---- circular.mean(), copied verbatim from 01OrganizeData.R lines 7-18 ---
# Do not modify. If 01OrganizeData.R is ever updated, re-copy this block so
# the two stay identical.
circular.mean <- function(x){
  ### x in degree according to the FIA dataset
  ### return a degree
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

## ---- read only the columns needed -----------------------------------------
# COND.csv is 682 MB; select= avoids loading the other ~140 columns.
cond_path <- file.path(DIR_FIA_SERVER, "COND.csv")
stopifnot("COND.csv not found -- check DIR_FIA_SERVER in 00_config.R" =
            file.exists(cond_path))

cond <- fread(cond_path, select = c("PLT_CN", "SLOPE", "ASPECT"))

## ---- aggregate to plot level, exactly matching 01OrganizeData.R:75-76 ----
plotTopo <- cond[, .(
  SLOPE  = mean(SLOPE, na.rm = TRUE),
  ASPECT = circular.mean(ASPECT)
), by = PLT_CN]

## ---- sanity checks ---------------------------------------------------------
message(sprintf("plots with slope/aspect: %d", nrow(plotTopo)))
message(sprintf("SLOPE  range: %.1f to %.1f (percent)",
                min(plotTopo$SLOPE, na.rm = TRUE), max(plotTopo$SLOPE, na.rm = TRUE)))
message(sprintf("ASPECT range: %.0f to %.0f (degrees azimuth)",
                min(plotTopo$ASPECT, na.rm = TRUE), max(plotTopo$ASPECT, na.rm = TRUE)))
message(sprintf("NA SLOPE: %d   NA ASPECT: %d",
                sum(is.na(plotTopo$SLOPE)), sum(is.na(plotTopo$ASPECT))))

if (max(plotTopo$SLOPE, na.rm = TRUE) < 1)
  warning("Max SLOPE < 1 -- check this isn't a fraction (0-1) instead of ",
         "percent (0-100+). Original COND.csv field is percent.")
if (max(plotTopo$ASPECT, na.rm = TRUE) > 360 || min(plotTopo$ASPECT, na.rm = TRUE) < 0)
  warning("ASPECT outside 0-360 -- something upstream has changed units.")

saveRDS(plotTopo, file.path(DIR_DATA, "plotTopo_slope_aspect.rds"))
message("wrote data/plotTopo_slope_aspect.rds (plot-level, keyed by PLT_CN)")
message("NEXT: 01b_link_slope_aspect_to_xdata_key.R links this to the (state, plot) ",
       "key used in xdata.rdata.")

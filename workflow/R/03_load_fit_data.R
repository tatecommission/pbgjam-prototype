# =============================================================================
# 03_load_fit_data.R -- load the lab's existing fitting data and assemble
# ONE combined xdata containing the fixed base terms plus BOTH models'
# climate columns. Model 1's formulas reference the precipGPM/LST/
# diffLST_airT columns (vegDRI removed); Model 2's formulas reference meanTemp/
# winTmin/annualPrec -- both live in the same data frame so nothing has to
# be duplicated or kept in sync across two separate objects.
#
# This script:
#   1. loads xdataCluster / baCluster exactly as before
#   2. merges in DEM columns (SLOPE, ASPECT) once 01a/01b/01c have produced
#      them -- UNCHANGED from the pre-rewrite version of this script
#   3. decomposes SLOPE/ASPECT into north_component/east_component --
#      UNCHANGED
#   4. merges in Model 2's meanTemp/winTmin/annualPrec (02a) -- REPLACES
#      whatever xdataCluster originally shipped with for these three
#      (originally TerraClimate from years ago; briefly NEX-GDDP-CMIP6
#      during the abandoned 2026-08 migration; now TerraClimate again,
#      re-extracted fresh over FIT_YEARS 2001-2020 rather than reusing
#      either old version)
#   5. merges in Model 1's precipGPM/LST/diffLST_airT (02b) -- new columns,
#      xdataCluster never had these. vegDRI is NOT part of this rewrite --
#      removed per Tong's instruction (see 00_config.R header)
#   6. drops annualDef entirely, unconditionally -- unchanged decision from
#      before, still final
#   7. trims species with gjamTrimY(), matching the original script
# =============================================================================

source("R/00_config.R")
Rcpp::sourceCpp(file.path(GJAM_SRC, "src", "cppFns.cpp"))
for (f in list.files(file.path(GJAM_SRC, "R"), full.names = TRUE, pattern = "\\.[Rr]$"))
  source(f)

## ---- 1. load exactly as the original workflow does ------------------------
load(XDATA_CLUSTER_RDATA)   # -> xdataCluster
load(BA_CLUSTER_RDATA)      # -> baCluster

xdata <- xdataCluster
xdata[, "soil_type"] <- as.factor(xdata[, "soil_type"]) |> relevel("reference")

ydataAdult <- baCluster

## ---- 2. DEM columns (SLOPE, ASPECT) -- unchanged from before -------------
if (DEM_COLS_READY) {
  topo_path <- file.path(DIR_DATA, "xdataCluster_with_slope_aspect.rds")
  stopifnot("DEM_COLS_READY is TRUE but xdataCluster_with_slope_aspect.rds not found -- run 01a/01b/01c first" =
              file.exists(topo_path))
  xdata_topo <- readRDS(topo_path)
  stopifnot(all(DEM_COLS_PLACEHOLDER %in% names(xdata_topo)))
  stopifnot("row alignment mismatch between xdataCluster and xdataCluster_with_slope_aspect.rds -- do not proceed" =
              identical(rownames(xdata), rownames(xdata_topo)))
  xdata[, DEM_COLS_PLACEHOLDER] <- xdata_topo[, DEM_COLS_PLACEHOLDER]
} else {
  message("DEM columns not yet merged -- inserting NA placeholders: ",
          paste(DEM_COLS_PLACEHOLDER, collapse = ", "),
          ". Run 01a -> 01b -> 01c to build them.")
  for (col in DEM_COLS_PLACEHOLDER) xdata[[col]] <- NA_real_
}

## ---- 3. topographic decomposition -- unchanged from before ---------------
# north_component = sin(slope)*cos(aspect), east_component = sin(slope)*sin(aspect).
# xdataCluster's SLOPE is PERCENT GRADE (FIA field convention) -- atan(), not
# a direct degrees-to-radians conversion. See the wall-to-wall side
# (gee/wall_to_wall_static.py) for the DIFFERENT conversion used there
# (GEE's ee.Terrain.slope() already returns degrees) -- do not swap these.
slope_rad_fit  <- atan(xdata$SLOPE / 100)
aspect_rad_fit <- xdata$ASPECT * pi / 180

xdata$north_component <- sin(slope_rad_fit) * cos(aspect_rad_fit)
xdata$east_component  <- sin(slope_rad_fit) * sin(aspect_rad_fit)

## ---- 4. Model 2 climate: meanTemp, winTmin, annualPrec (TerraClimate) ----
MODEL2_CLIMATE_COLS <- c("meanTemp", "winTmin", "annualPrec")

if (MODEL2_COLS_READY) {
  m2_path <- file.path(DIR_DATA, "xdataCluster_model2_climate.rds")
  stopifnot("MODEL2_COLS_READY is TRUE but xdataCluster_model2_climate.rds not found -- run 02a first" =
              file.exists(m2_path))
  xdata_m2 <- readRDS(m2_path)
  stopifnot(all(MODEL2_CLIMATE_COLS %in% names(xdata_m2)))
  stopifnot("row alignment mismatch between xdataCluster and xdataCluster_model2_climate.rds -- do not proceed" =
              identical(rownames(xdata), rownames(xdata_m2)))
  xdata[, MODEL2_CLIMATE_COLS] <- xdata_m2[, MODEL2_CLIMATE_COLS]
  message("merged Model 2 climate columns (TerraClimate, ", min(FIT_YEARS), "-", max(FIT_YEARS), ")")
} else {
  message("Model 2 climate columns not yet ready -- NAing out meanTemp/winTmin/annualPrec. ",
          "Run 02a_get_terraclimate_m2.R to build them.")
  for (col in MODEL2_CLIMATE_COLS) xdata[[col]] <- NA_real_
}

## ---- 5. Model 1 climate: precipGPM, LST, diffLST_airT --------------------
# vegDRI REMOVED (Tong's call -- vegDRI and diffLST_airT "measure a
# similar thing"; see 00_config.R header). All three remaining columns now
# share the SAME FIT_YEARS window -- no more sub-window inconsistency.
MODEL1_CLIMATE_COLS <- c("precipGPM", "LST", "diffLST_airT")

if (MODEL1_COLS_READY) {
  m1_path <- file.path(DIR_DATA, "xdataCluster_model1_climate.rds")
  stopifnot("MODEL1_COLS_READY is TRUE but xdataCluster_model1_climate.rds not found -- run 02b first" =
              file.exists(m1_path))
  xdata_m1 <- readRDS(m1_path)
  stopifnot(all(MODEL1_CLIMATE_COLS %in% names(xdata_m1)))
  stopifnot("row alignment mismatch between xdataCluster and xdataCluster_model1_climate.rds -- do not proceed" =
              identical(rownames(xdata), rownames(xdata_m1)))
  xdata[, MODEL1_CLIMATE_COLS] <- xdata_m1[, MODEL1_CLIMATE_COLS]
  message("merged Model 1 remote-sensing columns (", min(FIT_YEARS), "-", max(FIT_YEARS), ")")
} else {
  message("Model 1 climate columns not yet ready -- NAing out precipGPM/LST/diffLST_airT. ",
          "Run 02b_get_remote_sensing_m1.R to build them.")
  for (col in MODEL1_CLIMATE_COLS) xdata[[col]] <- NA_real_
}

## ---- 6. drop annualDef entirely, unconditionally --------------------------
if ("annualDef" %in% names(xdata)) {
  xdata$annualDef <- NULL
  message("dropped annualDef (deficit) entirely -- unchanged decision from before")
}

## ---- 7. species trimming (matches original gjamTrimY call) ---------------
species    <- gjamTrimY(ydataAdult, maxCols = MAX_SPECIES, OTHER = FALSE)$y |> colnames()
ydataAdult <- ydataAdult[, species]

stopifnot(nrow(xdata) == nrow(ydataAdult))

saveRDS(xdata,      file.path(DIR_DATA, "xdata.rds"))
saveRDS(ydataAdult, file.path(DIR_DATA, "ydata_aligned.rds"))
saveRDS(species,    file.path(DIR_DATA, "species.rds"))

message(sprintf("\nloaded xdataCluster (n = %d, %d columns) and ydata (S = %d species)",
                nrow(xdata), ncol(xdata), length(species)))
message("xdata columns: ", paste(names(xdata), collapse = ", "))

if (!DEM_COLS_READY)
  message("\nNOTE: DEM columns are placeholders (all NA). Both models' formulas ",
          "reference north_component/east_component and will fail at fit time ",
          "until 01a/01b/01c have run.")
if (!MODEL2_COLS_READY)
  message("\nNOTE: Model 2 climate columns are placeholders (all NA) pending 02a.")
if (!MODEL1_COLS_READY)
  message("\nNOTE: Model 1 climate columns are placeholders (all NA) pending 02b.")

# =============================================================================
# 00_config.R  --  project settings. Source this first in every script.
#
# REVISION 2026-09 (v2), on top of the 2026-09 two-model rewrite. Three
# changes from that version:
#
#   1. VegDRI REMOVED from Model 1 entirely. Tong: VegDRI and diffLST_airT
#      "measure a similar thing" -- rather than carry two overlapping
#      drought/moisture-stress terms, Model 1 now relies on diffLST_airT
#      alone. This also removes the VegDRI-specific 2009-coverage-start
#      constraint that previously forced Model 1's drought term onto a
#      shorter sub-window than its other three terms -- Model 1 is now
#      internally consistent, ALL THREE remaining terms (precipGPM, LST,
#      diffLST_airT) over the same FIT_YEARS window.
#   2. Wall-to-wall prediction years now step every 5 years, 1950-2100
#      (Tate: "for now" -- every single year is a planned future upgrade,
#      not implemented here). PRED_YEARS_5YR is the single source of
#      truth; CMIP6_HISTORICAL_YEARS and PRED_YEARS_FUTURE are both
#      derived from it by splitting at CMIP6's historical/scenario
#      boundary (2014/2015), rather than maintained as two separately
#      hand-picked lists.
#   3. GCMS switched from a single model (ACCESS-CM2) to a 16-model
#      ensemble mean (16MME), per Wu, F., Jiao, D., Yang, X. et al. (2023)
#      "Evaluation of NEX-GDDP-CMIP6 in simulation performance and drought
#      capture utility over China -- based on DISO." Hydrology Research,
#      54(5), 703-721. https://doi.org/10.2166/nh.2023.140 (open access).
#      That paper's own conclusion: the 16-model ensemble mean outperforms
#      any single constituent model, including ACCESS-CM2, for both
#      precipitation and temperature simulation.
#      IMPORTANT CAVEAT, not silently glossed over: the paper's specific
#      performance rankings (which model is "best") were evaluated over
#      CHINA, not CONUS -- there is no guarantee the same ranking holds
#      over the US. The general finding that multi-model ensembles reduce
#      uncertainty relative to any single model is a broadly-supported,
#      domain-general result in the climate literature (the paper itself
#      cites Yang et al. 2020 for the same conclusion), so using this
#      specific 16-model list as a reasonable ensemble is a defensible
#      default -- but it is not the same as this paper having validated
#      16MME specifically for North America. Flag with Tong if a
#      CONUS-specific ensemble selection is wanted instead.
#      The paper describes the 16MME as built "by the weighted average
#      method" without stating the actual per-model weights anywhere in
#      the text available when this was written. GCM_ENSEMBLE_WEIGHTS
#      below defaults to EQUAL weights (simple arithmetic mean across all
#      16) -- NOT a reproduction of whatever weighting the paper actually
#      used, since that wasn't recoverable. Flagged, not silently assumed
#      to match.
#      ALSO NOT VERIFIED: that all 16 of these exact model-name strings
#      exist in NASA/GDDP-CMIP6's GEE "model" property. GEE silently
#      returns an empty ImageCollection for a filter that matches nothing
#      -- it does not error. Verify each name against the live catalog
#      (e.g. print the distinct "model" values in the collection) before
#      trusting a wall-to-wall export that includes all 16.
#
# REVISION 2026-09 (v3):
#   4. standAge REMOVED from both formulas entirely (Tong's instruction).
#      Nothing to grid or hold constant for it anymore.
#   5. Final wall-to-wall output is per-species-per-tag GeoTIFFs ONLY --
#      no combined/weighted total-community-BA layer. See
#      R/08_maps_per_species.R (renamed from 08_maps_and_weighted_ba.R).
#
# REVISION 2026-09 (v4):
#   6. moist REMOVED from both formulas entirely. The v3 fix (TreeMap
#      plot-ID linkage) turned out, on closer inspection, not to actually
#      provide a moist attribute at all -- TreeMap carries forest
#      structure/composition attributes (species, DBH, biomass, etc.) via
#      its FIA plot linkage, not a moisture variable, so the whole
#      TreeMap detour (R/07a_build_moist_wall_to_wall.R) is deleted, not
#      just deprioritized. Both formulas' base terms are now just
#      north_component + east_component + cec30 + [texture] -- no
#      standAge, no moist. Neither model has any remaining predictor
#      without a confirmed gridded wall-to-wall source.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(terra); library(sf); library(magrittr)
})

## ---- paths ----------------------------------------------------------------
PROJ     <- normalizePath("..", mustWork = FALSE)
DIR_DATA <- file.path(PROJ, "v4_gjam_workflow", "data")
DIR_OUT  <- file.path(PROJ, "out")
DIR_GEE  <- file.path(PROJ, "gee")

for (d in c(DIR_DATA, DIR_OUT, DIR_GEE,
            file.path(DIR_DATA, "fia"),
            file.path(DIR_DATA, "gee_exports"),
            file.path(DIR_OUT,  "fits"),
            file.path(DIR_OUT,  "pred_tables"),
            file.path(DIR_OUT,  "pred_maps"),
            file.path(DIR_OUT,  "formula_comparison")))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- gjam library location ------------------------------------------------
GJAM_SRC <- "~/Documents/PBGJAM-data-explorer/data4Tate/traitsMast/gjam_2.5.1"

## ---- existing fitting data (provided by lab, unchanged) -------------------
DIR_FIA_SERVER <- "/Volumes/Research/clark/clark.unix/fia/FIA2021"
XDATA_CLUSTER_RDATA <- file.path(DIR_FIA_SERVER, "xdataCluster_4_5_2021.rdata")
BA_CLUSTER_RDATA    <- file.path(DIR_FIA_SERVER, "baCluster_4_5_2021.rdata")

## ---- DEM columns -- UNCHANGED, still gated on 01a/01b/01c -----------------
DEM_COLS_PLACEHOLDER <- c("SLOPE", "ASPECT")
DEM_COLS_READY <- TRUE

## ---- fitting climate window ------------------------------------------------
# 2001-2020, unchanged. With VegDRI removed, Model 1 no longer has a
# sub-window inconsistency -- ALL of its terms (precipGPM, LST,
# diffLST_airT) now use this same window, same as Model 2.
FIT_YEARS <- 2001:2020

WINTER_MONTHS <- 1:3   # Model 2's winTmin only, unchanged

MODEL1_COLS_READY <- TRUE
MODEL2_COLS_READY <- TRUE

## ---- FIA scope --------------------------------------------------------------
STATES      <- c("NC")
MIN_DIA_IN  <- 5
MIN_PLOTS   <- 50
MAX_SPECIES <- 100

## ---- wall-to-wall grid -----------------------------------------------------
GRID_RES_M <- 1000

## ---- wall-to-wall prediction years: every 5 years, 1950-2100 (Tate) ------
# Single source of truth -- CMIP6_HISTORICAL_YEARS and PRED_YEARS_FUTURE
# below are BOTH derived from this by splitting at CMIP6's own
# historical/scenario boundary (2014/2015), not maintained separately.
# "For now" -- every single year is a planned future upgrade, not this one.
PRED_YEARS_5YR <- seq(1950, 2100, by = 5)

## ---- Model 1 historical wall-to-wall years --------------------------------
# Model 1 is observed-only (never touches CMIP6), so it's bound by actual
# data availability, not 1950. Core window: the 5-year grid points that
# fall within [min(FIT_YEARS), 2014] -- the only years where a genuine
# Model-1-vs-Model-2-CMIP6-historical-baseline comparison is possible (the
# CMIP6-historical constraint, not a VegDRI constraint anymore -- VegDRI's
# 2009 floor no longer applies to anything in this workflow).
MODEL1_WALLTOWALL_YEARS     <- PRED_YEARS_5YR[PRED_YEARS_5YR >= min(FIT_YEARS) & PRED_YEARS_5YR <= 2014]
# Extended range for a "just NASA" demonstration beyond 2014 -- upper bound
# set to match FIT_YEARS' own end (2020); verify each product's actual
# current end-of-record before trusting the top of this range.
MODEL1_WALLTOWALL_YEARS_EXT <- PRED_YEARS_5YR[PRED_YEARS_5YR >= min(FIT_YEARS) & PRED_YEARS_5YR <= 2020]

## ---- Model 2 / CMIP6 settings ----------------------------------------------
# 16-model ensemble (16MME), per Wu et al. (2023) -- see header. Exact
# Source ID strings from that paper's Table 1; NOT yet verified against
# the live NASA/GDDP-CMIP6 GEE catalog (see header caveat).
GCMS <- c("ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CESM2", "CESM2-WACCM",
         "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1", "IITM-ESM", "MIROC6",
         "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
         "NorESM2-MM", "TaiESM1")
GCM_ENSEMBLE_LABEL   <- "16MME"   # used in file naming, wall-to-wall exports
GCM_ENSEMBLE_WEIGHTS <- "equal"   # NOT verified to match the paper's own
                                  # (unstated) weighting -- see header

SCENARIOS <- c("ssp245", "ssp585")

# CMIP6 historical scenario covers 1950-2014; ssp scenarios take over at
# 2015 -- both derived from PRED_YEARS_5YR above, not hand-maintained.
CMIP6_HISTORICAL_YEARS <- PRED_YEARS_5YR[PRED_YEARS_5YR <= 2014]
PRED_YEARS_FUTURE      <- PRED_YEARS_5YR[PRED_YEARS_5YR >= 2015]

## ---- MCMC settings ----------------------------------------------------------
NG        <- 30000
BURNIN    <- 15000
REDUCT_R  <- 8
REDUCT_N  <- 25
NSIM_PRED <- 500

NG_EXPLORE     <- 1000
BURNIN_EXPLORE <- 500

## ---- CRS --------------------------------------------------------------------
CRS_WGS <- "EPSG:4326"

set.seed(20260901)
options(timeout = 3600)

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- v7: decade label = mean of yearly climate around it ----
W2W_WINDOW_BEFORE <- 5      # label-5 ... label+4, e.g. 2010 = 2005-2014
W2W_WINDOW_AFTER  <- 4

## ---- wall-to-wall pipeline: consolidated settings (later lines override earlier ones) ----
GRID_RES_M         <- 4000
FOREST_FRAC_MIN    <- 0
LAND_FRAC_MIN      <- 0.25
PRED_CHUNK         <- 20000L
NSIM_PRED          <- 100
GCM_ENSEMBLE_LABEL <- "13MME"
GCMS <- c("ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-CM2-SR5", "CMCC-ESM2",
          "CNRM-CM6-1", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0",
          "NorESM2-LM", "NorESM2-MM", "TaiESM1")
PRED_YEARS_5YR  <- seq(1950, 2100, by = 10)
MODEL1_WALLTOWALL_YEARS     <- PRED_YEARS_5YR[PRED_YEARS_5YR >= min(FIT_YEARS) & PRED_YEARS_5YR <= 2014]
MODEL1_WALLTOWALL_YEARS_EXT <- PRED_YEARS_5YR[PRED_YEARS_5YR >= min(FIT_YEARS) & PRED_YEARS_5YR <= 2020]
CMIP6_HISTORICAL_YEARS <- PRED_YEARS_5YR[PRED_YEARS_5YR <= 2014]
PRED_YEARS_FUTURE      <- PRED_YEARS_5YR[PRED_YEARS_5YR >= 2015]
W2W_WINDOW_BEFORE  <- 5
W2W_WINDOW_AFTER   <- 4
W2W_TEST           <- FALSE      # TRUE while testing; append "W2W_TEST <- FALSE" for the real run
W2W_OVERWRITE      <- FALSE

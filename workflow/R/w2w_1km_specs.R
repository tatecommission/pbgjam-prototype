# =============================================================================
# w2w_1km_specs.R -- the list of 26 predictions for the 1 km run, shared by
# 07_predict_1km_tiled.R and 07c_assemble_1km.R.  Same labels as the 4 km run:
#   1      = Model 1, historical observed (MODEL1_WALLTOWALL_YEARS, i.e. 2010)
#   2..8   = Model 2 baseline 1950..2010
#   9..17  = Model 2 ssp245 2020..2100
#   18..26 = Model 2 ssp585 2020..2100
# Needs 00_config.R (MODEL1_WALLTOWALL_YEARS, CMIP6_HISTORICAL_YEARS, SCENARIOS,
# PRED_YEARS_FUTURE) and w2w_helpers.R (W2W_TAG, W2W_SUFFIX) already loaded.
# =============================================================================
w2w_specs <- function() {
  s1 <- data.frame(model = 1L, scen = "historical", yr = MODEL1_WALLTOWALL_YEARS,
                   stem = sprintf("model1_historical_%d", MODEL1_WALLTOWALL_YEARS),
                   tag = sprintf("m1_historical_%d", MODEL1_WALLTOWALL_YEARS),
                   stringsAsFactors = FALSE)
  s2 <- data.frame(model = 2L, scen = "historical", yr = CMIP6_HISTORICAL_YEARS,
                   stem = sprintf("model2_baseline_%d", CMIP6_HISTORICAL_YEARS),
                   tag = sprintf("baseline_%d", CMIP6_HISTORICAL_YEARS),
                   stringsAsFactors = FALSE)
  s3 <- do.call(rbind, lapply(SCENARIOS, function(s)
    data.frame(model = 2L, scen = s, yr = PRED_YEARS_FUTURE,
               stem = sprintf("model2_%s_%d", s, PRED_YEARS_FUTURE),
               tag = sprintf("%s_%d", s, PRED_YEARS_FUTURE), stringsAsFactors = FALSE)))
  sp <- rbind(s1, s2, s3)
  sp$stem4 <- sp$stem                                       # name of the matching 4 km table
  sp$stem  <- paste0(sp$stem, "_", W2W_TAG, W2W_SUFFIX)     # e.g. model2_baseline_2010_1km
  sp$idx   <- seq_len(nrow(sp))
  sp
}
w2w_tile_file <- function(stem, tile, nt)
  file.path(DIR_OUT, "pred_tiles_1km", sprintf("%s_t%03d_of%03d.rds", stem, tile, nt))

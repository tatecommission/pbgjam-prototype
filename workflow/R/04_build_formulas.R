# =============================================================================
# 04_build_formulas.R -- build the fixed formula shortlist for EACH model.
# No VIF/correlation diagnostics step.
#
# NAMING: gjam rejects '_' and '-' in formula variable names ("reserved
# characters '_' or '-' in formula variables, rename" -- discovered when
# 05_compare_formulas.R failed on every single fit). All predictor names
# below are therefore underscore-free camelCase: north_component ->
# northComponent, east_component -> eastComponent, diffLST_airT ->
# diffLSTairT, and the centered-quadratic suffix is "Ctr" instead of "_c"
# (e.g. precipGPM_c -> precipGPMCtr). xdata's column names are renamed to
# match immediately after loading, so everything downstream -- formulas,
# saved xdata_centered.rds, predictor_centers.rds -- is consistently
# underscore-free.
#
# MODEL 1 -- fixed base (vegDRI removed, standAge removed, moist removed --
# see 00_config.R header for all three):
#   ~ northComponent + eastComponent + cec30 +
#     [TEXTURE] + diffLSTairT
#   diffLSTairT is LINEAR ONLY and present in ALL formulas -- not a
#   toggled axis. Toggled: texture (2) x precipitation=GPM (2) x
#   temperature=LST (2) = 8.
#
# MODEL 2 -- fixed base (standAge removed, moist removed):
#   ~ northComponent + eastComponent + cec30 +
#     [TEXTURE]
#   Toggled: texture (2) x temperature (4) x precipitation=annualPrec (2) = 16
# =============================================================================

source("R/00_config.R")

xdata <- readRDS(file.path(DIR_DATA, "xdata.rds"))

## ---- rename underscore-containing columns (gjam-safe names) ---------------
rename_map <- c(
  north_component = "northComponent",
  east_component  = "eastComponent",
  diffLST_airT    = "diffLSTairT"
)
matched <- intersect(names(rename_map), names(xdata))
if (length(matched)) {
  names(xdata)[match(matched, names(xdata))] <- rename_map[matched]
}
message("renamed for gjam compatibility: ",
        paste(matched, "->", rename_map[matched], collapse = ", "))

REQUIRED_BASE <- c("northComponent", "eastComponent", "cec30",
                   "sand30", "clay30")
missing_base <- setdiff(REQUIRED_BASE, names(xdata))
if (length(missing_base))
  stop("xdata is missing expected base columns: ", paste(missing_base, collapse = ", "),
       " -- these should already exist in xdataCluster; check 03_load_fit_data.R output.")

if (!DEM_COLS_READY)
  message("WARNING: DEM_COLS_READY is FALSE -- both formula sets reference ",
          "northComponent/eastComponent, currently NA.")
if (!MODEL1_COLS_READY)
  message("WARNING: MODEL1_COLS_READY is FALSE -- Model 1 formulas reference ",
          "columns currently NA.")
if (!MODEL2_COLS_READY)
  message("WARNING: MODEL2_COLS_READY is FALSE -- Model 2 formulas reference ",
          "columns currently NA.")

## ---- centering (quadratic terms only) -------------------------------------
# NOTE: suffix is "Ctr", not "_c" -- gjam rejects underscores in formula vars.
center_vars <- function(df, vars) {
  centers <- setNames(sapply(vars, function(v) mean(df[[v]], na.rm = TRUE)), vars)
  for (v in vars) df[[paste0(v, "Ctr")]] <- df[[v]] - centers[[v]]
  list(df = df, centers = centers)
}

# vegDRI REMOVED from this list -- no longer exists as a predictor.
QUAD_VARS <- c("meanTemp", "winTmin", "annualPrec", "precipGPM", "LST")
cx      <- center_vars(xdata, QUAD_VARS)
xdata_c <- cx$df
CENTERS <- cx$centers
saveRDS(CENTERS, file.path(DIR_DATA, "predictor_centers.rds"))
saveRDS(xdata_c, file.path(DIR_DATA, "xdata_centered.rds"))

BASE_TERMS <- c("northComponent", "eastComponent", "cec30")
texture_choices <- list(sand30 = "sand30", clay30 = "clay30")

## =============================================================================
## MODEL 1 -- 8 formulas (texture x precip x temp; drought axis removed)
## =============================================================================
m1_precip_choices <- list(
  precipGPM     = "precipGPM",
  precipGPMquad = c("precipGPMCtr", "I(precipGPMCtr^2)")
)
m1_temp_choices <- list(
  LST     = "LST",
  LSTquad = c("LSTCtr", "I(LSTCtr^2)")
)

FORMULAS_M1 <- list()
i <- 0
for (tx_n in names(texture_choices)) {
  for (pc_n in names(m1_precip_choices)) {
    for (tp_n in names(m1_temp_choices)) {
      i <- i + 1
      terms <- c(BASE_TERMS, texture_choices[[tx_n]], "diffLSTairT",
                 m1_precip_choices[[pc_n]], m1_temp_choices[[tp_n]])
      form_name <- sprintf("M1_%02d_%s_%s_%s", i, tx_n, pc_n, tp_n)
      FORMULAS_M1[[form_name]] <- as.formula(paste("~", paste(terms, collapse = " + ")))
    }
  }
}
stopifnot(length(FORMULAS_M1) == 8)

## =============================================================================
## MODEL 2 -- 16 formulas, UNCHANGED (aside from Ctr suffix rename)
## =============================================================================
m2_temp_choices <- list(
  meanTemp     = "meanTemp",
  meanTempquad = c("meanTempCtr", "I(meanTempCtr^2)"),
  winTmin      = "winTmin",
  winTminquad  = c("winTminCtr", "I(winTminCtr^2)")
)
m2_precip_choices <- list(
  annualPrec     = "annualPrec",
  annualPrecquad = c("annualPrecCtr", "I(annualPrecCtr^2)")
)

FORMULAS_M2 <- list()
i <- 0
for (tx_n in names(texture_choices)) {
  for (tp_n in names(m2_temp_choices)) {
    for (pc_n in names(m2_precip_choices)) {
      i <- i + 1
      terms <- c(BASE_TERMS, texture_choices[[tx_n]], m2_temp_choices[[tp_n]],
                 m2_precip_choices[[pc_n]])
      form_name <- sprintf("M2_%02d_%s_%s_%s", i, tx_n, tp_n, pc_n)
      FORMULAS_M2[[form_name]] <- as.formula(paste("~", paste(terms, collapse = " + ")))
    }
  }
}
stopifnot(length(FORMULAS_M2) == 16)

message(sprintf("built %d Model 1 formulas (vegDRI removed) and %d Model 2 formulas",
                length(FORMULAS_M1), length(FORMULAS_M2)))

saveRDS(FORMULAS_M1, file.path(DIR_DATA, "candidate_formulas_m1.rds"))
saveRDS(FORMULAS_M2, file.path(DIR_DATA, "candidate_formulas_m2.rds"))

cat("\n--- Model 1 (remote sensing), example formula ---\n")
print(FORMULAS_M1[[1]])
cat("\n--- Model 2 (climate-only), example formula ---\n")
print(FORMULAS_M2[[1]])
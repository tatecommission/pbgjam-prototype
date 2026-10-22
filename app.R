# =============================================================================
# PBGJAM v2 — Decision Making Tool (Shiny)
# =============================================================================

library(shiny)
library(leaflet)
library(sf)
library(tigris)
library(terra)
library(raster)
library(jsonlite)
library(tidyverse)

load("BA.rdata")
load("xdata.rdata")

# Wildlife / pest / riparian table (v2-wildlife-pests-workflow output)
WILDLIFE_FILE <- "BA_wildlife_traits_v2.csv"
# Army Corps wetland region used for riparian status (match 00_config.R)
WETLAND_COL   <- "wl_ncne"

# One entry per wildlife slider; `col` is an equal-step 0-1 score in the CSV.
WILD_COMPONENTS <- list(
  list(id="food",    col="wild_food_seeds_fruit", label="Fruit, nut & seed food",
       note="USDA PLANTS", on=TRUE,  w=0.165),
  list(id="pollin",  col="wild_pollinators",      label="Pollinators",
       note="PA DCNR Native Tree Selector", on=TRUE,  w=0.165),
  list(id="deer",    col="wild_deer_forage",      label="Deer forage / resistance",
       note="0 = resistant, 1 = forage.", on=FALSE, w=0.11)
)

pest_source_label <- function(x) {
  ifelse(is.na(x), NA_character_,
  ifelse(x == "PA Native Tree Selector", "PA DCNR Native Tree Selector",
         "Pennsylvania pest and pathogen ratings"))
}

wildlife_data <- read_csv(WILDLIFE_FILE, show_col_types = FALSE) |>
  dplyr::rename(ba_code = code, scientific_name = binomial, do_not_plant = client_do_not_plant) |>
  dplyr::mutate(
    dplyr::across(dplyr::any_of(c(
      vapply(WILD_COMPONENTS, function(x) x$col, character(1)),
      "pest_resist_score", "do_not_plant", "pa_common_riparian", WETLAND_COL
    )), as.numeric),
    riparian_recommended =
      (!is.na(.data[[WETLAND_COL]]) & .data[[WETLAND_COL]] >= 0.75) |
      (!is.na(pa_common_riparian) & pa_common_riparian == 1),
    has_pest_flag =
      (!is.na(do_not_plant) & do_not_plant == 1) |
      (!is.na(pest_resist_score) & pest_resist_score == 0),
    pest_text_source = pest_source_label(pest_text_source)
  ) |>
  dplyr::select(ba_code, scientific_name,
                dplyr::all_of(vapply(WILD_COMPONENTS, function(x) x$col, character(1))),
                pest_text_display, pest_text_source,
                riparian_recommended, has_pest_flag)


coords_lon <- xdata$lon
coords_lat <- xdata$lat
BA_mat <- BA

local({
  df <- as.data.frame(BA_mat)
  ok <- vapply(df, function(col) is.numeric(col) || is.integer(col), logical(1))
  if (!all(ok)) {
    message(sprintf(
      "BA: dropping %d non-numeric column(s) before matrix conversion: %s",
      sum(!ok), paste(names(df)[!ok], collapse=", ")
    ))
  }
  m <- as.matrix(df[, ok, drop = FALSE])
  mode(m) <- "double"
  m[is.na(m)] <- 0
  BA_mat <<- m
})

timber <- read_csv(
  "timber_price_by_gjam_species_wide.csv",
  show_col_types = FALSE
)

# How each price cell was resolved, same shape/order as `timber`. Used to label cells that
# come from a regional aggregate rather than the species itself.
timber_src <- read_csv(
  "timber_price_source_by_gjam_species_wide.csv",
  show_col_types = FALSE
)

fmt_species_fallback <- function(nm) {
  parts <- strsplit(nm, "(?<=[a-z])(?=[A-Z])", perl = TRUE)[[1]]
  if (length(parts) < 2)
    return(paste0(toupper(substr(nm,1,1)), tolower(substr(nm,2,nchar(nm)))))
  genus   <- paste0(toupper(substr(parts[1],1,1)), tolower(substr(parts[1],2,nchar(parts[1]))))
  epithet <- tolower(paste(parts[-1], collapse=""))
  paste(genus, epithet)
}

species_names <- colnames(BA_mat)
wl_lookup <- setNames(as.character(wildlife_data$scientific_name), as.character(wildlife_data$ba_code))

species_display <- setNames(
  vapply(species_names, function(nm) {
    sn <- as.character(wl_lookup[nm])
    if (length(sn) == 1L && !is.na(sn) && nzchar(trimws(sn))) {
      raw <- trimws(sn)
      parts <- strsplit(raw, "\\s+")[[1]]
      parts <- parts[nzchar(parts)]
      if (length(parts) >= 2) {
        paste0(toupper(substr(parts[1],1,1)), tolower(substr(parts[1],2,nchar(parts[1]))),
               " ", paste(tolower(parts[-1]), collapse=" "))
      } else if (length(parts) == 1) {
        paste0(toupper(substr(parts[1],1,1)), tolower(substr(parts[1],2,nchar(parts[1]))))
      } else {
        fmt_species_fallback(nm)
      }
    } else {
      fmt_species_fallback(nm)
    }
  }, FUN.VALUE = character(1), USE.NAMES = FALSE),
  species_names
)


# ── maps_v4: 1 km basal-area maps (Model 2) — Visualizations tab only ───────────
# One single-band GeoTIFF per species per year/scenario, made once by prep_maps_v4.R
# from the 5 big files in maps_v4/. These maps drive the map view AND the species ranking.
.ba4_dir <- path.expand("~/Documents/PBGJAM-data-explorer/maps_v4/by_species")
ba4_path <- function(key, yr, scen) {
  yr <- as.character(yr)
  tag <- if (identical(yr, "2000")) "baseline_2000" else paste0(scen, "_", yr)
  file.path(.ba4_dir, tag, paste0(key, ".tif"))
}
ba4_keys <- local({
  d <- file.path(.ba4_dir, "baseline_2000")
  if (dir.exists(d)) sub("\\.tif$", "", list.files(d, pattern = "\\.tif$")) else character(0)
})
if (length(ba4_keys) == 0)
  message("maps_v4 per-species maps not found at: ", .ba4_dir, "  (run prep_maps_v4.R once)")
# Species without their own map (the model covers 100 of the ~490 species in the data) use the
# AVERAGE of the modeled species in the same genus (the genus = the leading lower-case letters of
# the species code, e.g. quercusPalustri -> quercus). A genus with no modeled species stays NA.
ba4_genus <- function(k) sub("^([a-z]+).*$", "\\1", k)
ba4_has_map <- function(k) k %in% ba4_keys
ba4_genus_members <- function(k) ba4_keys[ba4_genus(ba4_keys) == ba4_genus(k)]
# Climate window behind each map label L = years (L - BA_W_BEFORE) .. (L + BA_W_AFTER).
# UNVERIFIED: must equal W2W_WINDOW_BEFORE / W2W_WINDOW_AFTER in the 00_config.R that the 1 km
# run used on the DCC (check: grep -n W2W_WINDOW /work/tvc11/PBGJAM/v4_gjam_workflow/R/00_config.R).
# 5 / 4 = 10-year windows (2060 -> 2055-2064). 10 / 9 would be 20-year windows (2060 -> 2050-2069).
# If you change it here, change WIN_BEFORE / WIN_AFTER in prep_predictors.R too and re-run it.
BA_W_BEFORE <- 5
BA_W_AFTER  <- 4
ba_period_text <- function(yr) {
  L <- suppressWarnings(as.integer(as.character(yr)))
  if (length(L) != 1L || is.na(L)) return("")
  sprintf("%d (%d–%d)", L, L - BA_W_BEFORE, L + BA_W_AFTER)
}
# Fixed scale per species, identical for every period and scenario (made by prep_ba_ranges.R):
# cap = the highest 99.5th percentile of that species' non-zero cells over the 5 maps.
ba4_caps <- local({
  f <- file.path(.ba4_dir, "ba_ranges.csv")
  if (file.exists(f)) {
    d <- utils::read.csv(f, stringsAsFactors = FALSE); stats::setNames(d$cap, d$species)
  } else {
    message("ba_ranges.csv not found - run prep_ba_ranges.R. Until then maps use per-map colour ranges and the species ranking is empty.")
    numeric(0)
  }
})
ba4_cap <- function(key) {
  v <- unname(ba4_caps[key])
  if (length(v) == 1L && is.finite(v) && v > 0) v else NA_real_
}
ba4_clamped_pal <- function(pal0, rng) function(v) pal0(pmin(pmax(v, rng[1]), rng[2]))

# Period / scenario selector. Every copy on the page shares the inputs vis_year / vis_scen,
# so Location, Decision-Making and Visualizations always agree.
ba4_selector <- function(dark = FALSE) {
  div(class = if (dark) "ba4-ctl ba4-dark" else "ba4-ctl",
      div(class = "ba4-group", div(class = "ba4-label", "Period"),
          div(class = "ba4-seg-row", `data-inputid` = "vis_year",
              div(class = "ba4-seg", `data-val` = "2000", "2000"),
              div(class = "ba4-seg active", `data-val` = "2060", "2060"),
              div(class = "ba4-seg", `data-val` = "2090", "2090"))),
      div(class = "ba4-group", div(class = "ba4-label", "Scenario"),
          div(class = "ba4-seg-row ba4-scen", `data-inputid` = "vis_scen",
              div(class = "ba4-seg active", `data-val` = "ssp245", "SSP245"),
              div(class = "ba4-seg", `data-val` = "ssp585", "SSP585"))))
}

# Species scores: mean (polygon) or cell (point) value of each species' map at the site.
.gjam4_cache <- new.env(hash = TRUE, parent = emptyenv())
get_gjam_vals_v4 <- function(site, keys, yr, scen) {
  poly_sig <- if (identical(site$type, "polygon")) paste(round(c(site$poly_lng, site$poly_lat), 5), collapse = ",") else ""
  ck <- paste(site$type, round(site$lng, 5), round(site$lat, 5), poly_sig,
              paste0(as.character(yr), "_", scen), length(keys), keys[1], sep = "|")
  if (exists(ck, envir = .gjam4_cache, inherits = FALSE))
    return(get(ck, envir = .gjam4_cache))
  pt   <- terra::vect(matrix(c(site$lng, site$lat), ncol = 2), crs = "EPSG:4326")
  poly <- NULL
  if (identical(site$type, "polygon") && length(site$poly_lng) >= 3)
    poly <- tryCatch(site_polygon_vect(site), error = function(e) NULL)
  extract_one <- function(f) tryCatch({
    r  <- terra::rast(f); rc <- terra::crs(r)
    use_pt   <- terra::project(pt, rc)
    use_poly <- if (!is.null(poly)) terra::project(poly, rc) else NULL
    v <- NA_real_
    if (!is.null(use_poly))
      v <- suppressWarnings(as.numeric(terra::extract(r, use_poly, fun = mean, na.rm = TRUE)[1, -1]))
    if (!is.finite(v)) v <- as.numeric(terra::extract(r, use_pt)[1, -1])
    v
  }, error = function(e) NA_real_)
  own <- function(k) { f <- ba4_path(k, yr, scen); if (file.exists(f)) extract_one(f) else NA_real_ }
  gmemo <- list()                                  # genus -> average of its modeled species here
  genus_avg <- function(k) {
    g <- ba4_genus(k)
    if (is.null(gmemo[[g]])) {
      v <- vapply(ba4_genus_members(k), own, numeric(1)); v <- v[is.finite(v)]
      gmemo[[g]] <<- if (length(v)) mean(v) else NA_real_
    }
    gmemo[[g]]
  }
  fb   <- character(0)
  vals <- vapply(keys, function(k) {
    if (ba4_has_map(k)) return(own(k))
    v <- genus_avg(k)
    if (is.finite(v)) fb <<- c(fb, k)
    v
  }, numeric(1))
  res <- list(values = vals, fallback_keys = fb)
  if (length(ls(.gjam4_cache)) > 300L) rm(list = ls(.gjam4_cache), envir = .gjam4_cache)
  assign(ck, res, envir = .gjam4_cache)
  res
}

ba4_title <- function(yr, scen) {
  if (is.null(yr) || identical(as.character(yr), "2000"))
    return(paste0("Basal area, ", ba_period_text(2000), " (m²/ha)"))
  sc <- if (identical(scen, "ssp585")) "SSP585" else "SSP245"
  paste0("Predicted BA, ", ba_period_text(yr), ", ", sc, " (m²/ha)")
}

# Draw only what is on screen: crop the single-species map to the current view (plus a
# margin), keep native 1 km unless that would be more than BA4_MAX_CELLS cells (whole-US
# view), then project to Web Mercator once with terra (fast) for leaflet.
BA4_MAX_CELLS <- 1.5e6
.ba4_rng <- new.env(hash = TRUE, parent = emptyenv())
# Bounded cache for the large in-memory rasters (whole-US copies, predictor maps): keeps the
# 24 most recently used and drops the oldest, so memory stays flat however long you browse.
.ba4_lru <- new.env(parent = emptyenv()); .ba4_lru$keys <- character(0)
BA4_CACHE_MAX <- 24L
ba4_cache_put <- function(key, val) {
  assign(key, val, envir = .ba4_rng)
  ord <- c(setdiff(.ba4_lru$keys, key), key)
  while (length(ord) > BA4_CACHE_MAX) {
    if (exists(ord[1], envir = .ba4_rng, inherits = FALSE)) rm(list = ord[1], envir = .ba4_rng)
    ord <- ord[-1]
  }
  .ba4_lru$keys <- ord
  invisible(val)
}
ba4_cache_get <- function(key) {
  if (!exists(key, envir = .ba4_rng, inherits = FALSE)) return(NULL)
  val <- get(key, envir = .ba4_rng); ba4_cache_put(key, val); val
}
ba4_view_raster <- function(r, b) {
  e <- terra::ext(r)
  if (!is.null(b) && all(c("north", "south", "east", "west") %in% names(b))) {
    dx <- (b$east - b$west) * 0.25; dy <- (b$north - b$south) * 0.25
    g  <- expand.grid(x = seq(max(b$west - dx, -179), min(b$east + dx, 179), length.out = 15),
                      y = seq(max(b$south - dy, -85), min(b$north + dy, 85), length.out = 15))
    pts <- terra::project(terra::vect(as.matrix(g), crs = "EPSG:4326"), terra::crs(r))
    e   <- tryCatch(terra::intersect(terra::ext(pts), terra::ext(r)), error = function(err) NULL)
    if (is.null(e) || terra::xmin(e) >= terra::xmax(e) || terra::ymin(e) >= terra::ymax(e)) return(NULL)
  }
  full_area <- (terra::xmax(r) - terra::xmin(r)) * (terra::ymax(r) - terra::ymin(r))
  if ((terra::xmax(e) - terra::xmin(e)) * (terra::ymax(e) - terra::ymin(e)) > 0.6 * full_area) {
    ck <- paste0("coarse_", terra::sources(r)[1])
    hit <- ba4_cache_get(ck)
    if (!is.null(hit)) return(hit)
    cr  <- terra::aggregate(r, fact = 4L, fun = "mean", na.rm = TRUE)
    out <- raster::raster(terra::project(cr, "EPSG:3857", method = "bilinear"))
    ba4_cache_put(ck, out)
    return(out)
  }
  rc  <- terra::crop(r, e)
  fct <- max(1L, ceiling(sqrt(terra::ncell(rc) / BA4_MAX_CELLS)))
  assign("last_fct", fct, envir = .ba4_rng)
  if (fct > 1L) rc <- terra::aggregate(rc, fact = fct, fun = "mean", na.rm = TRUE)
  raster::raster(terra::project(rc, "EPSG:3857", method = "bilinear"))
}

# ── Predictor maps (4 km, made once by prep_predictors.R) ──────────────────────────
# Climate: same window rule as the model; winTmin = coldest value in the window, annualPrec =
# mean over the window; one shared colour range per variable. Soil / topography are static.
.pred_dir <- path.expand("~/Documents/PBGJAM-data-explorer/maps_v4/predictors")
PRED_INFO <- list(
  pred_wintmin = list(var = "winTmin", static = FALSE,
                      title = "Winter minimum temperature (°C)",
                      cols  = c("#1E3A8A", "#3B82F6", "#93C5FD", "#FDE68A", "#F97316")),
  pred_precip  = list(var = "annualPrec", static = FALSE,
                      title = "Annual precipitation (mm)",
                      cols  = c("#EFF6FF", "#BFDBFE", "#60A5FA", "#2563EB", "#1E3A8A")),
  pred_cec     = list(var = "cec30", static = TRUE,
                      title = "Soil cation exchange capacity, 0–30 cm (cmol(c)/kg)",
                      cols  = c("#2E1065", "#6D28D9", "#C026D3", "#FB7185", "#FDE68A")),
  pred_sand    = list(var = "sand30", static = TRUE,
                      title = "Soil sand content, 0–30 cm (%)",
                      cols  = c("#2E1065", "#6D28D9", "#C026D3", "#FB7185", "#FDE68A")),
  pred_north   = list(var = "north_component", static = TRUE, diverging = TRUE,
                      title = "North–south slope index (+ north-facing, − south-facing)",
                      cols  = c("#2E1065", "#6D28D9", "#C026D3", "#FB7185", "#FDE68A")),
  pred_east    = list(var = "east_component", static = TRUE, diverging = TRUE,
                      title = "East–west slope index (+ east-facing, − west-facing)",
                      cols  = c("#2E1065", "#6D28D9", "#C026D3", "#FB7185", "#FDE68A"))
)
pred_title <- function(layer, yr) {
  info <- PRED_INFO[[layer]]
  if (isTRUE(info$static)) info$title else paste0(info$title, ", ", ba_period_text(yr))
}
pred_tag <- function(yr, scen) if (identical(as.character(yr), "2000")) "baseline_2000" else paste0(scen, "_", yr)
pred_common_range <- function(var) {
  ck <- paste0("common_", var)
  if (exists(ck, envir = .ba4_rng, inherits = FALSE)) return(get(ck, envir = .ba4_rng))
  fs <- file.path(.pred_dir, paste0(var, "_", c("baseline_2000", "ssp245_2060", "ssp245_2090",
                                                  "ssp585_2060", "ssp585_2090"), ".tif"))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) return(c(NA_real_, NA_real_))
  q <- do.call(rbind, lapply(fs, function(f) {
    x <- terra::values(terra::rast(f), na.rm = TRUE)
    x <- x[is.finite(x) & x > -9000]                 # -9999 = fill value outside the data
    stats::quantile(x, c(0.01, 0.99), names = FALSE)
  }))
  v <- c(min(q[, 1]), max(q[, 2])); assign(ck, v, envir = .ba4_rng); v
}
load_pred_raster <- function(layer, yr, scen) {
  info <- PRED_INFO[[layer]]
  none <- list(raster = NULL, range = c(NA_real_, NA_real_))
  f <- if (isTRUE(info$static)) file.path(.pred_dir, paste0("static_", info$var, ".tif"))
       else file.path(.pred_dir, paste0(info$var, "_", pred_tag(yr, scen), ".tif"))
  if (!file.exists(f)) return(none)
  ck <- paste0("pred_", f)
  hit <- ba4_cache_get(ck)
  if (!is.null(hit)) return(hit)
  tryCatch({
    r0 <- terra::rast(f)
    r0[r0 < -9000] <- NA                  # -9999 = fill / 'no data' code in these files
    if (isTRUE(info$static)) {            # 1st-99th percentile stretch so outliers don't wash it out
      q <- stats::quantile(terra::values(r0, na.rm = TRUE), c(0.01, 0.99), names = FALSE)
      if (isTRUE(info$diverging)) { m <- max(abs(q)); q <- c(-m, m) }
      r0  <- terra::clamp(r0, q[1], q[2]); rng <- q
    } else {
      rng <- pred_common_range(info$var)
      r0  <- terra::clamp(r0, rng[1], rng[2])
    }
    r   <- terra::project(r0, "EPSG:3857", method = "bilinear")
    res <- list(raster = raster::raster(r), range = rng)
    ba4_cache_put(ck, res); res
  }, error = function(e) none)
}

# ── GEE Canopy Height tile URL ─────────────────────────────────────────────────
# HOW TO REGENERATE (URL expires ~24 hours after generation):
#   1. Go to https://code.earthengine.google.com with qt
#   2. Paste and run initial_chm_script
#   3. Copy the printed URL (decode %7B/%7D → {z}/{x}/{y}) and paste below.

canopy_tile_url <- "https://earthengine.googleapis.com/v1/projects/pbgjam-v2/maps/77ab79e36978d8bde0608e8370a50a02-ca3d1aeb2c0f4747df36e5178ad31bdc/tiles/{z}/{x}/{y}"

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat/2)^2 + cos(lat1*to_rad) * cos(lat2*to_rad) * sin(dlon/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# ── Polygon helpers (base R only; no extra packages) ──────────────────────────
# Vectorised even-odd ray-casting point-in-polygon. Safe for self-intersecting
# user-drawn polygons (unlike s2/GEOS, which can reject invalid rings).
point_in_poly <- function(x, y, px, py) {
  n <- length(px)
  inside <- rep(FALSE, length(x))
  j <- n
  for (i in seq_len(n)) {
    cond <- ((py[i] > y) != (py[j] > y)) &
      (x < (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i]) + px[i])
    cond[is.na(cond)] <- FALSE
    inside <- xor(inside, cond)
    j <- i
  }
  inside
}

# Planar (lon/lat) area-weighted centroid; falls back to the vertex mean for
# degenerate polygons. Accurate enough at the scale of a site polygon.
polygon_centroid <- function(plng, plat) {
  tryCatch({
    cx <- mean(plng); cy <- mean(plat)
    x0 <- plng - cx;  y0 <- plat - cy
    n  <- length(x0)
    nx <- c(x0[-1], x0[1]); ny <- c(y0[-1], y0[1])
    cr <- x0 * ny - nx * y0
    A  <- sum(cr) / 2
    if (!is.finite(A) || abs(A) < 1e-14) stop("degenerate polygon")
    c(cx + sum((x0 + nx) * cr) / (6 * A),
      cy + sum((y0 + ny) * cr) / (6 * A))
  }, error = function(e) c(mean(plng), mean(plat)))
}

site_polygon_vect <- function(site) {
  px <- c(site$poly_lng, site$poly_lng[1])
  py <- c(site$poly_lat, site$poly_lat[1])
  terra::vect(cbind(px, py), type = "polygons", crs = "EPSG:4326")
}

# Species-scoring limits for polygon selections
POLY_MIN_PLOTS   <- 3    # fewer plots than this inside the polygon -> use 20 nearest to centroid
POLY_MAX_SPECIES <- 40   # cap on species scored for a polygon (ranked by summed BA inside it)

.usfs_centroids <- data.frame(
  region = c("R1","R2","R3","R4","R5","R6","R8","R9","R10"),
  lon    = c(-113.5,-105.5,-108.5,-113.0,-119.5,-121.5,-89.0,-84.0,-153.0),
  lat    = c(47.5,40.0,34.5,41.5,37.5,45.5,33.5,42.5,64.0),
  stringsAsFactors = FALSE
)
.usfs_regions_sf <- tryCatch(
  readRDS("usfs_regions.rds"),
  error = function(e) { message("USFS RDS not found - centroid fallback."); NULL }
)

get_usfs_region <- function(lng, lat) {
  region_col <- NULL
  if (!is.null(.usfs_regions_sf)) {
    pt  <- sf::st_sfc(sf::st_point(c(lng, lat)), crs = 4326)
    idx <- suppressMessages(sf::st_within(pt, .usfs_regions_sf, sparse = FALSE))
    if (any(idx)) {
      row <- .usfs_regions_sf[which(idx)[1], ]
      for (cn in c("REGION","REGIONNAME","region","ADMIN_REGI")) {
        if (cn %in% names(row)) { region_col <- as.character(row[[cn]][1]); break }
      }
      if (!is.null(region_col)) {
        region_col <- sub("^0*(\\d+)$","R\\1",region_col)
        region_col <- sub("^Region\\s*(\\d+)$","R\\1",region_col,ignore.case=TRUE)
        if (!region_col %in% .usfs_centroids$region) region_col <- NULL
      }
    }
  }
  if (is.null(region_col)) {
    d_km <- haversine_km(lng, lat, .usfs_centroids$lon, .usfs_centroids$lat)
    region_col <- .usfs_centroids$region[which.min(d_km)]
  }
  region_col
}

# Returns list(price, source, region) so the UI can label a price with the region it
# actually came from. `source` is the price_source of the cell that supplied the price,
# which matters because a fallback cell is a regional average shared by many species.
get_timber_price_info <- function(lng, lat, species_key) {
  none <- list(price = NA_real_, source = NA_character_, region = NA_character_)
  row  <- timber[timber$gjam_species == species_key, ]
  if (nrow(row) == 0) return(none)
  srow <- timber_src[timber_src$gjam_species == species_key, ]
  region <- get_usfs_region(lng, lat)
  price_for_region <- function(reg) {
    if (!reg %in% names(row)) return(NA_real_)
    v <- suppressWarnings(as.numeric(row[[reg]][1]))
    if (is.na(v) || !is.finite(v)) NA_real_ else v
  }
  source_for_region <- function(reg) {
    if (nrow(srow) == 0 || !reg %in% names(srow)) return(NA_character_)
    as.character(srow[[reg]][1])
  }
  p <- price_for_region(region)
  if (!is.na(p)) return(list(price = p, source = source_for_region(region), region = region))
  d_km   <- haversine_km(lng, lat, .usfs_centroids$lon, .usfs_centroids$lat)
  ranked <- .usfs_centroids$region[order(d_km)]
  ranked <- ranked[ranked != region]
  for (reg in ranked) {
    pt <- price_for_region(reg)
    if (!is.na(pt)) return(list(price = pt, source = source_for_region(reg), region = reg))
  }
  none
}

get_timber_price <- function(lng, lat, species_key) {
  get_timber_price_info(lng, lat, species_key)$price
}

# Short human-readable qualifier for a price cell. Aggregate-sourced cells are a regional
# mean shared across many species, so they must not read as a species-specific price.
timber_source_label <- function(src) {
  if (is.na(src)) return(NA_character_)
  switch(src,
         region_hardwood_fallback = "hardwood mean",
         region_softwood_fallback = "softwood mean",
         NA_character_)
}

# One scoring rule for all three components (basal area, wildlife value, timber price):
# divide each species' value by the largest value among the species scored at this site, so the
# best species at the site = 1 and the rest are the same fraction of it. NA stays NA.
norm_site <- function(x) {
  m <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(m) || m <= 0) return(rep(NA_real_, length(x)))
  pmax(x, 0) / m
}

# N-component weighted mean over switched-on components. NA scores are skipped and
# the remaining weights re-normalised, so unrated species are not penalised.
weighted_composite_n <- function(vals, wts, active) {
  vals <- as.matrix(vals)
  n <- nrow(vals)
  score <- rep(NA_real_, n)
  n_components <- integer(n)
  for (i in seq_len(n)) {
    ok <- active & !is.na(vals[i, ])
    n_components[i] <- sum(ok)
    wt_sum <- sum(wts[ok])
    score[i] <- if (!any(ok) || wt_sum == 0) NA_real_ else sum(vals[i, ok] * wts[ok]) / wt_sum
  }
  list(score = score, n_components = n_components)
}

options(tigris_use_cache=TRUE)
state_sf <- tigris::states(cb=TRUE, resolution="500k", progress_bar=FALSE)
state_sf <- state_sf[!state_sf$STUSPS %in% c("AK","HI","PR","VI","GU","MP","AS","DC"),]
state_sf <- sf::st_transform(state_sf, 4326)

# ── Pennsylvania view (all maps open fitted to PA) ─────────────────────────────
# Bounds come from the state polygons already loaded above, not hard-coded values.
.pa_bb <- sf::st_bbox(state_sf[state_sf$STUSPS == "PA", ])
PA_VIEW <- list(
  west  = unname(as.numeric(.pa_bb["xmin"])),
  south = unname(as.numeric(.pa_bb["ymin"])),
  east  = unname(as.numeric(.pa_bb["xmax"])),
  north = unname(as.numeric(.pa_bb["ymax"]))
)
fit_pa <- function(map, padding = c(20, 20)) {
  leaflet::fitBounds(map, PA_VIEW$west, PA_VIEW$south, PA_VIEW$east, PA_VIEW$north,
                     options = list(padding = padding))
}

.pin_svg <- paste0(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 50">',
  '<defs>',
  '<radialGradient id="hg" cx="38%" cy="32%" r="65%">',
  '<stop offset="0%" stop-color="#FF5050"/><stop offset="100%" stop-color="#8B0000"/></radialGradient>',
  '<linearGradient id="sg" x1="0%" y1="0%" x2="100%" y2="0%">',
  '<stop offset="0%" stop-color="#999"/><stop offset="50%" stop-color="#DDD"/><stop offset="100%" stop-color="#999"/></linearGradient>',
  '</defs>',
  '<ellipse cx="16" cy="48" rx="5" ry="2" fill="rgba(0,0,0,0.22)"/>',
  '<rect x="14.5" y="23" width="3" height="23" rx="1.5" fill="url(#sg)"/>',
  '<circle cx="16" cy="16" r="15" fill="url(#hg)"/>',
  '<ellipse cx="11" cy="9" rx="6" ry="4" fill="rgba(255,255,255,0.28)"/>',
  '</svg>'
)
thumbtack_icon <- makeIcon(
  iconUrl=paste0("data:image/svg+xml,",URLencode(.pin_svg,reserved=TRUE)),
  iconWidth=32, iconHeight=50, iconAnchorX=16, iconAnchorY=48
)

# `site` is list(type = "point"|"polygon", lng, lat, poly_lng, poly_lat).
# For a polygon, lng/lat hold the polygon centroid.
get_all_species <- function(site, k = 20) {
  idx <- NULL
  using_polygon <- FALSE
  if (identical(site$type, "polygon") && length(site$poly_lng) >= 3) {
    inside <- which(point_in_poly(coords_lon, coords_lat, site$poly_lng, site$poly_lat))
    if (length(inside) >= POLY_MIN_PLOTS) { idx <- inside; using_polygon <- TRUE }
  }
  if (is.null(idx)) {
    d_km <- haversine_km(site$lng, site$lat, coords_lon, coords_lat)
    idx  <- order(d_km)[seq_len(min(k, nrow(BA_mat)))]
  }
  sub <- BA_mat[idx, , drop = FALSE]
  if (!is.matrix(sub) || typeof(sub) == "list") {
    sub <- matrix(as.numeric(unlist(sub)), nrow = length(idx),
                  dimnames = list(NULL, colnames(BA_mat)))
  }
  tot     <- colSums(sub, na.rm = TRUE)
  present <- names(tot)[tot > 0]
  if (using_polygon && length(present) > POLY_MAX_SPECIES) {
    present <- names(sort(tot[present], decreasing = TRUE))[seq_len(POLY_MAX_SPECIES)]
  }
  data.frame(key          = present,
             display_name = as.character(species_display[present]),
             stringsAsFactors = FALSE)
}

# maps_v4 loader: opens ONE species' map lazily (nothing is read until it is drawn);
# only the min/max for the legend/colour scale is computed, once per file.
load_ba4_raster <- function(key, yr, scen) {
  none <- list(raster = NULL, fallback = FALSE, range = c(NA_real_, NA_real_))
  f <- ba4_path(key, yr, scen)
  if (!file.exists(f)) return(none)
  tryCatch({
    r   <- terra::rast(f)
    cap <- ba4_cap(key)
    rng <- if (!is.na(cap)) c(0, cap) else {      # fallback: this map's own range
      ck <- paste0("rng_", f)
      if (exists(ck, envir = .ba4_rng, inherits = FALSE)) get(ck, envir = .ba4_rng) else {
        g <- terra::global(r, c("min", "max"), na.rm = TRUE)
        v <- c(g[1, 1], g[1, 2]); assign(ck, v, envir = .ba4_rng); v
      }
    }
    list(raster = r, fallback = FALSE, range = rng)
  }, error = function(e) none)
}

all_vis_keys <- ba4_keys   # maps only for the species the model fitted

# After a species is drawn, quietly prepare its other period/scenario maps one at a time while
# the app is idle (reads the file into the disk cache and builds the whole-US copy), so
# switching period or scenario feels instant. A new selection cancels any unfinished preparation.
.ba4_pf <- new.env(parent = emptyenv()); .ba4_pf$gen <- 0L
ba4_prefetch <- function(key, yr, scen) {
  .ba4_pf$gen <- .ba4_pf$gen + 1L; gen <- .ba4_pf$gen
  yr <- as.character(yr)
  cb <- rbind(expand.grid(yr = c("2060", "2090"), scen = c("ssp245", "ssp585"), stringsAsFactors = FALSE),
              data.frame(yr = "2000", scen = "ssp245", stringsAsFactors = FALSE))
  cb <- cb[!(cb$yr == yr & (yr == "2000" | cb$scen == scen)), , drop = FALSE]
  cb <- cb[order((cb$yr != yr) + (cb$scen != scen)), , drop = FALSE]
  step <- function(i) {
    if (i > nrow(cb) || !identical(.ba4_pf$gen, gen)) return(invisible())
    f <- ba4_path(key, cb$yr[i], cb$scen[i])
    if (file.exists(f)) tryCatch({
      invisible(readBin(f, "raw", file.size(f)))
      ba4_view_raster(terra::rast(f), NULL)
    }, error = function(e) NULL)
    later::later(function() step(i + 1L), 0.5)
  }
  later::later(function() step(1L), 1.5)
  invisible()
}
all_vis_display <- sort(species_display[intersect(all_vis_keys, names(species_display))])

# ── Map JavaScript ────────────────────────────────────────────────────────────
# Location map: out-of-bounds popup + bottom-right zoom + polygon tool.
# Polygon tool: a small button above the zoom control toggles draw mode. While
# drawing, clicks add vertices (and do NOT drop a pin: the server ignores map
# clicks while `main_map_drawing` is TRUE). Finish with a double-click, a click
# on the first point, or Enter; cancel with Esc or by clicking the button again.
# The finished ring is sent to Shiny as `main_map_polygon` and drawn by the server.
main_map_js <- r"---(
function(el, x, data) {
  var map = this;
  var inBounds = function(lat, lng) { return !(lat < 24 || lat > 50 || lng < -125 || lng > -66); };
  var outMsg = '<span style="font-family:Lexend,sans-serif;font-size:12px;color:#555;">Select a location within the contiguous United States.</span>';

  L.control.zoom({position: 'bottomright'}).addTo(map);

  var drawing = false, pts = [], dots = [], tempLine = null, rubber = null, btn = null;

  if (!map.getPane('drawPane')) {
    map.createPane('drawPane');
    map.getPane('drawPane').style.zIndex = 450;
  }
  var hint = L.DomUtil.create('div', 'poly-hint', map.getContainer());
  hint.innerHTML = 'Click to add points &middot; double-click or click the first point to finish &middot; Esc to cancel';
  hint.style.display = 'none';

  map.on('click', function(e) {
    if (drawing) return;
    if (!inBounds(e.latlng.lat, e.latlng.lng)) {
      L.popup().setLatLng(e.latlng).setContent(outMsg).openOn(map);
    }
  });

  function clearTemp() {
    if (tempLine) { map.removeLayer(tempLine); tempLine = null; }
    if (rubber)   { map.removeLayer(rubber);   rubber = null; }
    dots.forEach(function(d) { map.removeLayer(d); });
    dots = [];
  }

  function redraw() {
    clearTemp();
    if (!pts.length) return;
    tempLine = L.polyline(pts, {color: '#FF5050', weight: 2, dashArray: '5,5',
                                interactive: false, pane: 'drawPane'}).addTo(map);
    pts.forEach(function(p, i) {
      dots.push(L.circleMarker(p, {radius: i === 0 ? 6 : 4, color: '#FF5050',
                                   fillColor: '#FFFFFF', fillOpacity: 1, weight: 2,
                                   interactive: false, pane: 'drawPane'}).addTo(map));
    });
  }

  function endDraw() {
    drawing = false;
    pts = [];
    clearTemp();
    L.DomUtil.removeClass(btn, 'poly-active');
    L.DomUtil.removeClass(map.getContainer(), 'poly-drawing');
    hint.style.display = 'none';
    // Give the trailing click/dblclick events time to pass before re-enabling
    // double-click zoom and letting the server accept pin clicks again.
    setTimeout(function() {
      map.doubleClickZoom.enable();
      Shiny.setInputValue('main_map_drawing', false);
    }, 400);
  }

  function startDraw() {
    drawing = true;
    pts = [];
    L.DomUtil.addClass(btn, 'poly-active');
    L.DomUtil.addClass(map.getContainer(), 'poly-drawing');
    hint.style.display = 'block';
    map.doubleClickZoom.disable();
    Shiny.setInputValue('main_map_drawing', true);
  }

  function finishDraw() {
    if (pts.length < 3) { endDraw(); return; }
    var lngs = pts.map(function(p) { return ((p.lng + 180) % 360 + 360) % 360 - 180; });
    var lats = pts.map(function(p) { return p.lat; });
    var mlng = lngs.reduce(function(a, b) { return a + b; }, 0) / lngs.length;
    var mlat = lats.reduce(function(a, b) { return a + b; }, 0) / lats.length;
    endDraw();
    if (!inBounds(mlat, mlng)) {
      L.popup().setLatLng([mlat, mlng]).setContent(outMsg).openOn(map);
      return;
    }
    Shiny.setInputValue('main_map_polygon', {lng: lngs, lat: lats}, {priority: 'event'});
  }

  map.on('click', function(e) {
    if (!drawing) return;
    var p = map.latLngToContainerPoint(e.latlng);
    if (pts.length >= 3 && p.distanceTo(map.latLngToContainerPoint(pts[0])) < 10) {
      finishDraw();
      return;
    }
    if (pts.length && p.distanceTo(map.latLngToContainerPoint(pts[pts.length - 1])) < 5) return;
    pts.push(e.latlng);
    redraw();
  });

  map.on('dblclick', function() { if (drawing) finishDraw(); });

  map.on('mousemove', function(e) {
    if (!drawing || !pts.length) return;
    var seg = [pts[pts.length - 1], e.latlng];
    if (!rubber) {
      rubber = L.polyline(seg, {color: '#FF5050', weight: 1.5, opacity: 0.6,
                                interactive: false, pane: 'drawPane'}).addTo(map);
    } else {
      rubber.setLatLngs(seg);
    }
  });

  document.addEventListener('keydown', function(ev) {
    if (!drawing) return;
    if (ev.key === 'Escape') endDraw();
    else if (ev.key === 'Enter') finishDraw();
  });

  var PolyCtl = L.Control.extend({
    options: {position: 'bottomright'},
    onAdd: function() {
      var div = L.DomUtil.create('div', 'leaflet-bar poly-ctl');
      btn = L.DomUtil.create('a', 'poly-btn', div);
      btn.href = '#';
      btn.title = 'Select an area (click points, double-click to finish)';
      btn.setAttribute('role', 'button');
      btn.innerHTML = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round" xmlns="http://www.w3.org/2000/svg"><polygon points="12,3 21,9 18,20 6,20 3,9"/></svg>';
      L.DomEvent.disableClickPropagation(div);
      L.DomEvent.on(btn, 'click', function(ev) {
        L.DomEvent.preventDefault(ev);
        if (drawing) endDraw(); else startDraw();
      });
      return div;
    }
  });
  map.addControl(new PolyCtl());

  // Safety net: once the container has its final size, refit to Pennsylvania.
  setTimeout(function() {
    map.invalidateSize();
    map.fitBounds([[data.south, data.west], [data.north, data.east]], {padding: [30, 30]});
  }, 250);
}
)---"

# Visualizations map: refit to Pennsylvania once the tab's container has its final size.
vis_map_js <- r"---(
function(el, x, data) {
  var map = this;
  setTimeout(function() {
    map.invalidateSize();
    map.fitBounds([[data.south, data.west], [data.north, data.east]], {padding: [20, 20]});
  }, 250);
}
)---"

# ── CSS ───────────────────────────────────────────────────────────────────────
app_css <- "
@import url('https://fonts.googleapis.com/css2?family=Lexend:wght@300;400;500;600;700&display=swap');

:root {
  --forest:   #1C3A28;
  --canopy:   #2A5240;
  --mid:      #3D7A58;
  --sage:     #6BAE88;
  --mist:     #C8DDD4;
  --cream:    #F7F9F6;
  --linen:    #EDF2EE;
  --border:   #C2D5C8;
  --accent:   #3A8A5A;
  --accent-d: #2A6B43;
  --gjam-c:   #2563EB;
  --wild-c:   #DC2626;
  --timb-c:   #CA8A04;
  --comp-c:   #166534;
  --danger:   #7F1D1D;
  --rip:      #0369A1;
  --ink:      #1A2E24;
  --ink2:     #4A6A58;
  --white:    #FFFFFF;
  --radius:   10px;
  --radius-sm:6px;
  --shadow:   0 2px 8px rgba(0,0,0,0.10), 0 1px 2px rgba(0,0,0,0.06);
  --shadow-lg:0 8px 24px rgba(0,0,0,0.14), 0 2px 6px rgba(0,0,0,0.08);
  --font: 'Lexend', sans-serif;
}

* { box-sizing: border-box; }
body { font-family: var(--font); background: var(--forest); color: var(--ink); margin: 0; padding: 0; }

/* ── Header ── */
.app-header {
  background: linear-gradient(135deg, #0F2218 0%, var(--forest) 60%, #1A3D2B 100%);
  padding: 0 28px;
  height: 58px;
  display: flex; align-items: center; justify-content: space-between;
  border-bottom: 1px solid rgba(255,255,255,0.08);
  box-shadow: 0 2px 12px rgba(0,0,0,0.3);
}
.header-left { display: flex; align-items: center; gap: 18px; }
.header-logo img { height: 32px; }
.header-title {
  font-family: var(--font);
  font-size: 15px; font-weight: 600; letter-spacing: 0.02em;
  color: #E0EEE6;
  border-left: 1px solid rgba(255,255,255,0.15);
  padding-left: 18px;
}
.header-sponsors { display: flex; align-items: center; gap: 20px; }
.header-sponsors img { height: 24px; opacity: 0.7; filter: brightness(0) invert(1); }
.header-sponsors .nasa-logo { height: 44px; }
.header-sponsors .pa_dcnr-logo { height: 38px; }

/* ── Tabs ── */
.nav-tabs {
  background: var(--canopy) !important;
  border-bottom: 1px solid rgba(0,0,0,0.25) !important;
  padding: 0 20px; margin: 0 !important;
}
.nav-tabs > li > a {
  font-family: var(--font) !important;
  font-size: 11px; font-weight: 500; letter-spacing: 0.04em; text-transform: uppercase;
  color: var(--mist) !important;
  border: none !important;
  border-radius: var(--radius-sm) var(--radius-sm) 0 0 !important;
  padding: 9px 20px !important; margin: 6px 2px 0 !important;
  background: transparent !important;
  transition: background 0.15s, color 0.15s;
}
.nav-tabs > li.active > a, .nav-tabs > li.active > a:focus {
  color: var(--ink) !important;
  background: var(--cream) !important;
  border: none !important;
  font-weight: 600 !important;
}
.nav-tabs > li > a:hover {
  color: var(--white) !important;
  background: rgba(255,255,255,0.1) !important;
}
.tab-content { padding: 0; }

/* ── Welcome ── */
.welcome-tab {
  height: calc(100vh - 106px);
  background-image: url('forest-background.jpg');
  background-size: cover; background-position: center 40%;
  display: flex; align-items: center; justify-content: center;
}
.welcome-panel {
  background: rgba(6,18,11,0.86);
  padding: 56px 68px;
  max-width: 680px; text-align: center;
  border-radius: 16px;
  border: 1px solid rgba(100,180,130,0.2);
  box-shadow: var(--shadow-lg);
  backdrop-filter: blur(6px);
}
.welcome-title {
  font-family: var(--font);
  font-size: 58px; font-weight: 700;
  color: #E8F5EC; letter-spacing: -0.01em;
  margin: 0 0 20px 0; line-height: 1;
}
.welcome-desc {
  font-family: var(--font);
  font-size: 17px; color: #A4C4AC; line-height: 1.7;
  font-weight: 300; margin: 0;
}

/* ── Location tab ── */
.tab1-layout { display: flex; height: calc(100vh - 106px); }
.sidebar-col {
  width: 300px; min-width: 300px;
  background: var(--canopy);
  padding: 18px 14px; overflow-y: auto;
  border-right: 1px solid rgba(0,0,0,0.2);
}
.sidebar-section-label {
  font-family: var(--font);
  font-size: 9px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--sage); margin: 0 0 10px 0; display: block;
}
.weight-group {
  background: rgba(0,0,0,0.2);
  border: 1px solid rgba(255,255,255,0.07);
  border-radius: var(--radius);
  padding: 12px 14px 10px; margin-bottom: 10px;
  box-shadow: inset 0 1px 3px rgba(0,0,0,0.15);
}
.weight-label {
  font-family: var(--font);
  font-size: 11px; font-weight: 500; color: white !important; margin-bottom: 6px;
}
.weight-note {
  font-family: var(--font);
  font-size: 9px; color: var(--sage); margin-top: 4px; line-height: 1.4;
}

/* Skeuomorphic slider */
.irs--shiny .irs-bar {
  background: linear-gradient(to bottom, #5AAA78, #3A8A5A);
  border-radius: 4px; height: 6px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.3), inset 0 1px 0 rgba(255,255,255,0.2);
}
.irs--shiny .irs-line {
  background: linear-gradient(to bottom, #0C1C16, #1A3028);
  border-radius: 4px; height: 6px;
  box-shadow: inset 0 1px 3px rgba(0,0,0,0.4);
}
.irs--shiny .irs-handle {
  background: linear-gradient(to bottom, #38A8C8, #0E6888) !important;
  border: 1px solid #8AB09A !important;
  border-radius: 50% !important;
  width: 18px !important; height: 18px !important;
  top: 20px !important;
  box-shadow: 0 2px 6px rgba(0,0,0,0.3), inset 0 1px 0 rgba(255,255,255,0.8) !important;
}
.irs--shiny .irs-handle:hover {
  background: linear-gradient(to bottom, #FFFFFF, #D8E8DC) !important;
}
.irs--shiny .irs-single {
  background: var(--accent-d);
  border-radius: var(--radius-sm);
  font-family: var(--font); font-size: 10px; font-weight: 500;
  box-shadow: var(--shadow);
}

/* Wildlife group */
.wg-gjam { border-left: 3px solid var(--gjam-c); }
.wg-timb { border-left: 3px solid var(--timb-c); }
.wildlife-group {
  background: rgba(0,0,0,0.2); border: 1px solid rgba(255,255,255,0.07);
  border-left: 3px solid var(--wild-c); border-radius: var(--radius);
  padding: 12px 14px 4px; margin-bottom: 10px;
  box-shadow: inset 0 1px 3px rgba(0,0,0,0.15);
}
.wildlife-head { font-family: var(--font); font-size: 12px; font-weight: 600; color: white; margin-bottom: 6px; }
.wild-item { padding: 8px 0 2px; border-top: 1px solid rgba(255,255,255,0.07); }
.wild-item-top { display: flex; align-items: center; justify-content: space-between; gap: 8px; }
.wild-item-top .weight-label { margin-bottom: 0; }
.wild-item .form-group, .wild-item .shiny-input-container { margin-bottom: 0; width: auto; }
.wild-item .checkbox { margin: 0; }
.wild-item .checkbox label { padding-left: 0; min-height: 0; }
.wild-item input[type=checkbox] {
  -webkit-appearance: none; appearance: none; position: relative; margin: 0;
  width: 28px; height: 15px; border-radius: 8px; cursor: pointer;
  background: rgba(255,255,255,0.18); transition: background 0.15s;
}
.wild-item input[type=checkbox]::after {
  content: ''; position: absolute; top: 2px; left: 2px; width: 11px; height: 11px;
  border-radius: 50%; background: white; transition: left 0.15s;
}
.wild-item input[type=checkbox]:checked { background: var(--wild-c); }
.wild-item input[type=checkbox]:checked::after { left: 15px; }
.wild-item .checkbox span { display: none; }
.wild-item.wild-off .irs, .wild-item.wild-off .weight-note { opacity: 0.3; pointer-events: none; }
.wild-item.wild-off .weight-label { opacity: 0.55; }
.dm-wild-part { display: flex; justify-content: space-between; gap: 6px;
  font-family: var(--font); font-size: 10px; color: var(--ink); margin-top: 5px; }
.dm-wild-part span:last-child { font-weight: 600; color: var(--wild-c); }
.dm-pest-note {
  margin-top: 8px; padding: 9px 12px; border-radius: var(--radius);
  border: 1px solid var(--border); background: white;
  font-family: var(--font); font-size: 10px; color: var(--ink); line-height: 1.45;
}
.dm-pest-note b { font-size: 9px; letter-spacing: 0.08em; text-transform: uppercase; color: var(--ink2); display: block; margin-bottom: 3px; }
.dm-pest-note .src { color: var(--ink2); font-style: italic; }
.about-inner:not(.about-sources) { padding-bottom: 8px; }
.about-sources { padding-top: 0; }
.about-sources .about-sources-title { font-size: 28px; font-weight: 700; margin: 0 0 10px; }
.about-sources h3 { font-family: var(--font); font-size: 15px; font-weight: 600; margin: 22px 0 6px; color: var(--ink); }
.about-sources li { font-family: var(--font); font-size: 12px; line-height: 1.6; color: var(--ink); }

/* Riparian slider — blue tint */
.riparian-group .irs--shiny .irs-bar {
  background: linear-gradient(to bottom, #38A8C8, #0E6888);
}
.riparian-group .irs--shiny .irs-handle {
  border-color: #38A8C8 !important;
}

.irs--shiny .irs-min,
.irs--shiny .irs-max {
  color: white !important;
  background: rgba(15,34,24,0.95) !important;
  border: 1px solid rgba(100,180,130,0.25) !important;
  border-radius: 4px !important;
  padding: 2px 6px !important;
  box-shadow: inset 0 1px 2px rgba(0,0,0,0.35) !important;
  text-shadow: none !important;
}

.riparian-group {
  background: rgba(0,0,0,0.2);
  border: 1px solid rgba(255,255,255,0.07);
  border-radius: var(--radius);
  padding: 12px 14px 10px; margin-bottom: 10px;
  box-shadow: inset 0 1px 3px rgba(0,0,0,0.15);
}
.riparian-toggle-note {
  font-family: var(--font); font-size: 9px; color: var(--sage);
  margin-top: 4px; line-height: 1.4; display: block;
}

.map-col { flex: 1; position: relative; overflow: hidden; }
.map-col #main_map { height: 100%; width: 100%; }
.coord-display {
  position: absolute; bottom: 18px; left: 18px; z-index: 1000;
  background: rgba(15,34,24,0.92);
  color: var(--mist);
  font-family: var(--font); font-size: 11px; font-weight: 500;
  padding: 8px 16px;
  border-radius: 20px;
  border: 1px solid rgba(100,180,130,0.25);
  box-shadow: var(--shadow);
  pointer-events: none;
}
.map-popups { position: absolute; top: 14px; left: 14px; z-index: 1000; display: flex; gap: 8px; pointer-events: none; }
.map-popup-panel {
  width: 185px; height: 190px;
  border-radius: var(--radius);
  border: 1px solid rgba(100,180,130,0.25);
  background: rgba(15,30,20,0.92);
  display: flex; flex-direction: column; pointer-events: auto;
  box-shadow: var(--shadow-lg);
  overflow: hidden;
}
.map-popup-header {
  padding: 6px 10px;
  background: rgba(20,50,32,0.95);
  border-bottom: 1px solid rgba(100,180,130,0.15);
  flex-shrink: 0;
}
.map-popup-title { font-family: var(--font); font-size: 9px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: var(--sage); }
.map-popup-body { flex: 1; position: relative; overflow: hidden; }

/* ── Polygon selection tool (map button + hint) ── */
.leaflet-bar a.poly-btn {
  display: flex; align-items: center; justify-content: center;
  color: var(--ink); cursor: pointer;
}
.leaflet-bar a.poly-btn svg { width: 16px; height: 16px; }
.leaflet-bar a.poly-btn.poly-active { background: var(--accent); color: var(--white); }
.poly-hint {
  position: absolute; bottom: 18px; left: 50%; transform: translateX(-50%);
  z-index: 1000; pointer-events: none; white-space: nowrap;
  background: rgba(15,34,24,0.92); color: var(--mist);
  font-family: var(--font); font-size: 10px; font-weight: 500;
  padding: 6px 14px; border-radius: 16px;
  border: 1px solid rgba(100,180,130,0.25);
  box-shadow: var(--shadow);
}
.leaflet-container.poly-drawing,
.leaflet-container.poly-drawing .leaflet-interactive { cursor: crosshair !important; }

.flag-tag {
  display: inline-block;
  font-family: var(--font); font-size: 8px; font-weight: 600;
  margin-left: 4px; vertical-align: middle;
  padding: 1px 5px; border-radius: 10px;
}
.flag-tag.flag-riparian { color: var(--rip); background: rgba(3,105,161,0.12); }

/* Time toggle */
.time-toggle {
  display: inline-flex;
  border-radius: 20px;
  background: rgba(0,0,0,0.25);
  padding: 3px;
  border: 1px solid rgba(255,255,255,0.1);
}
.time-seg {
  font-family: var(--font);
  font-size: 10px; font-weight: 500;
  padding: 5px 14px; cursor: pointer; user-select: none;
  color: var(--mist); border-radius: 16px;
  transition: background 0.15s, color 0.15s;
}
.time-seg.active {
  background: linear-gradient(135deg, var(--accent), var(--accent-d));
  color: white;
  box-shadow: 0 1px 4px rgba(0,0,0,0.3);
}

/* dark, colourless basemap for the Visualizations map (tiles only; data layers keep their colours) */
#vis_map.leaflet-container { background: #161616; }
#vis_map .dark-base { filter: invert(1) grayscale(1) brightness(0.62) contrast(0.85); }

/* maps_v4 period / scenario selector */
.ba4-ctl.off { opacity: 0.35; pointer-events: none; filter: grayscale(1); transition: opacity 0.15s; }
.ba4-ctl { display: flex; align-items: center; gap: 22px; margin-left: auto; }
.ba4-group { display: flex; align-items: center; gap: 8px; }
.ba4-label {
  font-family: var(--font); font-size: 10px; font-weight: 700; letter-spacing: 0.08em;
  text-transform: uppercase; color: var(--forest);
}
.ba4-seg-row {
  display: inline-flex; background: var(--white); border: 1.5px solid var(--forest);
  border-radius: 8px; overflow: hidden;
}
.ba4-seg {
  font-family: var(--font); font-size: 12px; font-weight: 600; color: var(--forest);
  padding: 6px 14px; cursor: pointer; user-select: none; white-space: nowrap;
  border-right: 1.5px solid var(--forest); transition: background 0.12s, color 0.12s;
}
.ba4-seg:last-child { border-right: none; }
.ba4-seg:hover { background: var(--linen); }
.ba4-seg.active { background: var(--forest); color: var(--white); }
.ba4-seg.disabled {
  color: #9AA6A0; background: #ECEFED; cursor: not-allowed; pointer-events: none;
}
.ba4-seg.disabled.active { background: #D5DBD8; color: #8A958F; }
.ba4-dark { margin-left: 0; flex-direction: column; align-items: flex-start; gap: 10px; }
.ba4-dark .ba4-label { color: var(--mist); }
.ba4-dark .ba4-seg-row { background: transparent; border-color: rgba(255,255,255,0.55); }
.ba4-dark .ba4-seg { color: #FFFFFF; border-right-color: rgba(255,255,255,0.55); }
.ba4-dark .ba4-seg:hover { background: rgba(255,255,255,0.12); }
.ba4-dark .ba4-seg.active { background: #FFFFFF; color: var(--forest); }
.ba4-dark .ba4-seg.disabled { color: rgba(255,255,255,0.3); background: transparent; }
.ba4-dark .ba4-seg.disabled.active { background: rgba(255,255,255,0.18); color: rgba(255,255,255,0.45); }

/* ── Visualizations tab ── */
.vis-outer { height: calc(100vh - 106px); display: flex; }
.vis-sidebar {
  width: 260px; min-width: 260px;
  background: var(--cream);
  border-right: 1px solid var(--border);
  display: flex; flex-direction: column;
  overflow: hidden;
}
.vis-sidebar-header {
  padding: 10px 12px 8px;
  border-bottom: 1px solid var(--border);
  background: var(--linen); flex-shrink: 0;
}
.vis-sidebar-header h4 {
  font-family: var(--font);
  font-size: 9px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--mid); margin: 0 0 6px 0;
}

/* Species scope toggle (Ranked / All) */
.vis-scope-toggle {
  display: flex;
  background: var(--border);
  border-radius: 20px;
  padding: 2px;
  margin-bottom: 7px;
}
.vis-scope-seg {
  flex: 1; text-align: center;
  font-family: var(--font); font-size: 10px; font-weight: 500;
  color: var(--ink2);
  padding: 4px 0; border-radius: 16px;
  cursor: pointer; user-select: none;
  transition: background 0.15s, color 0.15s;
}
.vis-scope-seg.active {
  background: white;
  color: var(--ink);
  font-weight: 600;
  box-shadow: 0 1px 3px rgba(0,0,0,0.12);
}

.vis-search {
  width: 100%; padding: 7px 10px;
  font-family: var(--font); font-size: 11px;
  border: 1px solid var(--border);
  background: white; color: var(--ink);
  border-radius: var(--radius-sm);
  outline: none; box-shadow: inset 0 1px 2px rgba(0,0,0,0.06);
}
.vis-search:focus { border-color: var(--accent); box-shadow: 0 0 0 3px rgba(58,138,90,0.15); }
.vis-species-list { flex: 1; overflow-y: auto; padding: 4px 0; }
.vis-sp-item {
  padding: 8px 12px;
  font-family: var(--font); font-style: italic; font-size: 11px;
  color: var(--ink); cursor: pointer;
  border-bottom: 1px solid var(--linen);
  transition: background 0.1s;
}
.vis-sp-item:hover { background: var(--linen); }
.vis-sp-item.active { background: var(--mid); color: white; }

/* Ranked badge shown in ranked-mode list items */
.vis-rank-badge {
  display: inline-block;
  font-family: var(--font); font-size: 9px; font-weight: 600; font-style: normal;
  color: var(--comp-c); background: rgba(22,101,52,0.10);
  border-radius: 8px; padding: 1px 6px; margin-right: 6px;
  vertical-align: middle;
}
.vis-sp-item.active .vis-rank-badge { color: rgba(255,255,255,0.8); background: rgba(255,255,255,0.15); }

/* No-site placeholder in ranked mode */
.vis-no-site-msg {
  padding: 20px 12px; text-align: center;
  font-family: var(--font); font-size: 11px; color: var(--ink2);
  line-height: 1.5;
}

.vis-map-area { flex: 1; display: flex; flex-direction: column; min-width: 0; }
.vis-topbar {
  display: flex; align-items: center; gap: 14px; padding: 8px 14px; flex-shrink: 0;
  background: var(--linen); border-bottom: 1px solid var(--border);
}
.vis-species-title {
  font-family: var(--font); font-style: italic; font-size: 13px;
  color: var(--ink); font-weight: 500; flex: 1; min-width: 0;
  white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}
.vis-fallback-note {
  font-family: var(--font); font-size: 9px; color: #92400E;
  background: #FEF3C7; padding: 3px 8px; border-radius: 10px;
  font-style: normal; flex-shrink: 0;
  border: 1px solid #F59E0B;
}
.time-toggle-wrap { display: flex; align-items: center; gap: 8px; }
.time-toggle-label {
  font-family: var(--font);
  font-size: 9px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase;
  color: var(--ink2); white-space: nowrap;
}
.vis-content { flex: 1; display: flex; min-height: 0; }
.vis-map-wrap { flex: 1; position: relative; }
.vis-map-wrap #vis_map { height: 100%; width: 100%; }

/* ── Layer select panel (single-select) ── */
.vis-layer-panel {
  width: 180px; min-width: 180px;
  background: var(--linen);
  border-left: 1px solid var(--border);
  display: flex; flex-direction: column;
  overflow: hidden;
}
.vis-layer-panel-header {
  padding: 10px 12px 8px;
  border-bottom: 1px solid var(--border);
  background: var(--cream); flex-shrink: 0;
}
.vis-layer-panel-header h4 {
  font-family: var(--font);
  font-size: 9px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--mid); margin: 0;
}
.vis-layer-list {
  flex: 1; padding: 10px 10px; display: flex; flex-direction: column; gap: 6px;
}

/* Single-select pill buttons */
.vis-layer-btn {
  display: flex; align-items: center;
  padding: 9px 16px;
  border-radius: 22px;
  border: 1.5px solid var(--border);
  background: white;
  cursor: pointer; user-select: none;
  transition: background 0.15s, border-color 0.15s, box-shadow 0.15s;
  box-shadow: var(--shadow);
}
.vis-layer-btn:hover {
  background: var(--linen);
  border-color: var(--sage);
}
.vis-layer-btn.active {
  background: linear-gradient(135deg, var(--accent), var(--accent-d));
  border-color: var(--accent-d);
  box-shadow: 0 2px 8px rgba(58,138,90,0.35);
}
.vis-layer-btn-label {
  font-family: var(--font); font-size: 11px; font-weight: 500;
  color: var(--ink); flex: 1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
  transition: color 0.15s;
}
.vis-layer-btn.active .vis-layer-btn-label { color: white; font-weight: 600; }

/* Vis map legend */
.vis-legend {
  position: absolute; bottom: 12px; left: 12px; z-index: 1000;
  background: rgba(247,249,246,0.97);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 10px 14px; min-width: 160px;
  pointer-events: none;
  box-shadow: var(--shadow);
}
.legend-title { font-family: var(--font); font-size: 9px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: var(--mid); display: block; margin-bottom: 6px; }
.legend-gradient { height: 10px; width: 100%; border-radius: 5px; }
.legend-ticks { display: flex; justify-content: space-between; margin-top: 4px; }
.legend-ticks span { font-family: var(--font); font-size: 9px; font-weight: 500; color: var(--ink2); }

/* ── Decision-Making tab ── */
.dm-outer { height: calc(100vh - 106px); display: flex; flex-direction: column; background: var(--linen); }
.dm-topbar {
  display: flex; align-items: center; gap: 14px; padding: 10px 18px; flex-shrink: 0;
  background: var(--cream); border-bottom: 1px solid var(--border);
  box-shadow: 0 1px 4px rgba(0,0,0,0.06);
}
.dm-topbar-title {
  font-family: var(--font);
  font-size: 12px; font-weight: 600; letter-spacing: 0.04em;
  color: var(--ink); flex: 1;
}
.dm-loc-text { font-family: var(--font); font-size: 10px; color: var(--ink2); }
.dm-time-controls { display: flex; align-items: center; gap: 10px; }

.dm-main { flex: 1; display: flex; min-height: 0; overflow: hidden; }

/* Left panel — species list */
.dm-left {
  width: 220px; min-width: 220px;
  background: var(--canopy);
  border-right: 1px solid rgba(0,0,0,0.2);
  display: flex; flex-direction: column; overflow: hidden;
}
.dm-left-section {
  padding: 12px 14px 10px;
  border-bottom: 1px solid rgba(255,255,255,0.08);
  flex-shrink: 0;
}
.dm-section-label {
  font-family: var(--font);
  font-size: 9px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--sage); display: block;
}
.dm-species-scroll { flex: 1; overflow-y: auto; }
.dm-sp-row {
  display: flex; align-items: center; gap: 8px;
  padding: 9px 14px;
  border-bottom: 1px solid rgba(255,255,255,0.05);
  cursor: pointer; transition: background 0.1s;
}
.dm-sp-row:hover { background: rgba(255,255,255,0.07); }
.dm-sp-row.selected {
  background: rgba(58,138,90,0.25);
  border-left: 3px solid var(--sage);
  padding-left: 11px;
}
.dm-sp-rank { font-family: var(--font); font-size: 9px; font-weight: 600; color: var(--sage); width: 16px; flex-shrink: 0; }
.dm-sp-name { font-family: var(--font); font-style: italic; font-size: 11px; color: var(--mist); flex: 1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.dm-sp-row.selected .dm-sp-name { color: white; }
.flag-tag.flag-notmod { background: #FEF3C7; color: #92400E; border: 1px solid #F59E0B; }
.dm-sp-score { font-family: var(--font); font-size: 10px; font-weight: 600; color: var(--sage); }
.dm-sp-row.selected .dm-sp-score { color: #A8D8C0; }
.dm-no-data { padding: 24px 14px; font-family: var(--font); font-size: 11px; color: var(--sage); text-align: center; line-height: 1.5; }

/* Center panel — score breakdown */
.dm-center {
  flex: 1; min-width: 0; display: flex; flex-direction: column;
  background: var(--cream); overflow: hidden;
}
.dm-panel-header {
  padding: 10px 16px; flex-shrink: 0;
  background: var(--linen); border-bottom: 1px solid var(--border);
  display: flex; align-items: baseline; gap: 10px;
}
.dm-panel-title {
  font-family: var(--font);
  font-size: 10px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase;
  color: var(--mid);
}
.dm-panel-subtitle { font-family: var(--font); font-size: 11px; font-style: italic; color: var(--ink2); }
.dm-panel-body { flex: 1; overflow-y: auto; padding: 16px; display: flex; flex-direction: column; gap: 14px; }

/* Score cards — 3 components + composite */
.dm-score-layout { display: flex; flex-direction: column; gap: 10px; }

/* Composite score — prominent pill */
.dm-composite-card {
  background: linear-gradient(135deg, var(--comp-c), #0F4A27);
  border-radius: var(--radius);
  padding: 14px 18px;
  display: flex; align-items: center; justify-content: space-between;
  box-shadow: var(--shadow);
}
.dm-composite-label {
  font-family: var(--font); font-size: 10px; font-weight: 600;
  letter-spacing: 0.08em; text-transform: uppercase; color: rgba(255,255,255,0.7);
}
.dm-composite-val {
  font-family: var(--font); font-size: 36px; font-weight: 700;
  color: white; line-height: 1; letter-spacing: -0.02em;
}
.dm-composite-sub {
  font-family: var(--font); font-size: 9px; color: rgba(255,255,255,0.55); margin-top: 2px;
}

/* Three component cards */
.dm-component-row { display: grid; grid-template-columns: repeat(3, 1fr); gap: 8px; }
.dm-comp-card {
  background: white;
  border-radius: var(--radius);
  border: 1px solid var(--border);
  padding: 12px 14px;
  box-shadow: var(--shadow);
  position: relative; overflow: hidden;
}
.dm-comp-card::before {
  content: '';
  position: absolute; top: 0; left: 0; right: 0; height: 3px;
  border-radius: var(--radius) var(--radius) 0 0;
}
.dm-comp-card.gjam-card::before  { background: var(--gjam-c); }
.dm-comp-card.wild-card::before  { background: var(--wild-c); }
.dm-comp-card.timb-card::before  { background: var(--timb-c); }
.dm-comp-label {
  font-family: var(--font); font-size: 9px; font-weight: 600;
  letter-spacing: 0.08em; text-transform: uppercase; color: var(--ink2);
  display: block; margin-bottom: 6px;
}
.dm-comp-val {
  font-family: var(--font); font-size: 26px; font-weight: 700;
  color: var(--ink); line-height: 1;
}
.dm-comp-card.gjam-card  .dm-comp-val { color: var(--gjam-c); }
.dm-comp-card.wild-card  .dm-comp-val { color: var(--wild-c); }
.dm-comp-card.timb-card  .dm-comp-val { color: var(--timb-c); }
.dm-comp-raw { font-family: var(--font); font-size: 9px; color: var(--ink2); margin-top: 4px; line-height: 1.3; }

/* Genus-average badge */
.genus-avg-badge {
  display: inline-block;
  font-family: var(--font); font-size: 11px; font-weight: 700;
  background: #FEF3C7; color: #92400E;
  border: 1px solid #F59E0B;
  border-radius: 8px; padding: 2px 9px; margin-top: 4px;
}

/* Bar charts */
.dm-chart-block {
  background: white;
  border-radius: var(--radius);
  border: 1px solid var(--border);
  padding: 12px 14px;
  box-shadow: var(--shadow);
}
.dm-chart-label {
  font-family: var(--font);
  font-size: 9px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase;
  color: var(--mid); display: block; margin-bottom: 10px;
}
.bar-row { display: flex; align-items: center; gap: 8px; margin-bottom: 6px; }
.bar-row:last-child { margin-bottom: 0; }
.bar-sp-name { font-family: var(--font); font-style: italic; font-size: 10px; color: var(--ink); width: 120px; flex-shrink: 0; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.bar-track { flex: 1; height: 10px; background: var(--linen); overflow: hidden; border-radius: 5px; }
.bar-fill { height: 100%; border-radius: 5px; transition: width 0.3s ease; }
.bar-val { font-family: var(--font); font-size: 9px; font-weight: 600; color: var(--ink2); width: 36px; flex-shrink: 0; text-align: right; }
.bar-row.hl .bar-sp-name { color: var(--comp-c); font-weight: 600; }
.bar-row.hl .bar-val { color: var(--comp-c); }

/* ── About tab ── */
.about-outer { height: calc(100vh - 106px); background: var(--cream); overflow-y: auto; }
.about-hero {
  background: linear-gradient(135deg, #0F2218, var(--canopy));
  padding: 48px 40px; text-align: center;
  border-bottom: 1px solid rgba(0,0,0,0.15);
}
.about-heading { font-family: var(--font); font-size: 30px; font-weight: 700; color: #E8F5EC; margin: 0 0 8px 0; }
.about-subheading { font-family: var(--font); font-size: 14px; color: var(--mist); margin: 0 auto; max-width: 600px; line-height: 1.6; font-weight: 300; }

.about-inner { max-width: 780px; margin: 0 auto; padding: 40px 24px 64px; }

/* Person card: image left, text right */
.about-person {
  display: flex; align-items: flex-start; gap: 24px;
  padding: 28px 0;
  border-bottom: 1px solid var(--border);
}
.about-person:last-child { border-bottom: none; }
.about-photo {
  width: 96px; height: 96px; border-radius: 50%;
  object-fit: cover; flex-shrink: 0;
  border: 3px solid var(--border);
  box-shadow: var(--shadow);
}
.about-person-info { flex: 1; }
.about-name { font-family: var(--font); font-size: 16px; font-weight: 600; color: var(--ink); margin: 0 0 2px 0; }
.about-role {
  font-family: var(--font); font-size: 9px; font-weight: 600;
  letter-spacing: 0.08em; text-transform: uppercase; color: var(--mid);
  margin: 0 0 10px 0;
}
.about-bio { font-family: var(--font); font-size: 12px; color: var(--ink2); line-height: 1.6; margin: 0; }

/* Misc */
.no-data-msg { padding: 24px 14px; text-align: center; color: var(--ink2); font-family: var(--font); font-size: 12px; }
.container-fluid { padding: 0 !important; }
.row { margin: 0 !important; }
.dm-species-scroll::-webkit-scrollbar,
.vis-species-list::-webkit-scrollbar,
.vis-layer-list::-webkit-scrollbar { width: 4px; }
.dm-species-scroll::-webkit-scrollbar-track,
.vis-species-list::-webkit-scrollbar-track,
.vis-layer-list::-webkit-scrollbar-track { background: rgba(0,0,0,0.1); }
.dm-species-scroll::-webkit-scrollbar-thumb,
.vis-species-list::-webkit-scrollbar-thumb,
.vis-layer-list::-webkit-scrollbar-thumb { background: var(--sage); border-radius: 2px; }

.riparian-info {
  background: white; border-radius: var(--radius);
  border: 1px solid var(--border); border-left: 3px solid var(--rip);
  padding: 10px 12px;
}
.riparian-info b {
  display: block; font-family: var(--font); font-size: 9px; font-weight: 600;
  letter-spacing: 0.07em; text-transform: uppercase; color: var(--rip); margin-bottom: 4px;
}
.riparian-info p { margin: 0; font-family: var(--font); font-size: 10px; color: var(--ink); line-height: 1.4; }
.riparian-info .ripar-no { color: var(--ink2); font-style: italic; }

.pest-alert {
  background: linear-gradient(135deg, #7F1D1D, #991B1B);
  border-radius: var(--radius);
  padding: 10px 12px; margin-top: 2px;
  box-shadow: var(--shadow);
}
.pest-alert b { display: block; font-family: var(--font); font-size: 9px; font-weight: 600; letter-spacing: 0.06em; color: #FCA5A5; margin-bottom: 4px; }
.pest-alert p { margin: 0; font-family: var(--font); font-size: 9px; color: #FCA5A5; line-height: 1.4; }
"

score_color_css <- function(norm_val, type="composite") {
  if (is.na(norm_val)) return("#D1D9D4")
  if (type == "gjam")     return(sprintf("rgba(37,99,235,%.2f)",   0.35 + norm_val*0.65))
  if (type == "wildlife") return(sprintf("rgba(220,38,38,%.2f)",   0.35 + norm_val*0.65))
  if (type == "timber")   return(sprintf("rgba(202,138,4,%.2f)",   0.35 + norm_val*0.65))
  if (norm_val >= 0.75) return("#166534")
  if (norm_val >= 0.5)  return("#15803D")
  if (norm_val >= 0.25) return("#16A34A")
  return("#4ADE80")
}

TIE_TOLERANCE <- 0.05

apply_riparian_mode <- function(df, riparian_mode) {
  if (is.null(riparian_mode) || riparian_mode == 0) return(df)
  is_riparian <- wildlife_data$riparian_recommended[match(df$key, wildlife_data$ba_code)]
  is_riparian <- ifelse(is.na(is_riparian), FALSE, is_riparian)
  df$is_riparian <- is_riparian
  if (riparian_mode == 2) return(df[is_riparian, , drop = FALSE])
  ord <- order(-df$score, !df$is_riparian, na.last = TRUE)
  df  <- df[ord, , drop = FALSE]
  if (nrow(df) > 1 && !is.na(df$score[1])) {
    near_top <- which(df$score >= df$score[1] - TIE_TOLERANCE)
    if (length(near_top) > 1) {
      tied_block <- df[near_top, , drop = FALSE]
      tied_block <- tied_block[order(!tied_block$is_riparian), , drop = FALSE]
      df[near_top, ] <- tied_block
    }
  }
  df
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML(app_css)),
    tags$script(HTML("
    
// ── Decision-Making species select ───────────────────────────────────────────
$(document).on('click', '.dm-sp-row', function() {
  $('.dm-sp-row').removeClass('selected'); $(this).addClass('selected');
  Shiny.setInputValue('dm_selected_species', $(this).data('key'), {priority:'event'});
});

// ── Visualizations species select ────────────────────────────────────────────
$(document).on('click', '.vis-sp-item', function() {
  $('.vis-sp-item').removeClass('active'); $(this).addClass('active');
  Shiny.setInputValue('vis_selected_species', $(this).data('key'), {priority:'event'});
});

$(document).on('input', '#vis_search', function() {
  var q = $(this).val().toLowerCase();
  $('.vis-sp-item').each(function() {
    $(this).toggle($(this).text().toLowerCase().indexOf(q) >= 0);
  });
});

// ── Visualizations species scope toggle (Ranked / All) ───────────────────────
$(document).on('click', '.vis-scope-seg', function() {
  $('.vis-scope-seg').removeClass('active');
  $(this).addClass('active');
  Shiny.setInputValue('vis_species_scope', $(this).data('scope'), {priority:'event'});
});

// ── maps_v4 period / scenario selector (all copies stay in sync) ─────────────
$(document).on('click', '.ba4-seg', function() {
  if ($(this).hasClass('disabled')) return;
  var id  = $(this).closest('.ba4-seg-row').data('inputid');
  var val = String($(this).data('val'));
  $('.ba4-seg-row').filter(function() { return $(this).data('inputid') === id; })
    .find('.ba4-seg').each(function() {
      $(this).toggleClass('active', String($(this).data('val')) === val);
    });
  Shiny.setInputValue(id, val, {priority:'event'});
  if (id === 'vis_year') $('.ba4-scen .ba4-seg').toggleClass('disabled', val === '2000');
});

// ── Layer single-select ───────────────────────────────────────────────────────
$(document).on('click', '.vis-layer-btn', function() {
  $('.vis-layer-btn').removeClass('active');
  $(this).addClass('active');
  Shiny.setInputValue('vis_active_layer', $(this).data('layer'), {priority:'event'});
  $('.vis-topbar .ba4-ctl').toggleClass('off', ['pred_cec','pred_sand','pred_north','pred_east','canopy'].indexOf($(this).data('layer')) >= 0);
});

// ── Wildlife part on/off switches ─────────────────────────────────────────────
function syncWildItems() {
  $('.wild-item').each(function() {
    $(this).toggleClass('wild-off', !$(this).find('input[type=checkbox]').prop('checked'));
  });
}
$(document).on('change', '.wild-item input[type=checkbox]', syncWildItems);
$(document).on('shiny:connected', function() { setTimeout(syncWildItems, 50); });

// ── Riparian slider label override ────────────────────────────────────────────
var riparianLabels = ['Disregard', 'Tie-breaker', 'Required'];
function relabelRiparianSlider() {
  $('#riparian_mode').closest('.riparian-group').find('.irs-single, .irs-from, .irs-to').each(function() {
    var v = parseInt($(this).text(), 10);
    if (!isNaN(v) && riparianLabels[v] !== undefined) $(this).text(riparianLabels[v]);
  });
}
$(document).on('shiny:value', function() { setTimeout(relabelRiparianSlider, 50); });
$(document).on('input change', '#riparian_mode', function() { setTimeout(relabelRiparianSlider, 10); });

$(document).ready(function() {
  Shiny.setInputValue('vis_active_layer', 'future_ba', {priority:'event'});
  Shiny.setInputValue('vis_year', '2060', {priority:'event'});
  Shiny.setInputValue('vis_scen', 'ssp245', {priority:'event'});
  // Default species scope to ranked
  Shiny.setInputValue('vis_species_scope', 'ranked', {priority:'event'});
  setTimeout(relabelRiparianSlider, 300);
});
    "))
  ),

div(class="app-header",
    div(class="header-left",
        div(class="header-logo", tags$img(src="logo.png", alt="PBGJAM")),
        div(class="header-title", "PBGJAM v2: Decision Making Tool")
    ),
    div(class="header-sponsors",
        tags$img(src="duke.svg", alt="Duke"),
        tags$img(src="nasa.png", alt="NASA", class="nasa-logo"),
        tags$img(src="pa_dcnr.png", alt="pa_dcnr", class="pa_dcnr-logo"))
),

tabsetPanel(id="main_tabs",
            
            tabPanel("Home",
                     div(class="welcome-tab",
                         div(class="welcome-panel",
                             tags$h1(class="welcome-title", "PBGJAM v2"),
                             tags$p(class="welcome-desc",
                                    "We apply the latest advancements in technology and statistics to ",
                                    "forecast the effects of a changing climate on the abundance and ",
                                    "distribution of America's trees."
                             )
                         )
                     )
            ),
            
            tabPanel("Location",
                     div(class="tab1-layout",
                         div(class="sidebar-col",
                             tags$span(class="sidebar-section-label", "Model Weights"),
                             div(class="weight-group wg-gjam",
                                 div(class="weight-label", "Predicted Future BA Score"),
                                 sliderInput("w1", NULL, min=0, max=1, value=0.34, step=0.01, ticks=FALSE, width="100%"),
                                 div(class="weight-note", "Predicted basal area at the site (selected period and scenario) in m\u00b2/ha; the best species at the site scores 1")
                             ),
                             div(class="wildlife-group",
                                 div(class="wildlife-head", "Wildlife"),
                                 lapply(WILD_COMPONENTS, function(cmp) {
                                   div(class=paste("wild-item", if (!cmp$on) "wild-off" else ""),
                                       div(class="wild-item-top",
                                           div(class="weight-label", cmp$label),
                                           checkboxInput(paste0("wt_on_", cmp$id), NULL, value=cmp$on)),
                                       sliderInput(paste0("wt_", cmp$id), NULL, min=0, max=1,
                                                   value=cmp$w, step=0.01, ticks=FALSE, width="100%"),
                                       div(class="weight-note", cmp$note))
                                 })
                             ),
                             div(class="weight-group wg-timb",
                                 div(class="weight-label", "Timber Market Value"),
                                 sliderInput("w3", NULL, min=0, max=1, value=0.33, step=0.01, ticks=FALSE, width="100%"),
                                 div(class="weight-note", "Regional price from USFS timber sales, 5-yr average (FY2021-FY2025)")
                             ),
                             div(class="riparian-group",
                                 div(class="weight-label", "Riparian Status"),
                                 sliderInput("riparian_mode", NULL, min=0, max=2, value=0, step=1,
                                             ticks=TRUE, width="100%", animate=FALSE),
                                 uiOutput("riparian_mode_note")
                             ),
                             tags$hr(style="border-color:rgba(255,255,255,0.1); margin: 14px 0;"),
                             tags$span(class="sidebar-section-label", "Prediction Period"),
                             ba4_selector(dark = TRUE),
                             tags$hr(style="border-color:rgba(255,255,255,0.1); margin: 14px 0;"),
                             uiOutput("map_sidebar_info")
                         ),
                         div(class="map-col",
                             leafletOutput("main_map", height="100%", width="100%"),
                             div(class="map-popups",
                                 div(class="map-popup-panel",
                                     div(class="map-popup-header", div(class="map-popup-title", "Satellite")),
                                     div(class="map-popup-body", style="padding:0;",
                                         leafletOutput("sat_minimap", height="163px", width="185px"))),
                                 div(class="map-popup-panel",
                                     div(class="map-popup-header", div(class="map-popup-title", "Terrain")),
                                     div(class="map-popup-body", style="padding:0;",
                                         leafletOutput("terrain_minimap", height="163px", width="185px")))
                             ),
                             uiOutput("coord_display")
                         )
                     )
            ),
            
            tabPanel("Decision-Making",
                     div(class="dm-outer",
                         div(class="dm-topbar",
                             div(class="dm-topbar-title", "Species Scoring"),
                             uiOutput("dm_loc_label"),
                             div(class="dm-time-controls",
                                 ba4_selector()
                             ),
                             uiOutput("dm_region_label")
                         ),
                         div(class="dm-main",
                             div(class="dm-left",
                                 div(class="dm-left-section",
                                     tags$span(class="dm-section-label", "Top 10 by composite score")
                                 ),
                                 div(class="dm-species-scroll", uiOutput("dm_species_list"))
                             ),
                             div(class="dm-center",
                                 div(class="dm-panel-header",
                                     div(class="dm-panel-title", "Score Breakdown"),
                                     uiOutput("dm_selected_name")
                                 ),
                                 div(class="dm-panel-body",
                                     uiOutput("dm_score_cards"),
                                     div(class="dm-chart-block",
                                         tags$span(class="dm-chart-label", "By composite score"),
                                         uiOutput("dm_composite_bars")
                                     ),
                                     div(class="dm-chart-block",
                                         tags$span(class="dm-chart-label", "By component"),
                                         uiOutput("dm_component_bars")
                                     )
                                 )
                             )
                         )
                     )
            ),
            
            # ── Visualizations tab ────────────────────────────────────────────────────
            tabPanel("Visualizations",
                     div(class="vis-outer",
                         # Species list sidebar
                         div(class="vis-sidebar",
                             div(class="vis-sidebar-header",
                                 tags$h4("Species"),
                                 # Ranked / All pill toggle — default "Ranked"
                                 div(class="vis-scope-toggle",
                                     div(class="vis-scope-seg active", `data-scope`="ranked", "Ranked"),
                                     div(class="vis-scope-seg",        `data-scope`="all",    "All")
                                 ),
                                 tags$input(id="vis_search", class="vis-search",
                                            type="text", placeholder="Filter species...")
                             ),
                             div(class="vis-species-list", uiOutput("vis_species_list_ui"))
                         ),
                         # Map + layer panel
                         div(class="vis-map-area",
                             div(class="vis-topbar",
                                 uiOutput("vis_species_title_ui"),
                                 ba4_selector()
                             ),
                             div(class="vis-content",
                                 div(class="vis-map-wrap",
                                     leafletOutput("vis_map", height="100%", width="100%"),
                                     div(class="vis-legend", uiOutput("vis_legend_ui"))
                                 ),
                                 # Layer single-select panel
                                 div(class="vis-layer-panel",
                                     div(class="vis-layer-panel-header",
                                         tags$h4("Layer")
                                     ),
                                     div(class="vis-layer-list",
                                         div(class="vis-layer-btn active", `data-layer`="future_ba",
                                             div(class="vis-layer-btn-label", "Predicted Basal Area")
                                         ),
                                         div(class="vis-layer-btn", `data-layer`="pred_wintmin",
                                             div(class="vis-layer-btn-label", "Winter Min Temp")
                                         ),
                                         div(class="vis-layer-btn", `data-layer`="pred_precip",
                                             div(class="vis-layer-btn-label", "Annual Precip")
                                         ),
                                         div(class="vis-layer-btn", `data-layer`="pred_cec",
                                             div(class="vis-layer-btn-label", "Soil CEC")
                                         ),
                                         div(class="vis-layer-btn", `data-layer`="pred_sand",
                                             div(class="vis-layer-btn-label", "Soil Sand")
                                         ),
                                         div(class="vis-layer-btn", `data-layer`="pred_north",
                                             div(class="vis-layer-btn-label", "Slope: North–South")
                                         ),
                                         div(class="vis-layer-btn", `data-layer`="pred_east",
                                             div(class="vis-layer-btn-label", "Slope: East–West")
                                         ),
                                         div(class="vis-layer-btn", `data-layer`="canopy",
                                             div(class="vis-layer-btn-label", "CHM (all species)")
                                         )
                                     )
                                 ) # end vis-layer-panel
                             ) # end vis-content
                         ) # end vis-map-area
                     ) # end vis-outer
            ), # end Visualizations tabPanel
            
            tabPanel("About",
                     div(class="about-outer",
                         div(class="about-hero",
                             div(class="about-heading", "About PBGJAM"),
                             div(class="about-subheading",
                                 "PBGJAM v2 was developed by student Tate Commission, Dr. Tong Qiu, and ",
                                 "Dr. James Clark of Duke University, in collaboration with",
                                 "Jim Julian of the Pennsylvania Department of Natural Resources and",
                                 "Dr. Marc Mcdill of The Pennsylvania State University"
                             )
                         ),
                         div(class="about-inner",
                             div(class="about-person",
                                 tags$img(class="about-photo", src="tatecommission.png", alt="Tate Commission"),
                                 div(class="about-person-info",
                                     div(class="about-name", "Tate Commission"),
                                     div(class="about-role", "Lead Developer"),
                                     tags$p(class="about-bio",
                                            "Tate Commission is an undergraduate at Duke University studying Mechanical Engineering.")
                                 )
                             ),
                             div(class="about-person",
                                 tags$img(class="about-photo", src="tongqiu.png", alt="Dr. Tong Qiu"),
                                 div(class="about-person-info",
                                     div(class="about-name", "Dr. Tong Qiu"),
                                     div(class="about-role", "Co-Principal Investigator"),
                                     tags$p(class="about-bio",
                                            "Dr. Tong Qiu is an Assistant Professor of Ecology at Duke University.")
                                 )
                             ),
                             div(class="about-person",
                                 tags$img(class="about-photo", src="marcmcdill.png", alt="Dr. Marc Mcdill"),
                                 div(class="about-person-info",
                                     div(class="about-name", "Dr. Marc Mcdill"),
                                     div(class="about-role", "Co-Principal Investigator"),
                                     tags$p(class="about-bio",
                                            "Dr. Marc Mcdill is an Associate Professor of Forest Management at The Pennsylvania State University.")
                                 )
                             ),
                             div(class="about-person",
                                 tags$img(class="about-photo", src="jamesclark.png", alt="Dr. James Clark"),
                                 div(class="about-person-info",
                                     div(class="about-name", "Dr. James Clark"),
                                     div(class="about-role", "Project Advisor"),
                                     tags$p(class="about-bio",
                                            "Dr. James Clark is the Nicholas Distinguished Professor of Environmental Science at Duke University.")
                                 )
                             ),
                             div(class="about-person",
                                 tags$img(class="about-photo", src="null.png", alt="Jim Julian"),
                                 div(class="about-person-info",
                                     div(class="about-name", "Jim Julian"),
                                     div(class="about-role", "Application Liason"),
                                     tags$p(class="about-bio",
                                            "Jim Julian is an Ecological Program Specialisty with the Pennsylvania Department of Conservation and Natural Resources Bureau of Forestry.")
                                 )
                             )
                         ),
                         div(class="about-inner about-sources",
                             tags$h3(class="about-sources-title", "Data sources"),
                             tags$ul(
                               tags$li("USDA Forest Service Forest Inventory and Analysis (FIA)"),
                               tags$li("GJAM (generalized joint attribute modeling) species list and predictions"),
                               tags$li("TerraClimate"),
                               tags$li("NASA NEX-GDDP-CMIP6 climate projections (16-model ensemble mean)"),
                               tags$li("SoilGrids (ISRIC) \u2013 soil cation exchange capacity and sand content"),
                               tags$li("USGS 3D Elevation Program \u2013 slope and aspect"),
                               tags$li("Meta / World Resources Institute canopy height"),
                               tags$li("USDA PLANTS Database (via plants.db)"),
                               tags$li("Pennsylvania DCNR Bureau of Forestry Native Tree Selector"),
                               tags$li("Pennsylvania pest and pathogen ratings"),
                               tags$li("US Army Corps of Engineers National Wetland Plant List"),
                               tags$li("USFS timber sale reports, FY2021\u2013FY2025"),
                               tags$li("US Census Bureau state boundaries (tigris)"),
                               tags$li("Base maps: Esri World Imagery, OpenTopoMap, OpenStreetMap")
                             )
                         )
                     )
            )
) # end tabsetPanel
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    click_lng          = NULL,   # pin location, or polygon centroid
    click_lat          = NULL,
    geom_type          = "point",   # "point" | "polygon"
    poly_lng           = NULL,      # polygon vertices (lon/lat), only when geom_type == "polygon"
    poly_lat           = NULL,
    # Visualizations tab — only the active raster is held in memory
    vis_raster         = NULL,
    vis_raster_range   = c(NA_real_, NA_real_),
    vis_fallback       = FALSE,
    usfs_region        = NULL
  )
  
  output$riparian_mode_note <- renderUI({
    mode <- if (is.null(input$riparian_mode)) 0 else input$riparian_mode
    txt <- switch(as.character(mode),
                  "0" = "Disregard riparian status entirely.",
                  "1" = "Use as a tie-breaker among near-equal scores.",
                  "2" = "Only show species recommended for riparian planting.",
                  "Disregard riparian status entirely."
    )
    tags$span(class="riparian-toggle-note", txt)
  })
  
  # The selected site (pin or polygon) in the shape the scoring helpers expect.
  current_site <- reactive({
    if (is.null(rv$click_lng)) return(NULL)
    list(type     = rv$geom_type,
         lng      = rv$click_lng,
         lat      = rv$click_lat,
         poly_lng = rv$poly_lng,
         poly_lat = rv$poly_lat)
  })
  
  wild_state <- reactive({
    do.call(rbind, lapply(WILD_COMPONENTS, function(cmp) {
      on <- input[[paste0("wt_on_", cmp$id)]]
      w  <- input[[paste0("wt_", cmp$id)]]
      if (is.null(on)) on <- cmp$on
      if (is.null(w))  w  <- cmp$w
      data.frame(id=cmp$id, label=cmp$label, col=cmp$col,
                 on=isTRUE(on), w=as.numeric(w), stringsAsFactors=FALSE)
    }))
  })
  
  composite_scores <- reactive({
    site <- current_site()
    if (is.null(site)) return(NULL)
    all_sp <- tryCatch(get_all_species(site), error=function(e) NULL)
    if (is.null(all_sp) || nrow(all_sp) == 0) return(NULL)
    keys <- all_sp$key
    
    yr_sel <- if (!is.null(input$vis_year)) as.character(input$vis_year) else "2060"
    sc_sel <- if (!is.null(input$vis_scen)) input$vis_scen else "ssp245"
    gjam_result        <- get_gjam_vals_v4(site, keys, yr_sel, sc_sel)
    gjam_raw           <- gjam_result$values
    gjam_n             <- norm_site(gjam_raw)   # predicted BA (m2/ha) / best species at the site
    
    # Wildlife parts are fixed equal-step 0-1 scores, used as-is (not site-normalised).
    ws  <- wild_state()
    idx <- match(keys, wildlife_data$ba_code)
    wild_mat <- vapply(ws$col, function(cn) as.numeric(wildlife_data[[cn]][idx]),
                       numeric(length(keys)))
    if (!is.matrix(wild_mat)) wild_mat <- matrix(wild_mat, nrow=length(keys))
    colnames(wild_mat) <- paste0("wp_", ws$id)
    
    # Timber price is regional: for a polygon, the USFS region of its centroid is used.
    timber_info    <- lapply(keys, function(k) get_timber_price_info(site$lng, site$lat, k))
    timber_prices  <- vapply(timber_info, function(x) x$price,  numeric(1))
    timber_sources <- vapply(timber_info, function(x) x$source, character(1))
    timber_n       <- norm_site(timber_prices)
    
    composite <- weighted_composite_n(
      vals   = cbind(gjam_n, wild_mat, timber_n),
      wts    = c(input$w1, ws$w, input$w3),
      active = c(TRUE, ws$on, TRUE)
    )
    
    df <- data.frame(
      key               = keys,
      display_name      = all_sp$display_name,
      gjam              = round(gjam_n, 4),
      timber            = round(timber_n, 4),
      score             = round(composite$score, 4),
      n_components      = composite$n_components,
      timber_raw        = round(timber_prices, 2),
      timber_source     = timber_sources,
      gjam_raw          = unname(gjam_raw),
      gjam_not_modeled  = !ba4_has_map(keys),
      gjam_is_genus_avg = keys %in% gjam_result$fallback_keys,
      stringsAsFactors  = FALSE
    )
    df <- cbind(df, as.data.frame(round(wild_mat, 4)))
    df <- df[order(df$score, decreasing=TRUE, na.last=TRUE), ]
    apply_riparian_mode(df, input$riparian_mode)
  })
  
  top10_scored <- reactive({
    df <- composite_scores()
    if (is.null(df)) return(NULL)
    head(df, 10)
  })
  
  POPUP_ZOOM <- 14
  
  output$main_map <- renderLeaflet({
    leaflet(options=leafletOptions(zoomControl=FALSE, maxBounds=list(c(24,-125),c(50,-66)))) %>%
      addProviderTiles("Esri.WorldImagery", options=providerTileOptions(maxZoom=18)) %>%
      addMapPane("states", zIndex=420) %>%
      addPolygons(data=state_sf, fillOpacity=0, color="white", weight=0.9, opacity=0.55,
                  options=pathOptions(pane="states")) %>%
      fit_pa(padding=c(30, 30)) %>%
      htmlwidgets::onRender(main_map_js, data=PA_VIEW)
  })
  
  # ── Drawing the selected site on every map ───────────────────────────────────
  show_point_on_maps <- function(lng, lat) {
    leafletProxy("main_map") %>% clearMarkers() %>% clearGroup("site_poly") %>%
      addMarkers(lng=lng, lat=lat, icon=thumbtack_icon)
    leafletProxy("vis_map") %>% clearGroup("vis_pin") %>%
      addMarkers(lng=lng, lat=lat, icon=thumbtack_icon, group="vis_pin")
    for (mm in c("sat_minimap", "terrain_minimap")) {
      leafletProxy(mm) %>% clearGroup("site_poly") %>%
        setView(lng=lng, lat=lat, zoom=POPUP_ZOOM) %>% clearMarkers() %>%
        addCircleMarkers(lng=lng, lat=lat, radius=4,
                         color="#FF5050", fillColor="#FF5050", fillOpacity=1, opacity=1,
                         weight=1.5, stroke=TRUE)
    }
  }
  
  show_polygon_on_maps <- function(plng, plat) {
    leafletProxy("main_map") %>% clearMarkers() %>% clearGroup("site_poly") %>%
      addPolygons(lng=plng, lat=plat, group="site_poly",
                  color="#FF5050", weight=2.5, opacity=1,
                  fillColor="#FF5050", fillOpacity=0.18)
    leafletProxy("vis_map") %>% clearGroup("vis_pin") %>%
      addPolygons(lng=plng, lat=plat, group="vis_pin",
                  color="#FF5050", weight=2.5, opacity=1,
                  fillColor="#FF5050", fillOpacity=0.15)
    for (mm in c("sat_minimap", "terrain_minimap")) {
      leafletProxy(mm) %>% clearMarkers() %>% clearGroup("site_poly") %>%
        addPolygons(lng=plng, lat=plat, group="site_poly",
                    color="#FF5050", weight=2, opacity=1,
                    fillColor="#FF5050", fillOpacity=0.12) %>%
        fitBounds(min(plng), min(plat), max(plng), max(plat),
                  options=list(maxZoom=POPUP_ZOOM, padding=c(8, 8)))
    }
  }
  
  # ── Pin drop ─────────────────────────────────────────────────────────────────
  observeEvent(input$main_map_click, {
    if (isTRUE(input$main_map_drawing)) return()   # clicks while drawing a polygon are vertices
    cl  <- input$main_map_click
    lng <- ((cl$lng+180)%%360)-180; lat <- cl$lat
    if (lat < 24 || lat > 50 || lng < -125 || lng > -66) return()
    rv$geom_type   <- "point"
    rv$poly_lng    <- NULL
    rv$poly_lat    <- NULL
    rv$click_lng   <- lng
    rv$click_lat   <- lat
    rv$usfs_region <- get_usfs_region(lng, lat)
    show_point_on_maps(lng, lat)
  })
  
  # ── Polygon selection ────────────────────────────────────────────────────────
  observeEvent(input$main_map_polygon, {
    pg   <- input$main_map_polygon
    plng <- ((as.numeric(unlist(pg$lng)) + 180) %% 360) - 180
    plat <- as.numeric(unlist(pg$lat))
    if (length(plng) < 3 || length(plng) != length(plat) ||
        any(!is.finite(plng)) || any(!is.finite(plat))) return()
    ctr <- polygon_centroid(plng, plat)
    if (ctr[2] < 24 || ctr[2] > 50 || ctr[1] < -125 || ctr[1] > -66) return()
    rv$geom_type   <- "polygon"
    rv$poly_lng    <- plng
    rv$poly_lat    <- plat
    rv$click_lng   <- ctr[1]
    rv$click_lat   <- ctr[2]
    rv$usfs_region <- get_usfs_region(ctr[1], ctr[2])
    show_polygon_on_maps(plng, plat)
  })
  
  output$coord_display <- renderUI({
    if (is.null(rv$click_lat))
      div(class="coord-display", "Click the map to drop a pin, or use the polygon button to select an area")
    else if (identical(rv$geom_type, "polygon"))
      div(class="coord-display",
          sprintf("Area (%d points) \u00b7 centroid %.4f\u00b0N  %.4f\u00b0W",
                  length(rv$poly_lng), rv$click_lat, abs(rv$click_lng)))
    else
      div(class="coord-display", sprintf("%.4f\u00b0N  %.4f\u00b0W", rv$click_lat, abs(rv$click_lng)))
  })
  
  output$map_sidebar_info <- renderUI({
    if (is.null(rv$click_lat)) return(
      div(style="color:var(--sage); font-family:var(--font); font-size:11px; line-height:1.5;",
          "Click the map to select a site, or use the polygon button (bottom right of the map) to outline an area. Top 10 species by composite score appear in the Decision-Making tab.")
    )
    df <- top10_scored()
    if (is.null(df)) return(NULL)
    div(
      tags$span(class="sidebar-section-label", style="margin-top:0;", "Top species at site"),
      tagList(lapply(seq_len(min(5, nrow(df))), function(i) {
        row <- df[i,]
        wd  <- wildlife_data[wildlife_data$ba_code == row$key, ]
        is_riparian <- nrow(wd) > 0 && isTRUE(wd$riparian_recommended[1])
        div(style="display:flex;align-items:center;gap:6px;padding:5px 0;border-bottom:1px solid rgba(255,255,255,0.07);",
            div(style="font-family:var(--font);font-size:9px;font-weight:600;color:var(--sage);width:14px;", i),
            div(style="font-family:var(--font);font-style:italic;font-size:11px;color:var(--mist);flex:1;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;",
                row$display_name,
                if (is_riparian) tags$span(class="flag-tag flag-riparian", "RP") else NULL),
            div(style="font-family:var(--font);font-size:10px;font-weight:600;color:#6BAE88;", sprintf("%.3f",row$score))
        )
      }))
    )
  })
  
  output$sat_minimap <- renderLeaflet({
    leaflet(options=leafletOptions(zoomControl=FALSE, attributionControl=FALSE,
                                   dragging=FALSE, scrollWheelZoom=FALSE, doubleClickZoom=FALSE, boxZoom=FALSE, keyboard=FALSE)) %>%
      addProviderTiles("Esri.WorldImagery", options=providerTileOptions(maxZoom=20)) %>%
      fit_pa(padding=c(4, 4))
  })
  output$terrain_minimap <- renderLeaflet({
    leaflet(options=leafletOptions(zoomControl=FALSE, attributionControl=FALSE,
                                   dragging=FALSE, scrollWheelZoom=FALSE, doubleClickZoom=FALSE, boxZoom=FALSE, keyboard=FALSE)) %>%
      addProviderTiles("OpenTopoMap", options=providerTileOptions(maxZoom=17)) %>%
      fit_pa(padding=c(4, 4))
  })
  
  # ── Visualizations tab map ──────────────────────────────────────────────────
  output$vis_map <- renderLeaflet({
    leaflet(options=leafletOptions(zoomControl=TRUE)) %>%
      addProviderTiles("OpenStreetMap.Mapnik", options=providerTileOptions(maxZoom=19, className="dark-base")) %>%
      addPolygons(data=state_sf, fillOpacity=0, color="#9FB8AA",
                  weight=0.8, opacity=0.5) %>%
      fit_pa(padding=c(20, 20)) %>%
      addScaleBar(position="bottomright", options=scaleBarOptions(metric=TRUE, imperial=FALSE)) %>%
      htmlwidgets::onRender(vis_map_js, data=PA_VIEW)
  })
  
  # ── Species list respects Ranked / All scope ───────────────────────────────
  output$vis_species_list_ui <- renderUI({
    scope     <- if (is.null(input$vis_species_scope)) "ranked" else input$vis_species_scope
    sel       <- input$vis_selected_species
    ranked_df <- top10_scored()
    
    if (scope == "ranked") {
      # Show only ranked species (top 10 at selected site)
      if (is.null(ranked_df) || nrow(ranked_df) == 0) {
        return(div(class="vis-no-site-msg",
                   "Select a site on the Location tab to see ranked species."))
      }
      tagList(lapply(seq_len(nrow(ranked_df)), function(i) {
        row <- ranked_df[i, ]
        key <- row$key
        nm  <- row$display_name
        cls <- if (!is.null(sel) && sel == key) "vis-sp-item active" else "vis-sp-item"
        div(class=cls, `data-key`=key,
            tags$span(class="vis-rank-badge", i),
            nm)
      }))
    } else {
      # Show all species (original behaviour)
      keys_available <- names(all_vis_display)
      tagList(lapply(seq_along(all_vis_display), function(i) {
        key <- keys_available[i]; nm <- all_vis_display[i]
        cls <- if (!is.null(sel) && sel == key) "vis-sp-item active" else "vis-sp-item"
        div(class=cls, `data-key`=key, nm)
      }))
    }
  })
  
  output$vis_species_title_ui <- renderUI({
    layer <- if (!is.null(input$vis_active_layer)) input$vis_active_layer else "future_ba"
    key   <- input$vis_selected_species
    if (layer == "canopy") {
      return(div(class="vis-species-title",
                 style="font-style:normal;color:var(--ink2);font-size:11px;",
                 "CHM (all species) — select a BA layer to filter by species"))
    }
    if (layer %in% names(PRED_INFO))
      return(div(class="vis-species-title", style="font-style:normal;", pred_title(layer, input$vis_year)))
    if (is.null(key) || !key %in% names(species_display))
      return(div(class="vis-species-title",
                 style="font-style:normal;color:var(--ink2);font-size:11px;",
                 "Select a species to view its basal area maps"))
    div(class="vis-species-title", species_display[key])
  })
  
  # CHM (canopy) layer: species-independent
  observeEvent(input$vis_active_layer, {
    layer <- input$vis_active_layer
    if (!is.null(layer) && layer == "canopy") {
      rv$vis_raster       <- NULL
      rv$vis_fallback     <- FALSE
      rv$vis_raster_range <- c(NA, NA)
      if (identical(input$main_tabs, "Visualizations"))
        session$onFlushed(function() { render_vis_layers() }, once=TRUE)
    }
  }, ignoreNULL=TRUE, ignoreInit=FALSE)
  
  # Basal-area and predictor layers
  observeEvent(list(input$vis_selected_species, input$vis_year, input$vis_scen, input$vis_active_layer), {
    key    <- input$vis_selected_species
    yr     <- if (!is.null(input$vis_year)) as.character(input$vis_year) else "2060"
    sc     <- if (!is.null(input$vis_scen)) input$vis_scen else "ssp245"
    layer  <- if (!is.null(input$vis_active_layer)) input$vis_active_layer else "future_ba"
    
    if (layer == "canopy") return()
    if (layer %in% names(PRED_INFO)) {
      result <- load_pred_raster(layer, yr, sc)
      rv$vis_raster       <- result$raster
      rv$vis_fallback     <- FALSE
      rv$vis_raster_range <- result$range
      return()
    }
    if (is.null(key) || key == "") return()
    
    if (layer == "future_ba") {
      result <- load_ba4_raster(key, yr, sc)
      rv$vis_raster       <- result$raster
      rv$vis_fallback     <- isTRUE(result$fallback)
      rv$vis_raster_range <- result$range
      if (!is.null(result$raster)) ba4_prefetch(key, yr, sc)   # warm this species' other periods
    }
  }, ignoreNULL=FALSE)
  
  # Green palette for BA layers
  orange_pal <- function(rng) {
    colorNumeric(
      palette  = c("#0B3B2E","#17803D","#2FD16B","#B8F25A","#FFF07A"),
      domain   = rng,
      na.color = "transparent"
    )
  }
  
  VIS_RASTER_ID <- "vis_ba_raster"
  VIS_TILE_ID   <- "vis_canopy_tile"
  VIS_MAX_BYTES <- 32L * 1024L * 1024L
  ba4_last      <- new.env()   # view the BA image was last drawn for
  draw_ba_overlay <- function() {
    r   <- isolate(rv$vis_raster)
    rng <- isolate(rv$vis_raster_range)
    if (is.null(r) || any(is.na(rng))) return(invisible())
    b    <- isolate(input$vis_map_bounds)
    clip <- tryCatch(ba4_view_raster(r, b),
                     error = function(e) { message("BA draw failed: ", conditionMessage(e)); NULL })
    ba4_last$b <- b
    if (is.null(clip)) return(invisible())
    if (rng[2] <= rng[1]) rng[2] <- rng[1] + 1e-6
    leafletProxy("vis_map") %>% addRasterImage(
      clip, colors = ba4_clamped_pal(orange_pal(rng), rng), opacity = 0.78, project = FALSE,
      layerId = VIS_RASTER_ID, maxBytes = VIS_MAX_BYTES)   # same id = replaces the old image
    invisible()
  }
  
  render_vis_layers <- function() {
    layer <- isolate(input$vis_active_layer)
    if (is.null(layer)) layer <- "future_ba"
    
    proxy <- leafletProxy("vis_map")
    
    proxy <- proxy %>%
      removeImage(VIS_RASTER_ID) %>%
      removeTiles(VIS_TILE_ID)
    
    # Re-draw the selected site (pin or polygon) on this map
    lng <- isolate(rv$click_lng)
    lat <- isolate(rv$click_lat)
    if (!is.null(lng)) {
      proxy <- proxy %>% clearGroup("vis_pin")
      if (identical(isolate(rv$geom_type), "polygon") && !is.null(isolate(rv$poly_lng))) {
        proxy <- proxy %>%
          addPolygons(lng=isolate(rv$poly_lng), lat=isolate(rv$poly_lat), group="vis_pin",
                      color="#FF5050", weight=2.5, opacity=1,
                      fillColor="#FF5050", fillOpacity=0.15)
      } else {
        proxy <- proxy %>%
          addMarkers(lng=lng, lat=lat, icon=thumbtack_icon, group="vis_pin")
      }
    }
    
    if (layer == "future_ba") {
      r   <- isolate(rv$vis_raster)
      rng <- isolate(rv$vis_raster_range)
      if (!is.null(r) && !any(is.na(rng))) draw_ba_overlay()
    } else if (layer %in% names(PRED_INFO)) {
      r   <- isolate(rv$vis_raster)
      rng <- isolate(rv$vis_raster_range)
      if (!is.null(r) && !any(is.na(rng))) {
        if (rng[2] <= rng[1]) rng[2] <- rng[1] + 1e-6
        pal <- ba4_clamped_pal(colorNumeric(PRED_INFO[[layer]]$cols, domain = rng, na.color = "transparent"), rng)   # tiny float overshoot after reprojection is clamped
        proxy <- proxy %>% addRasterImage(
          r, colors = pal, opacity = 0.85, project = FALSE,
          layerId = VIS_RASTER_ID, maxBytes = VIS_MAX_BYTES)
      }
    } else if (layer == "canopy") {
      if (!is.null(canopy_tile_url) && nzchar(canopy_tile_url)) {
        proxy <- proxy %>%
          addTiles(
            urlTemplate = canopy_tile_url,
            layerId     = VIS_TILE_ID,
            options     = tileOptions(opacity=0.78, maxZoom=18),
            attribution = "Meta/WRI Canopy Height \u2022 GEE Community Catalog"
          )
      }
    }
  }
  
  observeEvent(input$main_tabs, {
    if (identical(input$main_tabs, "Visualizations"))
      session$onFlushed(function() { render_vis_layers() }, once=TRUE)
  })
  
  observeEvent(rv$vis_raster, {
    if (identical(input$main_tabs, "Visualizations")) render_vis_layers()
  }, ignoreNULL=FALSE)
  
  observeEvent(input$vis_active_layer, {
    if (identical(input$main_tabs, "Visualizations"))
      session$onFlushed(function() { render_vis_layers() }, once=TRUE)
  }, ignoreNULL=TRUE)
  
  # Redraw the BA image when the user pans/zooms outside what was drawn
  vis_bounds_d <- debounce(reactive(input$vis_map_bounds), 400)
  observeEvent(vis_bounds_d(), {
    if (!identical(input$main_tabs, "Visualizations")) return()
    if (!isTRUE(input$vis_active_layer %in% "future_ba")) return()
    if (is.null(rv$vis_raster)) return()
    b <- vis_bounds_d(); lb <- ba4_last$b
    if (!is.null(b) && !is.null(lb)) {
      w  <- b$east - b$west; lw <- lb$east - lb$west
      pw <- lw * 0.25;       ph <- (lb$north - lb$south) * 0.25
      inside <- b$west >= lb$west - pw && b$east <= lb$east + pw &&
                b$south >= lb$south - ph && b$north <= lb$north + ph
      if (inside && w <= lw * 1.001 && w >= lw * 0.7) return()
    }
    draw_ba_overlay()
  }, ignoreNULL = TRUE)

  # ── Legend ────────────────────────────────────────────────────────────────────
  output$vis_legend_ui <- renderUI({
    layer <- input$vis_active_layer
    if (is.null(layer)) layer <- "future_ba"
    
    if (layer == "canopy") {
      return(tagList(
        tags$span(class="legend-title", "CHM — all species (m)"),
        div(class="legend-gradient",
            style="background:linear-gradient(to right,#0A1F0A,#2E6B2E,#7DB87D,#E8F4D0,#FFFFFF);"),
        div(class="legend-ticks",
            tags$span("0"), tags$span("7"), tags$span("15"), tags$span("22"), tags$span("30+"))
      ))
    }
    
    rng <- rv$vis_raster_range
    if (is.null(rng) || any(is.na(rng))) return(NULL)
    
    if (layer %in% names(PRED_INFO)) {
      ticks <- seq(rng[1], rng[2], length.out = 5)
      return(tagList(
        tags$span(class="legend-title", pred_title(layer, input$vis_year)),
        div(class="legend-gradient",
            style=paste0("background:linear-gradient(to right,", paste(PRED_INFO[[layer]]$cols, collapse=","), ");")),
        div(class="legend-ticks", lapply(formatC(ticks, format="g", digits=3), tags$span))
      ))
    }
    
    title <- ba4_title(input$vis_year, input$vis_scen)
    
    ticks <- pretty(rng, n=4)
    ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
    if (length(ticks) < 2) ticks <- round(c(rng[1], mean(rng), rng[2]), 2)
    tick_labels <- formatC(ticks, format="f", digits=2)
    
    tagList(
      tags$span(class="legend-title", title),
      div(class="legend-gradient",
          style="background:linear-gradient(to right,#0B3B2E,#17803D,#2FD16B,#B8F25A,#FFF07A);"),
      div(class="legend-ticks",
          lapply(tick_labels, tags$span))
    )
  })
  
  # ── Decision-Making tab ─────────────────────────────────────────────────────
  output$dm_loc_label <- renderUI({
    if (is.null(rv$click_lat))
      div(class="dm-loc-text", "No site selected — click the map on the Location tab.")
    else if (identical(rv$geom_type, "polygon"))
      div(class="dm-loc-text",
          sprintf("Area (%d points) \u00b7 centroid %.4f\u00b0N  %.4f\u00b0W",
                  length(rv$poly_lng), rv$click_lat, abs(rv$click_lng)))
    else
      div(class="dm-loc-text", sprintf("%.4f\u00b0N  %.4f\u00b0W", rv$click_lat, abs(rv$click_lng)))
  })
  
  output$dm_region_label <- renderUI({
    if (!is.null(rv$usfs_region))
      div(class="dm-loc-text", paste0("USFS ", rv$usfs_region))
    else NULL
  })
  
  output$dm_species_list <- renderUI({
    df  <- top10_scored()
    sel <- input$dm_selected_species
    if (is.null(df))
      return(div(class="dm-no-data", "Select a site on the Location tab to see species rankings."))
    
    items <- lapply(seq_len(nrow(df)), function(i) {
      row <- df[i,]
      cls <- if (!is.null(sel) && identical(row$key, sel)) "dm-sp-row selected" else "dm-sp-row"
      wd  <- wildlife_data[wildlife_data$ba_code == row$key, ]
      riparian_flag    <- nrow(wd) > 0 && isTRUE(wd$riparian_recommended[1])
      riparian_mode_on <- !is.null(input$riparian_mode) && input$riparian_mode > 0
      div(class=cls, `data-key`=row$key,
          div(class="dm-sp-rank", i),
          div(class="dm-sp-name", row$display_name,
              if (riparian_mode_on && riparian_flag)
                tags$span(class="flag-tag flag-riparian", "RP") else NULL,
              if (isTRUE(row$gjam_not_modeled)) tags$span(class="flag-tag flag-notmod", "Not modeled") else NULL),
          div(class="dm-sp-score", sprintf("%.3f", row$score))
      )
    })
    
    if (is.null(sel) || !(sel %in% df$key)) {
      top_key <- df$key[1]
      items <- c(items, list(tags$script(HTML(sprintf(
        "(function(){Shiny.setInputValue('dm_selected_species','%s',{priority:'event'});})();",
        top_key
      )))))
    }
    tagList(items)
  })
  
  output$dm_selected_name <- renderUI({
    key <- input$dm_selected_species
    if (is.null(key) || !key %in% names(species_display)) return(div(class="dm-panel-subtitle", "—"))
    div(class="dm-panel-subtitle", species_display[key])
  })
  
  output$dm_score_cards <- renderUI({
    df  <- top10_scored()
    key <- input$dm_selected_species
    fmt <- function(x) if (is.na(x)) "—" else sprintf("%.3f", x)
    ws_on <- wild_state()
    ws_on <- ws_on[ws_on$on, , drop = FALSE]
    
    if (is.null(df) || is.null(key)) {
      return(div(class="dm-score-layout",
                 div(class="dm-composite-card",
                     div(class="dm-composite-label", "Composite Score"),
                     div(class="dm-composite-val", "—")),
                 div(class="dm-component-row",
                     div(class="dm-comp-card gjam-card",
                         tags$span(class="dm-comp-label","Predicted Future BA"),
                         div(class="dm-comp-val","—")),
                     div(class="dm-comp-card wild-card",
                         tags$span(class="dm-comp-label","Wildlife"),
                         div(class="dm-comp-val","—")),
                     div(class="dm-comp-card timb-card",
                         tags$span(class="dm-comp-label","Timber Value"),
                         div(class="dm-comp-val","—"))
                 )
      ))
    }
    
    row <- df[df$key == key, ]
    if (nrow(row) == 0) return(NULL)
    
    mode <- if (is.null(input$riparian_mode)) 0 else input$riparian_mode
    comp_sub <- if (mode == 2) {
      "riparian-required mode"
    } else if (!is.na(row$n_components[1]) && row$n_components[1] < 2 + nrow(ws_on)) {
      paste0("based on ", row$n_components[1], " of ", 2 + nrow(ws_on), " components")
    } else paste0("weighted composite of ", 2 + nrow(ws_on), " components")
    
    gjam_raw_note <- tagList(
      if (isTRUE(row$gjam_not_modeled[1])) div(class="genus-avg-badge", "Not modeled"),
      div(class="dm-comp-raw",
        if (is.na(row$gjam_raw[1])) "no data"
        else paste0(formatC(row$gjam_raw[1], format="fg", digits=3, flag="#"), " m\u00b2/ha",
                    if (isTRUE(row$gjam_is_genus_avg[1])) ", genus average" else "")))
    
    wild_card_body <- if (nrow(ws_on) == 0) {
      div(class="dm-comp-raw", "All wildlife parts are off")
    } else {
      tagList(lapply(seq_len(nrow(ws_on)), function(j) {
        v <- row[[paste0("wp_", ws_on$id[j])]][1]
        div(class="dm-wild-part", span(ws_on$label[j]),
            span(if (is.na(v)) "not rated" else sprintf("%.2f", v)))
      }))
    }
    
    wd <- wildlife_data[wildlife_data$ba_code == key, ]
    pest_txt  <- if (nrow(wd) > 0) wd$pest_text_display[1] else NA_character_
    pest_src  <- if (nrow(wd) > 0) wd$pest_text_source[1] else NA_character_
    pest_block <- div(class="dm-pest-note",
        tags$b("Pests and pathogens"),
        if (!is.na(pest_txt) && nzchar(pest_txt)) div(pest_txt) else div("Not rated"),
        if (!is.na(pest_src) && nzchar(pest_src)) div(class="src", pest_src) else NULL)
    
    div(class="dm-score-layout",
        div(class="dm-composite-card",
            div(div(class="dm-composite-label", "Composite Score"),
                div(class="dm-composite-sub", comp_sub)),
            div(class="dm-composite-val", fmt(row$score[1]))
        ),
        div(class="dm-component-row",
            div(class="dm-comp-card gjam-card",
                tags$span(class="dm-comp-label", "Predicted Future BA"),
                div(class="dm-comp-val", fmt(row$gjam[1])),
                gjam_raw_note
            ),
            div(class="dm-comp-card wild-card",
                tags$span(class="dm-comp-label", "Wildlife"),
                wild_card_body
            ),
            div(class="dm-comp-card timb-card",
                tags$span(class="dm-comp-label", "Timber Value"),
                div(class="dm-comp-val", fmt(row$timber[1])),
                div(class="dm-comp-raw",
                    if (!is.na(row$timber_raw[1])) {
                      lbl <- timber_source_label(row$timber_source[1])
                      paste0("$", formatC(row$timber_raw[1], format="f", digits=0), "/MBF",
                             if (!is.na(lbl)) paste0(" (", lbl, ")") else "")
                    } else "no price data")
            )
        ),
        pest_block
    )
  })
  
  output$dm_composite_bars <- renderUI({
    df  <- top10_scored()
    sel <- input$dm_selected_species
    if (is.null(df)) return(div(class="no-data-msg", "—"))
    tagList(lapply(seq_len(nrow(df)), function(i) {
      row    <- df[i,]
      is_sel <- !is.null(sel) && identical(row$key, sel)
      pct    <- if (!is.na(row$score)) row$score * 100 else 0
      col    <- score_color_css(row$score, "composite")
      div(class=if(is_sel) "bar-row hl" else "bar-row",
          div(class="bar-sp-name", row$display_name),
          div(class="bar-track",
              div(class="bar-fill", style=sprintf("width:%.1f%%;background:%s;", pct, col))),
          div(class="bar-val", if (!is.na(row$score)) sprintf("%.3f",row$score) else "—")
      )
    }))
  })
  
  output$dm_component_bars <- renderUI({
    df  <- top10_scored()
    sel <- input$dm_selected_species
    if (is.null(df)) return(div(class="no-data-msg", "—"))
    
    ws_on <- wild_state()
    ws_on <- ws_on[ws_on$on, , drop = FALSE]
    header <- div(style="display:flex;gap:12px;margin-bottom:8px;",
                  div(style="width:120px;flex-shrink:0;"),
                  div(style="flex:1;display:flex;gap:4px;",
                      div(style="flex:1;font-family:var(--font);font-size:8px;font-weight:600;letter-spacing:0.06em;text-transform:uppercase;color:var(--gjam-c);text-align:center;","Pred. BA"),
                      div(style="flex:1;font-family:var(--font);font-size:8px;font-weight:600;letter-spacing:0.06em;text-transform:uppercase;color:var(--wild-c);text-align:center;","Wildlife"),
                      div(style="flex:1;font-family:var(--font);font-size:8px;font-weight:600;letter-spacing:0.06em;text-transform:uppercase;color:var(--timb-c);text-align:center;","Timber")
                  )
    )
    rows <- lapply(seq_len(nrow(df)), function(i) {
      row    <- df[i,]
      is_sel <- !is.null(sel) && identical(row$key, sel)
      g_pct  <- if (!is.na(row$gjam))    row$gjam*100    else 0
      wild_bars <- lapply(seq_len(nrow(ws_on)), function(j) {
        v <- row[[paste0("wp_", ws_on$id[j])]][1]
        div(class="bar-track", style="height:5px;margin:1px 0;flex:none;", title=ws_on$label[j],
            div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--wild-c);opacity:0.75;",
                                                if (is.na(v)) 0 else v*100)))
      })
      t_pct  <- if (!is.na(row$timber))  row$timber*100  else 0
      div(class=if(is_sel) "bar-row hl" else "bar-row", style="align-items:center;",
          div(class="bar-sp-name", row$display_name),
          div(style="flex:1;display:flex;gap:4px;",
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--gjam-c);opacity:0.75;", g_pct))),
              div(style="flex:1;display:flex;flex-direction:column;justify-content:center;", wild_bars),
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--timb-c);opacity:0.75;", t_pct)))
          )
      )
    })
    tagList(header, rows)
  })
  
}

shinyApp(ui, server)
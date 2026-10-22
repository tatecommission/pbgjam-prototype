# =============================================================================
# PBGJAM v2 — Decision Making Tool (Shiny)
#
# CHANGES IN THIS VERSION (everything else is unchanged):
#   1. Polygon selection tool added next to the pin drop. A small button sits
#      above the zoom control on the Location map (bottom-right). Click it,
#      click points on the map, then double-click (or click the first point, or
#      press Enter) to finish; Esc cancels. Clicking the map normally drops a
#      pin as before, and a new pin clears any polygon (and vice versa).
#        - Species list: BA plots INSIDE the polygon (>= POLY_MIN_PLOTS); if the
#          polygon contains too few plots, falls back to the 20 nearest plots
#          to the polygon centroid (same as the pin). At most POLY_MAX_SPECIES
#          species (ranked by summed BA inside the polygon) are scored.
#        - Predicted BA: mean of the raster cells inside the polygon (falls back
#          to the centroid cell if the polygon is smaller than a cell).
#        - Timber price: USFS region of the polygon centroid.
#   2. All four maps (Location, Satellite/Terrain minimaps, Visualizations)
#      now open fitted to Pennsylvania. PA bounds are computed from the
#      tigris state polygons, not hard-coded.
#   3. Raster extraction now projects the site geometry into each raster's
#      CRS first, so it keeps working when the new EPSG:5070 GEE/DCC rasters
#      replace the current ones.
#   4. WILDLIFE is now an umbrella group of separate sliders (no single
#      "Wildlife Value" slider). Each part has its own on/off switch and weight:
#      overall wildlife value, deer forage, deer resistance, fruit/seed food,
#      pollinators. Switched-off parts are ignored entirely. Species with no
#      rating for a part are not penalised: weights are re-normalised over the
#      parts that have data (same rule as before).
#   5. Wildlife / pest / riparian data now come from
#      BA_wildlife_traits_v2.csv (v2-wildlife-pests-workflow). Scores are on
#      equal-step 0-1 scales already, so they are used as-is.
#   6. About tab: data sources and method notes added.
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

# ── Wildlife / pest / riparian table (v2 workflow output) ─────────────────────
WILDLIFE_FILE <- "BA_wildlife_traits_v2.csv"
# Army Corps wetland-region column used for riparian status. Must match
# wetland_region in v2-wildlife-pests-workflow/00_config.R.
WETLAND_COL   <- "wl_ncne"

# One entry per wildlife slider. `col` is a 0-1 equal-step score in the CSV.
WILD_COMPONENTS <- list(
  list(id="overall",   col="wild_wildlife_overall",
       label="Overall wildlife value",
       note="PA DCNR Native Tree Selector wildlife rating",
       on=TRUE,  w=0.11),
  list(id="food",      col="wild_food_seeds_fruit",
       label="Fruit, nut & seed food",
       note="USDA PLANTS fruit/seed abundance and product",
       on=TRUE,  w=0.11),
  list(id="pollin",    col="wild_pollinators",
       label="Pollinators",
       note="PA tool: bees, butterflies, hummingbirds, other insects",
       on=TRUE,  w=0.11),
  list(id="deerfor",   col="wild_deer_forage",
       label="Deer forage (food for deer)",
       note="Higher = more palatable to deer. Opposite of deer resistance",
       on=FALSE, w=0.11),
  list(id="deerres",   col="wild_deer_resistance",
       label="Deer resistance",
       note="Higher = less browsed by deer. Opposite of deer forage",
       on=FALSE, w=0.11)
)
WILD_IDS <- vapply(WILD_COMPONENTS, function(x) x$id, character(1))

wildlife_data <- read_csv(WILDLIFE_FILE, show_col_types = FALSE) |>
  dplyr::rename(ba_code = code, scientific_name = binomial) |>
  dplyr::mutate(
    dplyr::across(dplyr::any_of(c(
      vapply(WILD_COMPONENTS, function(x) x$col, character(1)),
      "pest_resist_score", "client_do_not_plant", "pa_common_riparian", WETLAND_COL
    )), as.numeric),
    # Riparian: top two steps of the 5-step wetland scale (FACW/OBL-type) in the
    # chosen region, or flagged "commonly riparian" by the PA tool.
    riparian_recommended =
      (!is.na(.data[[WETLAND_COL]]) & .data[[WETLAND_COL]] >= 0.75) |
      (!is.na(pa_common_riparian) & pa_common_riparian == 1),
    # Pest flag: client says do not plant, or lowest step of the client scale.
    has_pest_flag =
      (!is.na(client_do_not_plant) & client_do_not_plant == 1) |
      (!is.na(pest_resist_score) & pest_resist_score == 0)
  ) |>
  dplyr::select(ba_code, scientific_name,
                dplyr::all_of(vapply(WILD_COMPONENTS, function(x) x$col, character(1))),
                pest_resist_score, client_do_not_plant,
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


# ── Scoring rasters: SCALED (PredRasScale) ──────────────────────────────────────────
# Used by gjam_range, get_gjam_vals, Decision-Making composite scores.
# Must remain on PredRasScale — do NOT change this path.
.gjam_base <- path.expand(
  "~/Documents/PBGJAM-data-explorer/data4Tate/gjamCreateMap/Trees/PredRasScale/rcp45"
)
gjam_files <- lapply(
  setNames(c("2040_2069","2070_2099"), c("2040_2069","2070_2099")),
  function(period) {
    dir  <- file.path(.gjam_base, period)
    tifs <- list.files(dir, pattern="\\.tif$", full.names=TRUE)
    keys <- sub(paste0("_",period,"_rcp45\\.tif$"), "", sub("^mean_","",basename(tifs)))
    setNames(tifs, keys)
  }
)

# ── Visualization rasters: RAW (PredRas) ───────────────────────────────────────────
# Used ONLY for map display in the Visualizations tab.
# Kept separate from gjam_files so scoring is completely unaffected.
.vis_gjam_base <- path.expand(
  "~/Documents/PBGJAM-data-explorer/data4Tate/gjamCreateMap/Trees/PredRas/rcp45"
)
vis_gjam_files <- lapply(
  setNames(c("2040_2069","2070_2099"), c("2040_2069","2070_2099")),
  function(period) {
    dir  <- file.path(.vis_gjam_base, period)
    tifs <- list.files(dir, pattern="\\.tif$", full.names=TRUE)
    keys <- sub(paste0("_",period,"_rcp45\\.tif$"), "", sub("^mean_","",basename(tifs)))
    setNames(tifs, keys)
  }
)

# ── Historical BA rasters — raw (visualization only) ───────────────────────────
.hist_base <- path.expand(
  "~/Documents/PBGJAM-data-explorer/data4Tate/gjamCreateMap/Trees/PredRas/hist"
)
hist_files <- local({
  tifs <- list.files(.hist_base, pattern="\\.tif$", full.names=TRUE)
  if (length(tifs) == 0) {
    message("Historical BA rasters not found at: ", .hist_base)
    list()
  } else {
    keys <- sub("_hist\\.tif$", "", sub("^mean_", "", basename(tifs)))
    keys <- sub("\\.tif$", "", keys)
    setNames(tifs, keys)
  }
})

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

gjam_range <- local({
  tifs <- unname(gjam_files[["2040_2069"]])
  tifs <- tifs[file.exists(tifs)]
  if (length(tifs) == 0) return(c(0, 1))
  sample_tifs <- tifs[seq(1, length(tifs), length.out = min(20, length(tifs)))]
  vals <- unlist(lapply(sample_tifs, function(f) {
    v <- terra::values(terra::rast(f), na.rm = TRUE)
    v[is.finite(v)]
  }))
  if (length(vals) == 0) c(0, 1) else c(min(vals), max(vals))
})

# Stumpage prices are strongly right-skewed (median ~$21/MBF, max ~$491/MBF), so a
# linear min-max squashed ~86% of species below 0.1 and the timber slider had far less
# influence than its weight implied. Normalize on a log scale instead: this preserves
# the exact ordering but spreads species across the usable range.
timber_scale <- function(p) log10(pmax(p, 0) + 1)

timber_range <- local({
  pc <- c("R1","R2","R3","R4","R5","R6","R8","R9","R10")
  ap <- suppressWarnings(as.numeric(unlist(timber[,intersect(pc,names(timber))],use.names=FALSE)))
  ap <- ap[is.finite(ap)]
  if (length(ap)==0) c(NA_real_,NA_real_) else c(min(timber_scale(ap)), max(timber_scale(ap)))
})

norm_global <- function(x, rng) {
  if (any(is.na(rng))||rng[2]==rng[1]) return(rep(NA_real_,length(x)))
  pmin(pmax((x-rng[1])/(rng[2]-rng[1]),0),1)
}

# N-component weighted mean. `vals` is n x k, `wts` length k, `active` length k.
# Inactive (switched-off) components are ignored; NA scores are skipped and the
# remaining weights are re-normalised, so unrated species are not penalised.
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

get_gjam_vals <- function(site, keys, period) {
  files <- gjam_files[[period]]
  if (is.null(files) || length(files) == 0)
    return(list(values=setNames(rep(NA_real_, length(keys)), keys), fallback_keys=character(0)))
  
  pt   <- terra::vect(matrix(c(site$lng, site$lat), ncol=2), crs="EPSG:4326")
  poly <- NULL
  if (identical(site$type, "polygon") && length(site$poly_lng) >= 3)
    poly <- tryCatch(site_polygon_vect(site), error=function(e) NULL)
  
  # Polygon: mean of the cells inside it (centroid cell if the polygon is smaller than
  # a cell). Point: the cell under the pin. The geometry is projected into each raster's
  # own CRS first so this works for lon/lat rasters and for the EPSG:5070 GEE rasters.
  extract_one <- function(fpath) {
    tryCatch({
      r  <- terra::rast(fpath)
      rc <- terra::crs(r)
      use_pt   <- if (nzchar(rc)) terra::project(pt, rc) else pt
      use_poly <- if (!is.null(poly) && nzchar(rc)) terra::project(poly, rc) else poly
      v <- NA_real_
      if (!is.null(use_poly))
        v <- suppressWarnings(as.numeric(terra::extract(r, use_poly, fun=mean, na.rm=TRUE)[1, -1]))
      if (!is.finite(v))
        v <- as.numeric(terra::extract(r, use_pt)[1, -1])
      v
    }, error = function(e) NA_real_)
  }
  
  all_keys      <- names(files)
  fallback_keys <- character(0)
  
  values <- sapply(keys, function(k) {
    if (k %in% all_keys && file.exists(files[[k]])) {
      v <- extract_one(files[[k]])
      if (is.finite(v)) return(v)
    }
    genus_lower <- tolower(regmatches(k, regexpr("^[a-z]+", k)))
    genus_mask  <- startsWith(tolower(all_keys), genus_lower) &
      file.exists(unlist(files[all_keys]))
    genus_keys  <- all_keys[genus_mask]
    gv <- vapply(genus_keys, function(gk) extract_one(files[[gk]]), numeric(1))
    gv <- gv[is.finite(gv)]
    if (length(gv) > 0) fallback_keys <<- c(fallback_keys, k)
    if (length(gv) == 0) NA_real_ else mean(gv)
  })
  list(values=values, fallback_keys=fallback_keys)
}

safe_vals <- function(r) { v <- raster::values(r); v[!is.na(v)&is.finite(v)] }

# ── Raster cache ──────────────────────────────────────────────────────────────
.raster_cache    <- new.env(hash=TRUE, parent=emptyenv())
.raster_cache_dir <- local({
  d <- file.path(tempdir(), "pbgjam_raster_cache")
  dir.create(d, showWarnings=FALSE, recursive=TRUE)
  d
})

.prep_raster <- function(terra_rast) {
  raster::raster(terra_rast)
}

.cached_load <- function(cache_key, load_fn) {
  if (exists(cache_key, envir=.raster_cache))
    return(get(cache_key, envir=.raster_cache))
  rds_path <- file.path(.raster_cache_dir, paste0(cache_key, ".rds"))
  if (file.exists(rds_path)) {
    result <- readRDS(rds_path)
    assign(cache_key, result, envir=.raster_cache)
    return(result)
  }
  result <- load_fn()
  tryCatch(saveRDS(result, rds_path), error=function(e) NULL)
  assign(cache_key, result, envir=.raster_cache)
  result
}

load_gjam_raster_for_key <- function(key, period) {
  cache_key <- paste0("v2_future_", period, "_", key)
  .cached_load(cache_key, function() {
    files <- vis_gjam_files[[period]]
    tryCatch({
      if (key %in% names(files) && file.exists(files[[key]])) {
        list(raster=.prep_raster(terra::rast(files[[key]])), fallback=FALSE)
      } else {
        genus_lower <- tolower(regmatches(key, regexpr("^[a-z]+", key)))
        genus_keys  <- names(files)[startsWith(tolower(names(files)), genus_lower)]
        genus_keys  <- genus_keys[vapply(genus_keys, function(gk) file.exists(files[[gk]]), logical(1))]
        if (length(genus_keys) == 0) return(list(raster=NULL, fallback=TRUE))
        stack_r <- terra::rast(lapply(genus_keys, function(gk) terra::rast(files[[gk]])))
        r_mean  <- terra::app(stack_r, fun="mean", na.rm=TRUE)
        list(raster=.prep_raster(r_mean), fallback=TRUE)
      }
    }, error=function(e) list(raster=NULL, fallback=TRUE))
  })
}

load_hist_raster_for_key <- function(key) {
  if (length(hist_files) == 0) return(list(raster=NULL, fallback=FALSE))
  cache_key <- paste0("v2_hist_", key)
  .cached_load(cache_key, function() {
    tryCatch({
      if (key %in% names(hist_files) && file.exists(hist_files[[key]])) {
        list(raster=.prep_raster(terra::rast(hist_files[[key]])), fallback=FALSE)
      } else {
        genus_lower <- tolower(regmatches(key, regexpr("^[a-z]+", key)))
        genus_keys  <- names(hist_files)[startsWith(tolower(names(hist_files)), genus_lower)]
        genus_keys  <- genus_keys[vapply(genus_keys, function(gk) file.exists(hist_files[[gk]]), logical(1))]
        if (length(genus_keys) == 0) return(list(raster=NULL, fallback=TRUE))
        stack_r <- terra::rast(lapply(genus_keys, function(gk) terra::rast(hist_files[[gk]])))
        r_mean  <- terra::app(stack_r, fun="mean", na.rm=TRUE)
        list(raster=.prep_raster(r_mean), fallback=TRUE)
      }
    }, error=function(e) list(raster=NULL, fallback=TRUE))
  })
}

all_vis_keys <- unique(c(
  names(vis_gjam_files[["2040_2069"]]),
  names(vis_gjam_files[["2070_2099"]]),
  names(hist_files)
))
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
  --wild-c:   #B45309;
  --timb-c:   #0E7490;
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

/* Wildlife umbrella group */
.wildlife-group {
  background: rgba(0,0,0,0.2);
  border: 1px solid rgba(255,255,255,0.07);
  border-left: 3px solid var(--wild-c);
  border-radius: var(--radius);
  padding: 12px 14px 4px; margin-bottom: 10px;
  box-shadow: inset 0 1px 3px rgba(0,0,0,0.15);
}
.wildlife-head {
  font-family: var(--font); font-size: 12px; font-weight: 600; color: white;
}
.wildlife-sub {
  font-family: var(--font); font-size: 9px; color: var(--sage);
  display: block; margin: 2px 0 10px; line-height: 1.4;
}
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
.dm-pest-note.warn { border-color: rgba(127,29,29,0.5); background: rgba(127,29,29,0.06); }
.dm-pest-note .src { color: var(--ink2); font-style: italic; }
.about-data { padding-top: 0; }
.about-data h3 { font-family: var(--font); font-size: 15px; font-weight: 600; margin: 22px 0 6px; color: var(--ink); }
.about-data p, .about-data li { font-family: var(--font); font-size: 12px; line-height: 1.6; color: var(--ink); }

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
.flag-tag.flag-pest { color: var(--danger); background: rgba(127,29,29,0.12); }

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
  font-family: var(--font); font-size: 8px; font-weight: 600;
  background: #FEF3C7; color: #92400E;
  border: 1px solid #F59E0B;
  border-radius: 8px; padding: 1px 6px; margin-top: 4px;
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
  if (type == "wildlife") return(sprintf("rgba(180,83,9,%.2f)",    0.35 + norm_val*0.65))
  if (type == "timber")   return(sprintf("rgba(14,116,144,%.2f)",  0.35 + norm_val*0.65))
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
    
// ── Shared time-range toggle ──────────────────────────────────────────────────
$(document).on('click', '.time-seg', function() {
  var $toggle = $(this).closest('.time-toggle');
  $toggle.find('.time-seg').removeClass('active');
  $(this).addClass('active');
  var inputId = $toggle.data('inputid') || 'time_range';
  Shiny.setInputValue(inputId, $(this).data('val'), {priority:'event'});
});

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

// ── Layer single-select ───────────────────────────────────────────────────────
$(document).on('click', '.vis-layer-btn', function() {
  $('.vis-layer-btn').removeClass('active');
  $(this).addClass('active');
  Shiny.setInputValue('vis_active_layer', $(this).data('layer'), {priority:'event'});
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
  // Set all time toggles to 2040-2069 by default
  Shiny.setInputValue('time_range_loc', '2040_2069', {priority:'event'});
  Shiny.setInputValue('time_range_dm',  '2040_2069', {priority:'event'});
  Shiny.setInputValue('time_range_vis', '2040_2069', {priority:'event'});
  Shiny.setInputValue('vis_active_layer', 'future_ba', {priority:'event'});
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
                               div(class="weight-group",
                                   div(class="weight-label", "Predicted Future BA Score"),
                                   sliderInput("w1", NULL, min=0, max=1, value=0.34, step=0.01, ticks=FALSE, width="100%"),
                                   div(class="weight-note", "Predicted basal area from model across US by given timeframe")
                               ),
                               div(class="wildlife-group",
                                   div(class="wildlife-head", "Wildlife"),
                                   tags$span(class="wildlife-sub",
                                             "Switch each part on or off and set its weight. Species with no rating for a part are not penalised."),
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
                               div(class="weight-group",
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
                               div(class="time-toggle", `data-inputid`="time_range_vis",
                                   div(class="time-seg active", `data-val`="2040_2069", "2040–2069"),
                                   div(class="time-seg",        `data-val`="2070_2099", "2070–2099")
                               ),
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
                                   div(class="time-toggle-label", "Prediction period"),
                                   div(class="time-toggle", `data-inputid`="time_range_loc",
                                       div(class="time-seg active", `data-val`="2040_2069", "2040–2069"),
                                       div(class="time-seg",        `data-val`="2070_2099", "2070–2099")
                                   )
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
                                   div(class="time-toggle-wrap",
                                       div(class="time-toggle-label", "Period"),
                                       # 2040-2069 marked active by default
                                       div(class="time-toggle", `data-inputid`="time_range_loc",
                                           div(class="time-seg active", `data-val`="2040_2069", "2040–2069"),
                                           div(class="time-seg",        `data-val`="2070_2099", "2070–2099")
                                       )
                                   )
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
                                               div(class="vis-layer-btn-label", "Predicted Future BA")
                                           ),
                                           div(class="vis-layer-btn", `data-layer`="hist_ba",
                                               div(class="vis-layer-btn-label", "Historical BA")
                                           ),
                                           div(class="vis-layer-btn", `data-layer`="fecundity",
                                               div(class="vis-layer-btn-label", "Fecundity")
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
                                   "PBGJAM v2 was developed by Tate Commission, Dr. Tong Qiu, and ",
                                   "Dr. James Clark from Duke University."
                               )
                           ),
                           div(class="about-inner about-data",
                               tags$h3("Data and methods"),
                               tags$p("Species are ranked at a site by a weighted average of the components you switch on. ",
                                      "Each component is scaled 0 to 1. Species with no rating for a component are not penalised: ",
                                      "the weights are re-normalised over the components that do have a value."),
                               tags$h3("Predicted basal area"),
                               tags$p("Future basal area per species comes from a generalized joint attribute model (GJAM) fit to USFS Forest Inventory and Analysis plots and climate, topography and soil predictors."),
                               tags$h3("Timber"),
                               tags$p("Regional stumpage price from USFS timber sale reports, 5-year average (FY2021-FY2025), on a log scale."),
                               tags$h3("Wildlife"),
                               tags$p("Wildlife is a group of separate, optional components rather than one combined score:"),
                               tags$ul(
                                 tags$li("Overall wildlife value: the wildlife-value rating in the Pennsylvania DCNR Bureau of Forestry Native Tree Selector."),
                                 tags$li("Fruit, nut and seed food: USDA PLANTS fruit/seed abundance and product traits."),
                                 tags$li("Pollinators: bees, butterflies, hummingbirds and other insects listed for the species in the PA Native Tree Selector. A blank means unknown, not none."),
                                 tags$li("Deer forage and deer resistance: PA Native Tree Selector deer palatability, with USDA PLANTS browse palatability as a fallback. These are opposites, so they are off by default; use one, not both.")
                               ),
                               tags$p("These sources do not contain bird, mammal, reptile, amphibian or caterpillar host-plant scores, so those are not part of the ranking."),
                               tags$h3("Pest and pathogen risk"),
                               tags$p("Pennsylvania pest and pathogen ratings supplied by the project client take priority; the PA Native Tree Selector pest notes are used where the client has no rating. ",
                                      "The Decision-Making tab shows the text and its source. Species with no information are shown as not rated."),
                               tags$h3("Riparian status"),
                               tags$p("A species counts as riparian if its US Army Corps National Wetland Plant List indicator is in the wettest two classes for the region, or the PA Native Tree Selector lists it as commonly riparian."),
                               tags$h3("Scoring rule"),
                               tags$p("Ordinal ratings are converted to equal-step linear scales (for example 0 / 0.5 / 1, or 0 / 0.25 / 0.5 / 0.75 / 1) based on each source's number of levels, without hand-chosen weights. Names are matched to the GJAM species list through USDA PLANTS."),
                               tags$h3("Limits"),
                               tags$p("The PA Native Tree Selector is a compilation; its authors note that sources conflicted and they used best judgement. Coverage is greatest for Pennsylvania natives and thinner for other species. Treat scores as screening aids, not site prescriptions.")
                           ),
                           div(class="about-inner",
                               div(class="about-person",
                                   tags$img(class="about-photo", src="tatecommission.png", alt="Tate Commission"),
                                   div(class="about-person-info",
                                       div(class="about-name", "Tate Commission"),
                                       div(class="about-role", "Lead Developer"),
                                       tags$p(class="about-bio",
                                              "Tate Commission is a sophomore at Duke University studying Statistical Science. His interests include creating beautiful data visualizations, developing statistical simulations to understand complex systems, and building deep learning for satellite imagery analysis. He has extensive experience developing interactive tools for large-scale geospatial analysis through R and R Shiny. His work has supported projects with American Forests, NASA, and the African Parks Foundation including through the Qiu Lab at Duke University.")
                                   )
                               ),
                               div(class="about-person",
                                   tags$img(class="about-photo", src="tongqiu.png", alt="Dr. Tong Qiu"),
                                   div(class="about-person-info",
                                       div(class="about-name", "Dr. Tong Qiu"),
                                       div(class="about-role", "Faculty Advisor"),
                                       tags$p(class="about-bio",
                                              "Dr. Tong Qiu is an Assistant Professor of Ecology at Duke University. He is interested in understanding the causes and consequences of biodiversity change at scales ranging from individual organisms to the entire biosphere. As the Principal Investigator of the Qiu lab, Qiu develops data-model synthesis frameworks that integrate remote sensing (e.g., LiDAR, hyperspectral imaging), field sampling, and ecological monitoring networks with Bayesian hierarchical models and Earth System models.")
                                   )
                               ),
                               div(class="about-person",
                                   tags$img(class="about-photo", src="jamesclark.png", alt="Dr. James Clark"),
                                   div(class="about-person-info",
                                       div(class="about-name", "Dr. James Clark"),
                                       div(class="about-role", "Faculty Advisor"),
                                       tags$p(class="about-bio",
                                              "Dr. James Clark is the Nicholas Distinguished Professor of Environmental Science at Duke University and a Fellow of the Ecological Society of America. Clark's lab uses long-term experiments and monitoring studies to understand disturbance and climate controls on ecosystem dynamics. He has authored more than 250 refereed scientific articles and published four books, and he pioneered the framework that serves as the backbone of this application.")
                                   )
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
    time_range         = "2040_2069",
    # Visualizations tab — only the active raster is held in memory
    vis_raster         = NULL,
    vis_raster_range   = c(NA_real_, NA_real_),
    vis_fallback       = FALSE,
    usfs_region        = NULL
  )
  
  observeEvent(input$time_range_loc, { rv$time_range <- input$time_range_loc })
  observeEvent(input$time_range_dm,  { rv$time_range <- input$time_range_dm  })
  observeEvent(input$time_range_vis, { rv$time_range <- input$time_range_vis })
  
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
  
  # Which wildlife parts are switched on, and their weights (defaults until inputs exist).
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
    
    gjam_result        <- get_gjam_vals(site, keys, rv$time_range)
    gjam_raw           <- gjam_result$values
    gjam_n             <- norm_global(gjam_raw, gjam_range)
    gjam_fallback_keys <- gjam_result$fallback_keys
    
    ws <- wild_state()
    idx <- match(keys, wildlife_data$ba_code)
    wild_mat <- vapply(ws$col, function(cn) as.numeric(wildlife_data[[cn]][idx]),
                       numeric(length(keys)))
    if (!is.matrix(wild_mat)) wild_mat <- matrix(wild_mat, nrow=length(keys))
    colnames(wild_mat) <- paste0("wp_", ws$id)
    
    # Timber price is regional: for a polygon, the USFS region of its centroid is used.
    timber_info    <- lapply(keys, function(k) get_timber_price_info(site$lng, site$lat, k))
    timber_prices  <- vapply(timber_info, function(x) x$price,  numeric(1))
    timber_sources <- vapply(timber_info, function(x) x$source, character(1))
    timber_n       <- norm_global(timber_scale(timber_prices), timber_range)
    
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
      gjam_is_genus_avg = keys %in% gjam_fallback_keys,
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
      addProviderTiles("CartoDB.Positron",
                       options=providerTileOptions(maxZoom=18)) %>%
      addPolygons(data=state_sf, fillOpacity=0, color="#2D5A48",
                  weight=0.8, opacity=0.65) %>%
      fit_pa(padding=c(20, 20)) %>%
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
    if (is.null(key) || !key %in% names(species_display))
      return(div(class="vis-species-title",
                 style="font-style:normal;color:var(--ink2);font-size:11px;",
                 "Select a species to view its basal area maps"))
    fb_note <- if (isTRUE(rv$vis_fallback))
      div(class="vis-fallback-note", "\u00f8 genus average")
    else NULL
    tagList(div(class="vis-species-title", species_display[key]), fb_note)
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
  
  # BA and fecundity layers: need a species key
  observeEvent(list(input$vis_selected_species, rv$time_range, input$vis_active_layer), {
    key    <- input$vis_selected_species
    period <- if (!is.null(rv$time_range)) rv$time_range else "2040_2069"
    layer  <- if (!is.null(input$vis_active_layer)) input$vis_active_layer else "future_ba"
    
    if (layer == "canopy") return()
    if (is.null(key) || key == "") return()
    
    if (layer == "future_ba") {
      result <- load_gjam_raster_for_key(key, period)
      rv$vis_raster       <- result$raster
      rv$vis_fallback     <- isTRUE(result$fallback)
      rv$vis_raster_range <- if (!is.null(result$raster)) {
        v <- safe_vals(result$raster); if (length(v)) range(v) else c(NA,NA)
      } else c(NA, NA)
    } else if (layer == "hist_ba") {
      result <- load_hist_raster_for_key(key)
      rv$vis_raster       <- result$raster
      rv$vis_fallback     <- isTRUE(result$fallback)
      rv$vis_raster_range <- if (!is.null(result$raster)) {
        v <- safe_vals(result$raster); if (length(v)) range(v) else c(NA,NA)
      } else c(NA, NA)
    } else if (layer == "fecundity") {
      # Fecundity uses same source data as future_ba (dummy data)
      # Replace with real fecundity rasters when available
      result <- load_gjam_raster_for_key(key, period)
      rv$vis_raster       <- result$raster
      rv$vis_fallback     <- isTRUE(result$fallback)
      rv$vis_raster_range <- if (!is.null(result$raster)) {
        v <- safe_vals(result$raster); if (length(v)) range(v) else c(NA,NA)
      } else c(NA, NA)
    }
  }, ignoreNULL=FALSE)
  
  # Orange palette for BA layers
  orange_pal <- function(rng) {
    colorNumeric(
      palette  = c("#FFF8EC","#FDDBA5","#F59E0B","#B45309","#7C2D12"),
      domain   = rng,
      na.color = "transparent"
    )
  }
  
  # Blue palette for fecundity layer
  blue_pal <- function(rng) {
    colorNumeric(
      palette  = c("#EFF6FF","#BFDBFE","#60A5FA","#2563EB","#1E3A8A"),
      domain   = rng,
      na.color = "transparent"
    )
  }
  
  VIS_RASTER_ID <- "vis_ba_raster"
  VIS_TILE_ID   <- "vis_canopy_tile"
  VIS_MAX_BYTES <- 32L * 1024L * 1024L
  
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
    
    if (layer %in% c("future_ba", "hist_ba")) {
      r   <- isolate(rv$vis_raster)
      rng <- isolate(rv$vis_raster_range)
      if (!is.null(r) && !any(is.na(rng))) {
        pal <- orange_pal(rng)
        proxy <- proxy %>% addRasterImage(
          r, colors=pal, opacity=0.78, project=TRUE,
          layerId=VIS_RASTER_ID, maxBytes=VIS_MAX_BYTES
        )
      }
    } else if (layer == "fecundity") {
      # render fecundity raster in blue
      r   <- isolate(rv$vis_raster)
      rng <- isolate(rv$vis_raster_range)
      if (!is.null(r) && !any(is.na(rng))) {
        pal <- blue_pal(rng)
        proxy <- proxy %>% addRasterImage(
          r, colors=pal, opacity=0.78, project=TRUE,
          layerId=VIS_RASTER_ID, maxBytes=VIS_MAX_BYTES
        )
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
    
    # fecundity uses blue gradient in legend
    if (layer == "fecundity") {
      ticks <- pretty(rng, n=4)
      ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
      if (length(ticks) < 2) ticks <- round(c(rng[1], mean(rng), rng[2]), 2)
      return(tagList(
        tags$span(class="legend-title", "Fecundity (m\u00b2/ha) \u2014 dummy data"),
        div(class="legend-gradient",
            style="background:linear-gradient(to right,#EFF6FF,#BFDBFE,#60A5FA,#2563EB,#1E3A8A);"),
        div(class="legend-ticks", lapply(formatC(ticks, format="f", digits=2), tags$span))
      ))
    }
    
    title <- if (layer == "future_ba") "Predicted Future BA (m\u00b2/ha)"
    else "Historical BA (m\u00b2/ha)"
    
    ticks <- pretty(rng, n=4)
    ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
    if (length(ticks) < 2) ticks <- round(c(rng[1], mean(rng), rng[2]), 2)
    tick_labels <- formatC(ticks, format="f", digits=2)
    
    tagList(
      tags$span(class="legend-title", title),
      div(class="legend-gradient",
          style="background:linear-gradient(to right,#FFF8EC,#FDDBA5,#F59E0B,#B45309,#7C2D12);"),
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
      pest_flag        <- nrow(wd) > 0 && isTRUE(wd$has_pest_flag[1])
      riparian_flag    <- nrow(wd) > 0 && isTRUE(wd$riparian_recommended[1])
      riparian_mode_on <- !is.null(input$riparian_mode) && input$riparian_mode > 0
      div(class=cls, `data-key`=row$key,
          div(class="dm-sp-rank", i),
          div(class="dm-sp-name", row$display_name,
              if (riparian_mode_on && riparian_flag)
                tags$span(class="flag-tag flag-riparian", "RP") else NULL,
              if (pest_flag) tags$span(class="flag-tag flag-pest", "PEST") else NULL),
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
    fmt <- function(x) if (is.na(x)) "\u2014" else sprintf("%.3f", x)
    ws  <- wild_state()
    ws_on <- ws[ws$on, , drop = FALSE]
    
    if (is.null(df) || is.null(key)) {
      return(div(class="dm-score-layout",
                 div(class="dm-composite-card",
                     div(class="dm-composite-label", "Composite Score"),
                     div(class="dm-composite-val", "\u2014")),
                 div(class="dm-component-row",
                     div(class="dm-comp-card gjam-card",
                         tags$span(class="dm-comp-label","Predicted Future BA"),
                         div(class="dm-comp-val","\u2014")),
                     div(class="dm-comp-card wild-card",
                         tags$span(class="dm-comp-label","Wildlife"),
                         div(class="dm-comp-val","\u2014")),
                     div(class="dm-comp-card timb-card",
                         tags$span(class="dm-comp-label","Timber Value"),
                         div(class="dm-comp-val","\u2014"))
                 )
      ))
    }
    
    row <- df[df$key == key, ]
    if (nrow(row) == 0) return(NULL)
    
    n_total <- 2 + nrow(ws_on)
    mode <- if (is.null(input$riparian_mode)) 0 else input$riparian_mode
    comp_sub <- if (mode == 2) {
      "riparian-required mode"
    } else if (!is.na(row$n_components[1]) && row$n_components[1] < n_total) {
      paste0("based on ", row$n_components[1], " of ", n_total, " components")
    } else paste0("weighted composite of ", n_total, " components")
    
    gjam_raw_note <- if (isTRUE(row$gjam_is_genus_avg[1])) {
      div(class="genus-avg-badge", "\u00f8 genus-level average")
    } else div(class="dm-comp-raw", "direct model prediction")
    
    wild_card_body <- if (nrow(ws_on) == 0) {
      div(class="dm-comp-raw", "All wildlife parts are switched off")
    } else {
      tagList(lapply(seq_len(nrow(ws_on)), function(j) {
        v <- row[[paste0("wp_", ws_on$id[j])]][1]
        div(class="dm-wild-part",
            span(ws_on$label[j]),
            span(if (is.na(v)) "not rated" else sprintf("%.2f", v)))
      }))
    }
    
    wd <- wildlife_data[wildlife_data$ba_code == key, ]
    pest_txt <- if (nrow(wd) > 0) wd$pest_text_display[1] else NA_character_
    pest_src <- if (nrow(wd) > 0) wd$pest_text_source[1] else NA_character_
    pest_warn <- nrow(wd) > 0 && isTRUE(wd$has_pest_flag[1])
    pest_score <- if (nrow(wd) > 0) wd$pest_resist_score[1] else NA_real_
    pest_block <- div(class=paste("dm-pest-note", if (pest_warn) "warn" else ""),
        tags$b("Pest and pathogen risk"),
        if (!is.na(pest_score))
          div(sprintf("Client resistance rating: %.2f (0 = lowest, 1 = highest)", pest_score)),
        if (!is.na(pest_txt) && nzchar(pest_txt)) div(pest_txt) else div("Not rated: no pest information in the available sources."),
        if (!is.na(pest_src) && nzchar(pest_src)) div(class="src", paste0("Source: ", pest_src)) else NULL)
    
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
        div(class="bar-track", style="height:5px;margin:1px 0;", title=ws_on$label[j],
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

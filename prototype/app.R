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

wildlife_data <- read_csv(
  "BA_wildlife_traits.csv",
  show_col_types = FALSE
) |>
  dplyr::select(
    ba_code, scientific_name,
    wildlife_value_index,
    cat_A, cat_B, cat_C, cat_D,
    lep_spp_supported,
    berry_nut_seed_product, fruit_seed_abundance, fruit_seed_persistence,
    palatable_browse_animal, bloom_period,
    has_pest_flag, max_pest_severity, pest_flag_list, pest_flag_notes,
    wetland_wettest, wetland_modal, riparian_category, riparian_recommended
  ) |>
  mutate(
    wildlife_value_index = as.numeric(wildlife_value_index),
    cat_A                = as.numeric(cat_A),
    cat_B                = as.numeric(cat_B),
    cat_C                = as.numeric(cat_C),
    cat_D                = as.numeric(cat_D),
    has_pest_flag_raw    = has_pest_flag,
    has_pest_flag        = case_when(
      is.na(has_pest_flag_raw) ~ FALSE,
      toupper(trimws(as.character(has_pest_flag_raw))) %in%
        c("TRUE","T","YES","Y","1") ~ TRUE,
      TRUE ~ as.logical(has_pest_flag_raw)
    ),
    riparian_recommended = case_when(
      is.na(riparian_recommended) ~ FALSE,
      toupper(trimws(as.character(riparian_recommended))) %in%
        c("TRUE","T","YES","Y","1") ~ TRUE,
      TRUE ~ as.logical(riparian_recommended)
    )
  ) |>
  dplyr::select(-has_pest_flag_raw)

wildlife_range <- range(wildlife_data$wildlife_value_index, na.rm = TRUE)

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
#   1. Go to https://code.earthengine.google.com with qt email
#   2. Paste and run initial_chm_script
#   3. Copy the printed URL (decode %7B/%7D → {z}/{x}/{y}) and paste below.
#
canopy_tile_url <- "https://earthengine.googleapis.com/v1/projects/pbgjam-v2/maps/77ab79e36978d8bde0608e8370a50a02-ad7d4b9e68fa58d78c57ba645be8be48/tiles/%7Bz%7D/%7Bx%7D/%7By%7D"

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat/2)^2 + cos(lat1*to_rad) * cos(lat2*to_rad) * sin(dlon/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

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

get_timber_price <- function(lng, lat, species_key) {
  row <- timber[timber$gjam_species == species_key, ]
  if (nrow(row) == 0) return(NA_real_)
  region <- get_usfs_region(lng, lat)
  price_for_region <- function(reg) {
    if (!reg %in% names(row)) return(NA_real_)
    v <- suppressWarnings(as.numeric(row[[reg]][1]))
    if (is.na(v) || !is.finite(v)) NA_real_ else v
  }
  p <- price_for_region(region)
  if (!is.na(p)) return(p)
  d_km   <- haversine_km(lng, lat, .usfs_centroids$lon, .usfs_centroids$lat)
  ranked <- .usfs_centroids$region[order(d_km)]
  ranked <- ranked[ranked != region]
  for (reg in ranked) {
    pt <- price_for_region(reg)
    if (!is.na(pt)) return(pt)
  }
  NA_real_
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

timber_range <- local({
  pc <- c("R1","R2","R3","R4","R5","R6","R8","R9","R10")
  ap <- suppressWarnings(as.numeric(unlist(timber[,intersect(pc,names(timber))],use.names=FALSE)))
  ap <- ap[is.finite(ap)]
  if (length(ap)==0) c(NA_real_,NA_real_) else c(min(ap),max(ap))
})

norm_global <- function(x, rng) {
  if (any(is.na(rng))||rng[2]==rng[1]) return(rep(NA_real_,length(x)))
  pmin(pmax((x-rng[1])/(rng[2]-rng[1]),0),1)
}

weighted_composite_full <- function(gjam_n, wild_n, timber_n, w1, w2, w3) {
  n <- length(gjam_n)
  score <- numeric(n)
  n_components <- integer(n)
  for (i in seq_len(n)) {
    vals <- c(gjam_n[i], wild_n[i], timber_n[i])
    wts  <- c(w1, w2, w3)
    ok   <- !is.na(vals)
    n_components[i] <- sum(ok)
    wt_sum <- sum(wts[ok])
    score[i] <- if (!any(ok) || wt_sum == 0) NA_real_ else sum(vals[ok]*wts[ok]) / wt_sum
  }
  list(score = score, n_components = n_components)
}

options(tigris_use_cache=TRUE)
state_sf <- tigris::states(cb=TRUE, resolution="500k", progress_bar=FALSE)
state_sf <- state_sf[!state_sf$STUSPS %in% c("AK","HI","PR","VI","GU","MP","AS","DC"),]
state_sf <- sf::st_transform(state_sf, 4326)

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

get_all_species <- function(lng, lat, k = 20) {
  d_km    <- haversine_km(lng, lat, coords_lon, coords_lat)
  idx     <- order(d_km)[seq_len(min(k, nrow(BA_mat)))]
  sub     <- BA_mat[idx, , drop = FALSE]
  if (!is.matrix(sub) || typeof(sub) == "list") {
    sub <- matrix(as.numeric(unlist(sub)), nrow = length(idx),
                  dimnames = list(NULL, colnames(BA_mat)))
  }
  present <- colnames(sub)[colSums(sub, na.rm = TRUE) > 0]
  data.frame(key          = present,
             display_name = as.character(species_display[present]),
             stringsAsFactors = FALSE)
}

get_gjam_vals <- function(lng, lat, keys, period) {
  files <- gjam_files[[period]]
  if (is.null(files) || length(files) == 0)
    return(list(values=setNames(rep(NA_real_, length(keys)), keys), fallback_keys=character(0)))
  
  pt <- terra::vect(matrix(c(lng, lat), ncol=2), crs="EPSG:4326")
  
  extract_one <- function(fpath) {
    tryCatch({
      r <- terra::rast(fpath)
      as.numeric(terra::extract(r, pt)[1, -1])
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
# In-memory cache keyed by "future_2040_2069_quealba" etc.
# On first load: read raw TIF, reproject to WGS84, save as .rds to a session
# temp dir so the projected version survives server restarts mid-session.
# On repeat access: load from .rds (fast) or in-memory (instant).
.raster_cache    <- new.env(hash=TRUE, parent=emptyenv())
.raster_cache_dir <- local({
  d <- file.path(tempdir(), "pbgjam_raster_cache")
  dir.create(d, showWarnings=FALSE, recursive=TRUE)
  d
})

# Keep raster in native CRS; addRasterImage(project=TRUE) reprojects correctly.
# Pre-projecting via terra then wrapping in raster::raster() mis-registers the
# extent, shifting the layer north and into Canada. Let leaflet handle it.
.prep_raster <- function(terra_rast) {
  raster::raster(terra_rast)
}

.cached_load <- function(cache_key, load_fn) {
  # 1. In-memory hit
  if (exists(cache_key, envir=.raster_cache))
    return(get(cache_key, envir=.raster_cache))
  # 2. Disk hit (projected .rds from a previous load this session)
  rds_path <- file.path(.raster_cache_dir, paste0(cache_key, ".rds"))
  if (file.exists(rds_path)) {
    result <- readRDS(rds_path)
    assign(cache_key, result, envir=.raster_cache)
    return(result)
  }
  # 3. Cold load: read TIF, project, cache
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
  padding: 12px 12px 8px;
  border-bottom: 1px solid var(--border);
  background: var(--linen); flex-shrink: 0;
}
.vis-sidebar-header h4 {
  font-family: var(--font);
  font-size: 9px; font-weight: 600; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--mid); margin: 0 0 6px 0;
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
  $('.time-seg').removeClass('active'); $(this).addClass('active');
  Shiny.setInputValue('time_range', $(this).data('val'), {priority:'event'});
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

// ── Layer single-select ───────────────────────────────────────────────────────
$(document).on('click', '.vis-layer-btn', function() {
  $('.vis-layer-btn').removeClass('active');
  $(this).addClass('active');
  Shiny.setInputValue('vis_active_layer', $(this).data('layer'), {priority:'event'});
});

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
  Shiny.setInputValue('time_range', '2040_2069', {priority:'event'});
  Shiny.setInputValue('vis_active_layer', 'future_ba', {priority:'event'});
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
                               div(class="weight-group",
                                   div(class="weight-label", "Wildlife Value"),
                                   sliderInput("w2", NULL, min=0, max=1, value=0.33, step=0.01, ticks=FALSE, width="100%"),
                                   div(class="weight-note", "Combined scores from fruit, insects, pollinators, and deer")
                               ),
                               div(class="weight-group",
                                   div(class="weight-label", "Timber Market Value"),
                                   sliderInput("w3", NULL, min=0, max=1, value=0.33, step=0.01, ticks=FALSE, width="100%"),
                                   div(class="weight-note", "Regional price based on official 2025 timber sales with USFS")
                               ),
                               div(class="riparian-group",
                                   div(class="weight-label", "Riparian Status"),
                                   sliderInput("riparian_mode", NULL, min=0, max=2, value=0, step=1,
                                               ticks=TRUE, width="100%", animate=FALSE),
                                   uiOutput("riparian_mode_note")
                               ),
                               tags$hr(style="border-color:rgba(255,255,255,0.1); margin: 14px 0;"),
                               tags$span(class="sidebar-section-label", "Prediction Period"),
                               div(class="time-toggle",
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
                                   div(class="time-toggle",
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
              
              # ── Visualizations tab (formerly GJAM) ────────────────────────────────────
              tabPanel("Visualizations",
                       div(class="vis-outer",
                           # Species list sidebar
                           div(class="vis-sidebar",
                               div(class="vis-sidebar-header",
                                   tags$h4("Species"),
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
                                       div(class="time-toggle",
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
    click_lng          = NULL,
    click_lat          = NULL,
    time_range         = "2040_2069",
    # Visualizations tab — only the active raster is held in memory
    vis_raster         = NULL,
    vis_raster_range   = c(NA_real_, NA_real_),
    vis_fallback       = FALSE,
    usfs_region        = NULL
  )
  
  observeEvent(input$time_range, { rv$time_range <- input$time_range })
  
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
  
  composite_scores <- reactive({
    if (is.null(rv$click_lng)) return(NULL)
    all_sp <- tryCatch(get_all_species(rv$click_lng, rv$click_lat), error=function(e) NULL)
    if (is.null(all_sp) || nrow(all_sp) == 0) return(NULL)
    keys <- all_sp$key
    
    gjam_result        <- get_gjam_vals(rv$click_lng, rv$click_lat, keys, rv$time_range)
    gjam_raw           <- gjam_result$values
    gjam_n             <- norm_global(gjam_raw, gjam_range)
    gjam_fallback_keys <- gjam_result$fallback_keys
    
    wild_raw  <- wildlife_data$wildlife_value_index[match(keys, wildlife_data$ba_code)]
    wild_n    <- norm_global(wild_raw, wildlife_range)
    
    timber_prices <- sapply(keys, function(k) get_timber_price(rv$click_lng, rv$click_lat, k))
    timber_n      <- norm_global(timber_prices, timber_range)
    
    composite <- weighted_composite_full(gjam_n, wild_n, timber_n, input$w1, input$w2, input$w3)
    
    df <- data.frame(
      key               = keys,
      display_name      = all_sp$display_name,
      gjam              = round(gjam_n, 4),
      wildlife          = round(wild_n, 4),
      timber            = round(timber_n, 4),
      score             = round(composite$score, 4),
      n_components      = composite$n_components,
      timber_raw        = round(timber_prices, 2),
      wildlife_raw      = round(wild_raw, 1),
      gjam_is_genus_avg = keys %in% gjam_fallback_keys,
      stringsAsFactors  = FALSE
    )
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
      setView(lng=-96, lat=38, zoom=4) %>%
      htmlwidgets::onRender("function(el,x){
        L.control.zoom({position:'bottomright'}).addTo(this);
        this.on('click', function(e) {
          var lat=e.latlng.lat, lng=e.latlng.lng;
          if (lat<24||lat>50||lng<-125||lng>-66) {
            L.popup().setLatLng(e.latlng)
              .setContent('<span style=\\'font-family:Lexend,sans-serif;font-size:12px;color:#555;\\'>Select a location within the contiguous United States.</span>')
              .openOn(this);
          }
        });
      }")
  })
  
  observeEvent(input$main_map_click, {
    cl  <- input$main_map_click
    lng <- ((cl$lng+180)%%360)-180; lat <- cl$lat
    if (lat < 24 || lat > 50 || lng < -125 || lng > -66) return()
    rv$click_lng   <- lng
    rv$click_lat   <- lat
    rv$usfs_region <- get_usfs_region(lng, lat)
    
    leafletProxy("main_map") %>% clearMarkers() %>%
      addMarkers(lng=lng, lat=lat, icon=thumbtack_icon)
    leafletProxy("vis_map") %>% clearGroup("vis_pin") %>%
      addMarkers(lng=lng, lat=lat, icon=thumbtack_icon, group="vis_pin")
    
    leafletProxy("sat_minimap") %>%
      setView(lng=lng, lat=lat, zoom=POPUP_ZOOM) %>% clearMarkers() %>%
      addCircleMarkers(lng=lng, lat=lat, radius=4,
                       color="#FF5050", fillColor="#FF5050", fillOpacity=1, opacity=1,
                       weight=1.5, stroke=TRUE)
    leafletProxy("terrain_minimap") %>%
      setView(lng=lng, lat=lat, zoom=POPUP_ZOOM) %>% clearMarkers() %>%
      addCircleMarkers(lng=lng, lat=lat, radius=4,
                       color="#FF5050", fillColor="#FF5050", fillOpacity=1, opacity=1,
                       weight=1.5, stroke=TRUE)
  })
  
  output$coord_display <- renderUI({
    if (is.null(rv$click_lat))
      div(class="coord-display", "Click anywhere on the map to select a site")
    else
      div(class="coord-display", sprintf("%.4f\u00b0N  %.4f\u00b0W", rv$click_lat, abs(rv$click_lng)))
  })
  
  output$map_sidebar_info <- renderUI({
    if (is.null(rv$click_lat)) return(
      div(style="color:var(--sage); font-family:var(--font); font-size:11px; line-height:1.5;",
          "Click the map to select a site. Top 10 species by composite score appear in the Decision-Making tab.")
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
      setView(lng=-96, lat=38, zoom=POPUP_ZOOM)
  })
  output$terrain_minimap <- renderLeaflet({
    leaflet(options=leafletOptions(zoomControl=FALSE, attributionControl=FALSE,
                                   dragging=FALSE, scrollWheelZoom=FALSE, doubleClickZoom=FALSE, boxZoom=FALSE, keyboard=FALSE)) %>%
      addProviderTiles("OpenTopoMap", options=providerTileOptions(maxZoom=17)) %>%
      setView(lng=-96, lat=38, zoom=POPUP_ZOOM)
  })
  
  # ── Visualizations tab map ──────────────────────────────────────────────────
  output$vis_map <- renderLeaflet({
    leaflet(options=leafletOptions(zoomControl=TRUE)) %>%
      addProviderTiles("CartoDB.Positron",
                       options=providerTileOptions(maxZoom=18)) %>%
      addPolygons(data=state_sf, fillOpacity=0, color="#2D5A48",
                  weight=0.8, opacity=0.65) %>%
      setView(lng=-96, lat=38, zoom=4)
  })
  
  output$vis_species_list_ui <- renderUI({
    keys_available <- names(all_vis_display)
    sel <- input$vis_selected_species
    tagList(lapply(seq_along(all_vis_display), function(i) {
      key <- keys_available[i]; nm <- all_vis_display[i]
      cls <- if (!is.null(sel) && sel == key) "vis-sp-item active" else "vis-sp-item"
      div(class=cls, `data-key`=key, nm)
    }))
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
  
  # CHM (canopy) layer: species-independent — render immediately when selected,
  # regardless of whether a species has been picked.
  observeEvent(input$vis_active_layer, {
    layer <- input$vis_active_layer
    if (!is.null(layer) && layer == "canopy") {
      rv$vis_raster       <- NULL
      rv$vis_fallback     <- FALSE
      rv$vis_raster_range <- c(NA, NA)
      # render_vis_layers is also called by the vis_active_layer observer below,
      # but we trigger it here too to be explicit and immediate.
      if (identical(input$main_tabs, "Visualizations"))
        session$onFlushed(function() { render_vis_layers() }, once=TRUE)
    }
  }, ignoreNULL=TRUE, ignoreInit=FALSE)
  
  # BA layers: need a species key — only load when one is selected.
  observeEvent(list(input$vis_selected_species, input$time_range, input$vis_active_layer), {
    key    <- input$vis_selected_species
    period <- if (!is.null(input$time_range)) input$time_range else "2040_2069"
    layer  <- if (!is.null(input$vis_active_layer)) input$vis_active_layer else "future_ba"
    
    # Canopy needs no species — handled by its own observer above
    if (layer == "canopy") return()
    # BA layers DO need a species
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
    }
  }, ignoreNULL=FALSE)
  
  # Orange gradient palette used for all raster BA layers
  orange_pal <- function(rng) {
    colorNumeric(
      palette  = c("#FFF8EC","#FDDBA5","#F59E0B","#B45309","#7C2D12"),
      domain   = rng,
      na.color = "transparent"
    )
  }
  
  # ── Named layer IDs so we swap only the data layer, never the basemap ─────────
  VIS_RASTER_ID <- "vis_ba_raster"
  VIS_TILE_ID   <- "vis_canopy_tile"
  # 32 MB budget lets leaflet render at full detail; it downsamples for the
  # current viewport automatically, giving full detail at zoom and good
  # continental coverage without manual aggregation.
  VIS_MAX_BYTES <- 32L * 1024L * 1024L
  
  render_vis_layers <- function() {
    layer <- isolate(input$vis_active_layer)
    if (is.null(layer)) layer <- "future_ba"
    
    proxy <- leafletProxy("vis_map")
    
    # Remove only the previous data layer(s), never the basemap or polygons
    proxy <- proxy %>%
      removeImage(VIS_RASTER_ID) %>%
      removeTiles(VIS_TILE_ID)
    
    # Re-add pin if needed (pin is in its own group, unaffected by above)
    lng <- isolate(rv$click_lng)
    lat <- isolate(rv$click_lat)
    if (!is.null(lng)) {
      proxy <- proxy %>%
        clearGroup("vis_pin") %>%
        addMarkers(lng=lng, lat=lat, icon=thumbtack_icon, group="vis_pin")
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
  
  # ── Legend: single layer, real numeric ticks ──────────────────────────────────
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
    
    title <- if (layer == "future_ba") "Predicted Future BA (m\u00b2/ha)"
    else "Historical BA (m\u00b2/ha)"
    
    # 5 evenly-spaced tick values
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
    fmt <- function(x) if (is.na(x)) "—" else sprintf("%.3f", x)
    
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
                         tags$span(class="dm-comp-label","Wildlife Value"),
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
    } else if (!is.na(row$n_components[1]) && row$n_components[1] < 3) {
      paste0("based on ", row$n_components[1], " of 3 components")
    } else "weighted composite of 3 components"
    
    gjam_raw_note <- if (isTRUE(row$gjam_is_genus_avg[1])) {
      div(class="genus-avg-badge", "\u00f8 genus-level average")
    } else div(class="dm-comp-raw", "direct model prediction")
    
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
                tags$span(class="dm-comp-label", "Wildlife Value"),
                div(class="dm-comp-val", fmt(row$wildlife[1])),
                div(class="dm-comp-raw",
                    if (!is.na(row$wildlife_raw[1])) paste0("raw: ", row$wildlife_raw[1], " / 100") else "no data")
            ),
            div(class="dm-comp-card timb-card",
                tags$span(class="dm-comp-label", "Timber Value"),
                div(class="dm-comp-val", fmt(row$timber[1])),
                div(class="dm-comp-raw",
                    if (!is.na(row$timber_raw[1])) paste0("$", formatC(row$timber_raw[1], format="f", digits=0), "/MBF") else "no price data")
            )
        )
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
      w_pct  <- if (!is.na(row$wildlife)) row$wildlife*100 else 0
      t_pct  <- if (!is.na(row$timber))  row$timber*100  else 0
      div(class=if(is_sel) "bar-row hl" else "bar-row", style="align-items:center;",
          div(class="bar-sp-name", row$display_name),
          div(style="flex:1;display:flex;gap:4px;",
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--gjam-c);opacity:0.75;", g_pct))),
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--wild-c);opacity:0.75;", w_pct))),
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--timb-c);opacity:0.75;", t_pct)))
          )
      )
    })
    tagList(header, rows)
  })
  
}

shinyApp(ui, server)
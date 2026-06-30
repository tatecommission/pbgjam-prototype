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
  "~/Documents/PBGJAM-data-explorer/raw_data/v2/BA_wildlife_traits.csv",
  show_col_types = FALSE
) |>
  select(
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
  select(-has_pest_flag_raw)

wildlife_range <- range(wildlife_data$wildlife_value_index, na.rm = TRUE)

coords_lon <- xdata$lon
coords_lat <- xdata$lat
BA_mat     <- BA

timber <- read_csv(
  "~/Documents/PBGJAM-data-explorer/USFScut-sold-reports/timber_price_by_gjam_species_wide.csv",
  show_col_types = FALSE
)

fmt_species <- function(nm) {
  parts <- strsplit(nm, "(?<=[a-z])(?=[A-Z])", perl = TRUE)[[1]]
  if (length(parts) < 2)
    return(paste0(toupper(substr(nm,1,1)), tolower(substr(nm,2,nchar(nm)))))
  genus   <- paste0(toupper(substr(parts[1],1,1)), tolower(substr(parts[1],2,nchar(parts[1]))))
  epithet <- tolower(paste(parts[-1], collapse=""))
  paste(genus, epithet)
}

species_names   <- colnames(BA_mat)
species_display <- setNames(sapply(species_names, fmt_species, USE.NAMES=FALSE), species_names)

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
  readRDS("~/Documents/PBGJAM-data-explorer/raw_data/v2/usfs_regions.rds"),
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
    v <- row[[reg]][1]; if (is.na(v)||!is.finite(v)) NA_real_ else as.numeric(v)
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

message("Building GJAM raster stacks...")
.gjam_stack <- lapply(c("2040_2069","2070_2099"), function(period) {
  tifs <- unname(gjam_files[[period]])
  tifs <- tifs[file.exists(tifs)]
  if (length(tifs) == 0) return(NULL)
  s   <- terra::rast(tifs)
  nms <- sub(paste0("_",period,"_rcp45\\.tif$"), "",
             sub("^mean_","", basename(terra::sources(s))))
  names(s) <- nms
  s
})
names(.gjam_stack) <- c("2040_2069","2070_2099")

message("  2040-2069: ", terra::nlyr(.gjam_stack[["2040_2069"]]), " layers")
message("  2070-2099: ", terra::nlyr(.gjam_stack[["2070_2099"]]), " layers")

gjam_range <- local({
  s    <- .gjam_stack[["2040_2069"]]
  vals <- unlist(lapply(seq_len(terra::nlyr(s)), function(i) {
    v <- terra::values(s[[i]], na.rm=TRUE); v[is.finite(v)]
  }))
  c(min(vals), max(vals))
})
message("  GJAM value range: ", round(gjam_range[1],3), " to ", round(gjam_range[2],3))

timber_range <- local({
  pc <- c("R1","R2","R3","R4","R5","R6","R8","R9","R10")
  ap <- as.numeric(unlist(timber[,intersect(pc,names(timber))],use.names=FALSE))
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
    score[i] <- if (!any(ok)) NA_real_ else sum(vals[ok]*wts[ok]) / sum(wts[ok])
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

SEARCH_RADIUS_KM <- 50
IDW_POWER <- 1

get_all_species <- function(lng, lat, radius_km = SEARCH_RADIUS_KM) {
  d_km <- haversine_km(lng, lat, coords_lon, coords_lat)
  idxs <- which(d_km <= radius_km)
  if (length(idxs) == 0) {
    idxs <- which.min(d_km)
  }
  w <- 1 / (d_km[idxs]^IDW_POWER)
  w[!is.finite(w)] <- max(w[is.finite(w)], na.rm = TRUE) * 1e6
  w <- w / sum(w)
  
  sub_mat <- BA_mat[idxs, , drop = FALSE]
  sub_mat[is.na(sub_mat)] <- 0
  ba_agg <- as.numeric(w %*% sub_mat)
  names(ba_agg) <- colnames(sub_mat)
  
  data.frame(key=names(ba_agg), display_name=species_display[names(ba_agg)],
             ba=ba_agg, stringsAsFactors=FALSE)
}

# Returns a list with $values (named numeric vector) and $fallback_keys
# (character vector of species whose value came from a genus average rather
# than their own raster, so the UI can flag imprecise scores instead of
# presenting them as equally precise to direct matches).
get_gjam_vals <- function(lng, lat, keys, period) {
  stack <- .gjam_stack[[period]]
  if (is.null(stack))
    return(list(values=setNames(rep(NA_real_, length(keys)), keys), fallback_keys=character(0)))
  pt       <- terra::vect(matrix(c(lng, lat), ncol=2), crs="EPSG:4326")
  extracted <- tryCatch(
    as.numeric(terra::extract(stack, pt)[1, -1]),
    error = function(e) rep(NA_real_, terra::nlyr(stack))
  )
  names(extracted) <- names(stack)
  all_keys <- names(stack)
  fallback_keys <- character(0)
  values <- sapply(keys, function(k) {
    if (k %in% all_keys && is.finite(extracted[[k]])) return(extracted[[k]])
    genus_lower <- tolower(regmatches(k, regexpr("^[a-z]+", k)))
    genus_keys  <- all_keys[startsWith(tolower(all_keys), genus_lower)]
    gv <- extracted[genus_keys]
    gv <- gv[is.finite(gv)]
    if (length(gv) > 0) fallback_keys <<- c(fallback_keys, k)
    if (length(gv) == 0) NA_real_ else mean(gv)
  })
  list(values=values, fallback_keys=fallback_keys)
}

safe_vals <- function(r) { v <- raster::values(r); v[!is.na(v)&is.finite(v)] }

load_display_raster <- function(fpath, agg_fact=1) {
  if (is.null(fpath)||is.na(fpath)||!file.exists(fpath)) return(NULL)
  r <- terra::rast(fpath)
  if (agg_fact>1) r <- terra::aggregate(r, fact=agg_fact, fun="mean")
  raster::raster(r)
}

load_gjam_raster_for_key <- function(key, period, agg_fact=3) {
  files <- gjam_files[[period]]
  if (key %in% names(files) && file.exists(files[[key]]))
    return(list(raster=load_display_raster(files[[key]], agg_fact), fallback=FALSE))
  genus_lower <- tolower(regmatches(key, regexpr("^[a-z]+", key)))
  genus_keys  <- names(files)[startsWith(tolower(names(files)), genus_lower)]
  genus_keys  <- genus_keys[sapply(genus_keys, function(gk) file.exists(files[[gk]]))]
  if (length(genus_keys) == 0) return(list(raster=NULL, fallback=TRUE))
  tryCatch({
    stack_r <- terra::rast(lapply(genus_keys, function(gk) terra::rast(files[[gk]])))
    r_mean  <- terra::app(stack_r, fun="mean", na.rm=TRUE)
    if (agg_fact > 1) r_mean <- terra::aggregate(r_mean, fact=agg_fact, fun="mean")
    list(raster=raster::raster(r_mean), fallback=TRUE, n_genus=length(genus_keys))
  }, error=function(e) list(raster=NULL, fallback=TRUE))
}

all_gjam_keys <- unique(c(names(gjam_files[["2040_2069"]]), names(gjam_files[["2070_2099"]])))
all_gjam_display <- sort(species_display[intersect(all_gjam_keys, names(species_display))])

# STYLE GUIDE - PBGJAM v2

app_css <- "
:root {
  --pb-forest:     #1C3A28;
  --pb-canopy:     #2A4F38;
  --pb-understory: #3D6B50;
  --pb-sage:       #7DB89A;
  --pb-mist:       #C8DDD4;
  --pb-parchment:  #F4F7F5;
  --pb-linen:      #E8EEE9;
  --pb-pine-line:  #B4C8BC;
  --pb-accent:     #3A8A5A;
  --pb-accent-hi:  #2E7D4A;
  --pb-timber:     #D9A300;
  --pb-wildlife:   #D43B3B;
  --pb-gjam:       #1F6FCC;
  --pb-danger:     #8B1A1A;
  --pb-riparian:   #2E8FA3;
  --pb-ink:        #1A2E24;
  --pb-ink-muted:  #5A7A68;
  --font-display: 'Gill Sans', 'Gill Sans MT', Calibri, Arial, sans-serif;
  --font-body:    Arial, 'Helvetica Neue', sans-serif;
}

* { box-sizing: border-box; }
body { font-family: var(--font-body); background: var(--pb-forest); color: var(--pb-ink); margin: 0; padding: 0; }

.app-header {
  background: var(--pb-forest);
  padding: 10px 24px;
  display: flex; align-items: center; justify-content: space-between;
  border-bottom: 2px solid #0E1F18;
}
.header-logo img { height: 36px; }
.header-sponsors { display: flex; align-items: center; gap: 22px; }
.header-sponsors img { height: 28px; opacity: 0.85; filter: brightness(0) invert(1); }
.header-sponsors .nasa-logo { height: 44px; }

.nav-tabs {
  background: var(--pb-canopy) !important;
  border-bottom: 2px solid #0C1C16 !important;
  padding: 0 16px; margin: 0 !important;
}
.nav-tabs > li > a {
  font-family: var(--font-display);
  font-size: 12px; font-weight: 700; letter-spacing: 0.06em; text-transform: uppercase;
  color: var(--pb-sage) !important;
  border: 1px solid rgba(0,0,0,0.3) !important; border-bottom: none !important;
  border-radius: 0 !important;
  padding: 8px 20px !important; margin: 5px 2px 0 !important;
  background: var(--pb-canopy) !important;
  transition: background 0.1s, color 0.1s;
}
.nav-tabs > li.active > a, .nav-tabs > li.active > a:focus {
  color: var(--pb-forest) !important;
  background: var(--pb-parchment) !important;
  border-color: var(--pb-pine-line) !important;
}
.nav-tabs > li > a:hover {
  color: var(--pb-mist) !important;
  background: var(--pb-understory) !important;
}
.tab-content { padding: 0; }

.welcome-tab {
  height: calc(100vh - 110px);
  background-image: url('forest-background.jpg');
  background-size: cover; background-position: center 40%;
  display: flex; align-items: center; justify-content: center;
}
.welcome-panel {
  background: rgba(6,18,11,0.82);
  padding: 56px 68px;
  max-width: 720px; text-align: center;
  border: 1px solid rgba(80,140,100,0.25);
}
.welcome-title {
  font-family: var(--font-display);
  font-size: 62px; font-weight: 700;
  color: #E8F5EC; letter-spacing: 0.04em;
  margin: 0 0 24px 0; line-height: 1;
  text-shadow: 0 2px 14px rgba(0,0,0,0.95);
}
.welcome-desc {
  font-family: var(--font-display);
  font-size: 19px; color: #A4C4AC; line-height: 1.7;
  font-weight: 400; margin: 0;
  text-shadow: 0 1px 5px rgba(0,0,0,0.75);
}

.tab1-layout { display: flex; height: calc(100vh - 110px); }
.sidebar-col {
  width: 260px; min-width: 260px;
  background: var(--pb-canopy);
  color: var(--pb-mist);
  padding: 18px 14px; overflow-y: auto;
  border-right: 2px solid #0C1C16;
}
.sidebar-section-label {
  font-family: var(--font-display);
  font-size: 10px; font-weight: 700; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--pb-sage); margin: 0 0 10px 0; display: block;
}
.weight-group {
  background: rgba(0,0,0,0.25);
  border: 1px solid rgba(0,0,0,0.3);
  padding: 10px 12px 8px; margin-bottom: 10px;
}
.weight-label {
  font-family: var(--font-body);
  font-size: 11px; font-weight: 700; color: var(--pb-mist); margin-bottom: 4px;
}
.weight-note {
  font-family: var(--font-body);
  font-size: 9px; color: var(--pb-sage); margin-top: 2px; line-height: 1.3;
}
.irs--shiny .irs-bar { background: var(--pb-accent); height: 5px; }
.irs--shiny .irs-line { background: rgba(0,0,0,0.35); height: 5px; }
.irs--shiny .irs-handle {
  background: var(--pb-mist); border-color: var(--pb-sage);
  border-radius: 0 !important; width: 14px !important; height: 14px !important;
  top: 22px !important;
}
.irs--shiny .irs-single {
  background: var(--pb-understory);
  border-radius: 0;
  font-family: var(--font-body); font-size: 10px;
}
.map-col { flex: 1; position: relative; overflow: hidden; }
.map-col #main_map { height: 100%; width: 100%; }
.coord-display {
  position: absolute; bottom: 18px; left: 18px; z-index: 1000;
  background: rgba(28,58,40,0.92);
  color: var(--pb-mist);
  font-family: var(--font-body); font-size: 11px; font-weight: 700;
  padding: 7px 14px;
  border: 1px solid var(--pb-understory);
  pointer-events: none;
}
.map-popups { position: absolute; top: 14px; left: 14px; z-index: 1000; display: flex; gap: 8px; pointer-events: none; }
.map-popup-panel {
  width: 185px; height: 190px;
  border: 1px solid rgba(80,160,110,0.4);
  background: rgba(18,32,24,0.92);
  display: flex; flex-direction: column; pointer-events: auto;
}
.map-popup-header {
  padding: 5px 9px;
  background: rgba(28,58,40,0.95);
  border-bottom: 1px solid rgba(80,160,110,0.25);
  flex-shrink: 0;
}
.map-popup-title { font-family: var(--font-display); font-size: 10px; font-weight: 700; letter-spacing: 0.07em; text-transform: uppercase; color: var(--pb-sage); }
.map-popup-body { flex: 1; position: relative; overflow: hidden; }

.riparian-toggle-group {
  background: rgba(0,0,0,0.25);
  border: 1px solid rgba(0,0,0,0.3);
  padding: 10px 12px; margin-bottom: 4px;
  display: flex; align-items: flex-start; gap: 8px;
}
.riparian-toggle-group input[type='checkbox'] {
  margin-top: 2px; accent-color: var(--pb-riparian); cursor: pointer;
}
.riparian-toggle-label {
  font-family: var(--font-body); font-size: 11px; font-weight: 700;
  color: var(--pb-mist); cursor: pointer; display: block;
}
.riparian-toggle-note {
  font-family: var(--font-body); font-size: 9px; color: var(--pb-sage);
  margin-top: 2px; line-height: 1.3; display: block;
}
.flag-tag {
  display: inline-block;
  font-family: var(--font-body); font-size: 8px; font-weight: 700;
  margin-left: 4px;
  vertical-align: middle;
  color: var(--pb-ink-muted);
}
.flag-tag.flag-riparian { color: var(--pb-riparian); }
.flag-tag.flag-pest { color: var(--pb-danger); }

.gjam-outer { height: calc(100vh - 110px); display: flex; }

.gjam-sidebar {
  width: 230px; min-width: 230px;
  background: var(--pb-parchment);
  border-right: 2px solid var(--pb-pine-line);
  display: flex; flex-direction: column;
}
.gjam-sidebar-header {
  padding: 12px 12px 8px;
  border-bottom: 1px solid var(--pb-pine-line);
  background: var(--pb-linen);
  flex-shrink: 0;
}
.gjam-sidebar-header h4 {
  font-family: var(--font-display);
  font-size: 10px; font-weight: 700; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--pb-understory); margin: 0 0 6px 0;
}
.gjam-search {
  width: 100%; padding: 5px 8px;
  font-family: var(--font-body); font-size: 11px;
  border: 1px solid var(--pb-pine-line);
  background: white; color: var(--pb-ink);
  border-radius: 0;
  outline: none;
}
.gjam-search:focus { border-color: var(--pb-accent); }
.gjam-species-list { flex: 1; overflow-y: auto; padding: 4px 0; }
.gjam-sp-item {
  padding: 7px 12px;
  font-family: var(--font-body); font-style: italic; font-size: 11px;
  color: var(--pb-ink);
  cursor: pointer; border-bottom: 1px solid var(--pb-pine-line);
  transition: background 0.1s;
  white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}
.gjam-sp-item:hover { background: var(--pb-linen); }
.gjam-sp-item.active { background: var(--pb-understory); color: white; font-style: italic; }
.gjam-sp-item.active:hover { background: var(--pb-understory); }

.gjam-map-area { flex: 1; display: flex; flex-direction: column; }
.gjam-topbar {
  display: flex; align-items: center; gap: 14px; padding: 8px 14px; flex-shrink: 0;
  background: var(--pb-linen);
  border-bottom: 1px solid var(--pb-pine-line);
}
.gjam-species-title {
  font-family: var(--font-body); font-style: italic; font-size: 13px;
  color: var(--pb-ink); font-weight: 700; flex: 1; min-width: 0;
  white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}
.gjam-fallback-note {
  font-family: var(--font-body); font-size: 10px; color: var(--pb-timber);
  font-style: normal; flex-shrink: 0;
}
.time-toggle-wrap { display: flex; align-items: center; gap: 8px; }
.time-toggle-label {
  font-family: var(--font-display);
  font-size: 9px; font-weight: 700; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--pb-ink-muted); white-space: nowrap;
}
.time-toggle {
  display: inline-flex;
  border: 1px solid var(--pb-pine-line);
  background: white;
}
.time-seg {
  font-family: var(--font-display);
  font-size: 10px; font-weight: 700; letter-spacing: 0.05em;
  padding: 5px 14px; cursor: pointer; user-select: none;
  color: var(--pb-ink-muted); transition: background 0.1s, color 0.1s;
  border: none; border-right: 1px solid var(--pb-pine-line);
}
.time-seg:last-child { border-right: none; }
.time-seg.active { background: var(--pb-understory); color: white; }
.gjam-map-wrap { flex: 1; position: relative; }
.gjam-map-wrap #gjam_map { height: 100%; width: 100%; }
.gjam-legend {
  position: absolute; bottom: 10px; left: 10px; z-index: 1000;
  background: rgba(244,247,245,0.95);
  border: 1px solid var(--pb-pine-line);
  padding: 7px 10px; min-width: 130px; pointer-events: none;
}
.legend-title { font-family: var(--font-display); font-size: 9px; font-weight: 700; letter-spacing: 0.08em; text-transform: uppercase; color: var(--pb-understory); display: block; margin-bottom: 4px; }
.legend-gradient { height: 9px; width: 100%; }
.legend-labels { display: flex; justify-content: space-between; margin-top: 2px; }
.legend-labels span { font-family: var(--font-body); font-size: 9px; color: var(--pb-ink-muted); }

.dm-outer { height: calc(100vh - 110px); display: flex; flex-direction: column; background: var(--pb-parchment); }
.dm-topbar {
  display: flex; align-items: center; gap: 14px; padding: 8px 14px; flex-shrink: 0;
  background: var(--pb-linen); border-bottom: 1px solid var(--pb-pine-line);
}
.dm-topbar-title {
  font-family: var(--font-display);
  font-size: 11px; font-weight: 700; letter-spacing: 0.09em; text-transform: uppercase;
  color: var(--pb-understory); flex: 1;
}
.dm-loc-text { font-family: var(--font-body); font-size: 10px; color: var(--pb-ink-muted); }
.dm-time-controls { display: flex; align-items: center; gap: 10px; }

.dm-main { flex: 1; display: flex; min-height: 0; overflow: hidden; }

.dm-left {
  width: 220px; min-width: 220px;
  background: var(--pb-canopy);
  border-right: 2px solid #0C1C16;
  display: flex; flex-direction: column; overflow: hidden;
}
.dm-left-section {
  padding: 12px 12px 10px;
  border-bottom: 1px solid rgba(0,0,0,0.25);
  flex-shrink: 0;
}
.dm-section-label {
  font-family: var(--font-display);
  font-size: 9px; font-weight: 700; letter-spacing: 0.1em; text-transform: uppercase;
  color: var(--pb-sage); display: block; margin-bottom: 8px;
}
.dm-species-scroll { flex: 1; overflow-y: auto; }
.dm-sp-row {
  display: flex; align-items: center; gap: 6px;
  padding: 7px 12px;
  border-bottom: 1px solid rgba(0,0,0,0.15);
  cursor: pointer; transition: background 0.1s;
}
.dm-sp-row:hover { background: rgba(255,255,255,0.06); }
.dm-sp-row.selected { background: rgba(61,107,80,0.6); border-left: 3px solid var(--pb-sage); padding-left: 9px; }
.dm-sp-rank { font-family: var(--font-body); font-size: 9px; font-weight: 700; color: var(--pb-sage); width: 16px; flex-shrink: 0; }
.dm-sp-name { font-family: var(--font-body); font-style: italic; font-size: 11px; color: var(--pb-mist); flex: 1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.dm-sp-row.selected .dm-sp-name { color: white; }
.dm-sp-score { font-family: var(--font-body); font-size: 10px; font-weight: 700; color: var(--pb-sage); }
.dm-sp-row.selected .dm-sp-score { color: #A8D8C0; }
.dm-no-data { padding: 20px 12px; font-family: var(--font-body); font-size: 11px; color: var(--pb-sage); text-align: center; }

.dm-center {
  flex: 1; min-width: 0; display: flex; flex-direction: column;
  border-right: 1px solid var(--pb-pine-line);
  background: var(--pb-parchment);
  overflow: hidden;
}
.dm-panel-header {
  padding: 9px 14px; flex-shrink: 0;
  background: var(--pb-linen); border-bottom: 1px solid var(--pb-pine-line);
  display: flex; align-items: baseline; gap: 10px;
}
.dm-panel-title {
  font-family: var(--font-display);
  font-size: 10px; font-weight: 700; letter-spacing: 0.09em; text-transform: uppercase;
  color: var(--pb-understory);
}
.dm-panel-subtitle { font-family: var(--font-body); font-size: 10px; color: var(--pb-ink-muted); }
.dm-panel-body { flex: 1; overflow-y: auto; padding: 12px 14px; display: flex; flex-direction: column; gap: 12px; }

.dm-metric-row { display: grid; grid-template-columns: repeat(4, 1fr); gap: 8px; }
.dm-metric-card {
  background: white;
  border: 1px solid var(--pb-pine-line);
  border-top: 3px solid var(--pb-pine-line);
  padding: 10px 12px;
}
.dm-metric-card.gjam-card  { border-top-color: var(--pb-gjam); }
.dm-metric-card.wild-card  { border-top-color: var(--pb-wildlife); }
.dm-metric-card.timb-card  { border-top-color: var(--pb-timber); }
.dm-metric-card.total-card { border-top-color: var(--pb-accent-hi); }
.dm-metric-label { font-family: var(--font-display); font-size: 9px; font-weight: 700; letter-spacing: 0.08em; text-transform: uppercase; color: var(--pb-ink-muted); display: block; margin-bottom: 4px; }
.dm-metric-val { font-family: var(--font-body); font-size: 22px; font-weight: 700; color: var(--pb-ink); line-height: 1; }
.dm-metric-card.gjam-card  .dm-metric-val { color: var(--pb-gjam); }
.dm-metric-card.wild-card  .dm-metric-val { color: var(--pb-wildlife); }
.dm-metric-card.timb-card  .dm-metric-val { color: var(--pb-timber); }
.dm-metric-card.total-card .dm-metric-val { color: var(--pb-accent-hi); }
.dm-metric-raw { font-family: var(--font-body); font-size: 10px; color: var(--pb-ink-muted); margin-top: 2px; }

.dm-chart-block {
  background: white;
  border: 1px solid var(--pb-pine-line);
  padding: 10px 12px;
}
.dm-chart-label {
  font-family: var(--font-display);
  font-size: 9px; font-weight: 700; letter-spacing: 0.08em; text-transform: uppercase;
  color: var(--pb-understory); display: block; margin-bottom: 8px;
}
.bar-row { display: flex; align-items: center; gap: 6px; margin-bottom: 5px; }
.bar-row:last-child { margin-bottom: 0; }
.bar-sp-name { font-family: var(--font-body); font-style: italic; font-size: 10px; color: var(--pb-ink); width: 110px; flex-shrink: 0; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.bar-track { flex: 1; height: 12px; background: var(--pb-linen); overflow: hidden; border: 1px solid var(--pb-pine-line); }
.bar-fill { height: 100%; transition: width 0.25s; }
.bar-val { font-family: var(--font-body); font-size: 9px; font-weight: 700; color: var(--pb-ink-muted); width: 34px; flex-shrink: 0; text-align: right; }
.bar-row.hl .bar-sp-name { color: var(--pb-accent-hi); font-weight: 700; }
.bar-row.hl .bar-val { color: var(--pb-accent-hi); }

.dm-right {
  width: 230px; min-width: 230px;
  background: var(--pb-parchment);
  display: flex; flex-direction: column; overflow: hidden;
}
.dm-right-body { flex: 1; overflow-y: auto; padding: 12px 12px; display: flex; flex-direction: column; gap: 10px; }
.subscore-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 7px; }
.subscore-card {
  background: white;
  border: 1px solid var(--pb-pine-line);
  padding: 9px 10px;
}
.subscore-name { font-family: var(--font-display); font-size: 9px; font-weight: 700; letter-spacing: 0.06em; text-transform: uppercase; color: var(--pb-wildlife); display: block; margin-bottom: 3px; }
.subscore-val { font-family: var(--font-body); font-size: 20px; font-weight: 700; color: var(--pb-ink); line-height: 1; }
.subscore-denom { font-family: var(--font-body); font-size: 10px; color: var(--pb-ink-muted); }
.subscore-desc { font-family: var(--font-body); font-size: 9px; color: var(--pb-ink-muted); margin-top: 3px; display: block; line-height: 1.3; }
.pest-alert {
  background: var(--pb-danger);
  padding: 8px 10px; margin-top: 2px;
}
.pest-alert b { display: block; font-family: var(--font-display); font-size: 10px; font-weight: 700; letter-spacing: 0.06em; color: #FFD0D0; margin-bottom: 3px; }
.pest-alert p { margin: 0; font-family: var(--font-body); font-size: 9px; color: #FFD0D0; line-height: 1.4; }

.riparian-info {
  background: white;
  border: 1px solid var(--pb-pine-line);
  border-left: 3px solid var(--pb-riparian);
  padding: 9px 10px;
}
.riparian-info b {
  display: block; font-family: var(--font-display); font-size: 9px; font-weight: 700;
  letter-spacing: 0.07em; text-transform: uppercase; color: var(--pb-riparian); margin-bottom: 4px;
}
.riparian-info p { margin: 0; font-family: var(--font-body); font-size: 10px; color: var(--pb-ink); line-height: 1.4; }
.riparian-info .ripar-no { color: var(--pb-ink-muted); font-style: italic; }

.about-outer { height: calc(100vh - 110px); background: var(--pb-parchment); overflow-y: auto; }
.about-inner { max-width: 900px; margin: 0 auto; padding: 48px 32px 64px; }
.about-heading {
  font-family: var(--font-display); font-size: 32px; font-weight: 700;
  color: var(--pb-forest); margin: 0 0 8px 0; text-align: center;
}
.about-subheading {
  font-family: var(--font-body); font-size: 13px; color: var(--pb-ink-muted);
  text-align: center; margin: 0 0 40px 0; line-height: 1.6;
}
.about-team-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 24px; }
.about-card {
  background: white; border: 1px solid var(--pb-pine-line);
  padding: 20px 18px; text-align: center;
}
.about-photo {
  width: 120px; height: 120px; border-radius: 50%;
  object-fit: cover; margin: 0 auto 14px;
  border: 3px solid var(--pb-linen); display: block;
}
.about-name {
  font-family: var(--font-display); font-size: 15px; font-weight: 700;
  color: var(--pb-forest); margin: 0 0 2px 0;
}
.about-role {
  font-family: var(--font-body); font-size: 10px; font-weight: 700;
  letter-spacing: 0.06em; text-transform: uppercase; color: var(--pb-understory);
  margin: 0 0 10px 0;
}
.about-bio {
  font-family: var(--font-body); font-size: 11.5px; color: var(--pb-ink-muted);
  line-height: 1.55; text-align: left;
}

.no-data-msg { padding: 24px 14px; text-align: center; color: var(--pb-ink-muted); font-family: var(--font-body); font-size: 12px; }
.container-fluid { padding: 0 !important; }
.row { margin: 0 !important; }

.dm-left ::-webkit-scrollbar, .gjam-species-list::-webkit-scrollbar { width: 5px; }
.dm-left ::-webkit-scrollbar-track, .gjam-species-list::-webkit-scrollbar-track { background: rgba(0,0,0,0.2); }
.dm-left ::-webkit-scrollbar-thumb, .gjam-species-list::-webkit-scrollbar-thumb { background: var(--pb-sage); }
"

score_color_css <- function(norm_val, type="composite") {
  if (is.na(norm_val)) return("#C8CCC8")
  if (type == "gjam")     return(sprintf("rgba(31,111,204,%.2f)",  0.3 + norm_val*0.7))
  if (type == "wildlife") return(sprintf("rgba(212,59,59,%.2f)",   0.3 + norm_val*0.7))
  if (type == "timber")   return(sprintf("rgba(217,163,0,%.2f)",   0.3 + norm_val*0.7))
  if (norm_val >= 0.75) return("#2E7D4A")
  if (norm_val >= 0.5)  return("#3A8A5A")
  if (norm_val >= 0.25) return("#5AAA72")
  return("#8AC888")
}

RIPARIAN_BOOST <- 0.12

apply_riparian_boost <- function(score, ba_codes, riparian_on) {
  if (!isTRUE(riparian_on)) return(score)
  is_riparian <- wildlife_data$riparian_recommended[match(ba_codes, wildlife_data$ba_code)]
  is_riparian <- ifelse(is.na(is_riparian), FALSE, is_riparian)
  ifelse(is_riparian, pmin(score + RIPARIAN_BOOST, 1), score)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML(app_css)),
    tags$script(HTML("
$(document).on('click', '.dm-sp-row', function() {
  $('.dm-sp-row').removeClass('selected'); $(this).addClass('selected');
  Shiny.setInputValue('dm_selected_species', $(this).data('key'), {priority:'event'});
});
$(document).on('click', '.time-seg', function() {
  $('.time-seg').removeClass('active'); $(this).addClass('active');
  Shiny.setInputValue('time_range', $(this).data('val'), {priority:'event'});
});
$(document).on('click', '.gjam-sp-item', function() {
  $('.gjam-sp-item').removeClass('active'); $(this).addClass('active');
  Shiny.setInputValue('gjam_selected_species', $(this).data('key'), {priority:'event'});
});
$(document).on('input', '#gjam_search', function() {
  var q = $(this).val().toLowerCase();
  $('.gjam-sp-item').each(function() {
    var name = $(this).text().toLowerCase();
    $(this).toggle(name.indexOf(q) >= 0);
  });
});
    "))
  ),
  
  div(class="app-header",
      div(class="header-logo", tags$img(src="logo.png", alt="PBGJAM")),
      div(class="header-sponsors",
          tags$img(src="duke.png", alt="Duke"),
          tags$img(src="nasa.png", alt="NASA", class="nasa-logo"),
          tags$img(src="neon.png", alt="NEON"),
          tags$img(src="nsf.png",  alt="NSF"))
  ),
  
  tabsetPanel(id="main_tabs",
              
              tabPanel("Home",
                       div(class="welcome-tab",
                           div(class="welcome-panel",
                               tags$h1(class="welcome-title", "PBGJAM v2.0"),
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
                                   div(class="weight-label", "GJAM Climate Score"),
                                   sliderInput("w1", NULL, min=0, max=1, value=0, step=0.01, ticks=FALSE, width="100%"),
                                   div(class="weight-note", "Predicted BA in plots across US by given timeframe")
                               ),
                               div(class="weight-group",
                                   div(class="weight-label", "Wildlife Value"),
                                   sliderInput("w2", NULL, min=0, max=1, value=0, step=0.01, ticks=FALSE, width="100%"),
                                   div(class="weight-note", "Combined scores from fruit, insects, pollinators, and deer")
                               ),
                               div(class="weight-group",
                                   div(class="weight-label", "Timber Market Value"),
                                   sliderInput("w3", NULL, min=0, max=1, value=0, step=0.01, ticks=FALSE, width="100%"),
                                   div(class="weight-note", "Regional price based on official 2025 timber sales directly with USFS.")
                               ),
                               div(class="riparian-toggle-group",
                                   tags$input(type="checkbox", id="riparian_priority"),
                                   tags$label(`for`="riparian_priority",
                                              tags$span(class="riparian-toggle-label", "Prioritize riparian species"),
                                              tags$span(class="riparian-toggle-note",
                                                        "Favors species recommended for streambank and wetland-edge planting.")
                                   )
                               ),
                               tags$hr(style="border-color:rgba(255,255,255,0.1); margin: 14px 0;"),
                               tags$span(class="sidebar-section-label", "Climate Period"),
                               div(class="time-toggle",
                                   div(class="time-seg active", `data-val`="2040_2069", "2040-2069"),
                                   div(class="time-seg",        `data-val`="2070_2099", "2070-2099")
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
              
              tabPanel("GJAM",
                       div(class="gjam-outer",
                           
                           div(class="gjam-sidebar",
                               div(class="gjam-sidebar-header",
                                   tags$h4("Species"),
                                   tags$input(id="gjam_search", class="gjam-search",
                                              type="text", placeholder="Filter species...")
                               ),
                               div(class="gjam-species-list",
                                   uiOutput("gjam_species_list_ui")
                               )
                           ),
                           
                           div(class="gjam-map-area",
                               div(class="gjam-topbar",
                                   uiOutput("gjam_species_title_ui"),
                                   div(class="time-toggle-wrap",
                                       div(class="time-toggle-label", "Period"),
                                       div(class="time-toggle",
                                           div(class="time-seg active", `data-val`="2040_2069", "2040-2069"),
                                           div(class="time-seg",        `data-val`="2070_2099", "2070-2099")
                                       )
                                   )
                               ),
                               div(class="gjam-map-wrap",
                                   leafletOutput("gjam_map", height="100%", width="100%"),
                                   div(class="gjam-legend", uiOutput("gjam_legend_ui"))
                               )
                           )
                       )
              ),
              
              tabPanel("Decision-Making",
                       div(class="dm-outer",
                           
                           div(class="dm-topbar",
                               div(class="dm-topbar-title", "Species Scoring Dashboard"),
                               uiOutput("dm_loc_label"),
                               div(class="dm-time-controls",
                                   div(class="time-toggle-label", "Climate period"),
                                   div(class="time-toggle",
                                       div(class="time-seg active", `data-val`="2040_2069", "2040-2069"),
                                       div(class="time-seg",        `data-val`="2070_2099", "2070-2099")
                                   )
                               ),
                               uiOutput("dm_region_label")
                           ),
                           
                           div(class="dm-main",
                               
                               div(class="dm-left",
                                   div(class="dm-left-section",
                                       tags$span(class="dm-section-label", "Top 10 by composite score")
                                   ),
                                   div(class="dm-species-scroll",
                                       uiOutput("dm_species_list")
                                   )
                               ),
                               
                               div(class="dm-center",
                                   div(class="dm-panel-header",
                                       div(class="dm-panel-title", "Score Breakdown"),
                                       uiOutput("dm_selected_name")
                                   ),
                                   div(class="dm-panel-body",
                                       uiOutput("dm_metric_cards"),
                                       div(class="dm-chart-block",
                                           tags$span(class="dm-chart-label", "All 10 species, composite score"),
                                           uiOutput("dm_composite_bars")
                                       ),
                                       div(class="dm-chart-block",
                                           tags$span(class="dm-chart-label", "All 10 species, by component"),
                                           uiOutput("dm_component_bars")
                                       )
                                   )
                               ),
                               
                               div(class="dm-right",
                                   div(class="dm-panel-header",
                                       div(class="dm-panel-title", "Wildlife Detail")
                                   ),
                                   div(class="dm-right-body",
                                       uiOutput("dm_wildlife_subscores"),
                                       uiOutput("dm_riparian_info"),
                                       uiOutput("dm_pest_info")
                                   )
                               )
                           )
                       )
              ),
              
              tabPanel("About",
                       div(class="about-outer",
                           div(class="about-inner",
                               div(class="about-heading", "About PBGJAM"),
                               div(class="about-subheading",
                                   "PBGJAM v2 was developed by Tate Commission, Dr. Tong Qiu, and ",
                                   "Dr. James Clark from Duke University."
                               ),
                               div(class="about-team-grid",
                                   div(class="about-card",
                                       tags$img(class="about-photo", src="tatecommission.png", alt="Tate Commission"),
                                       div(class="about-name", "Tate Commission"),
                                       div(class="about-role", "Lead Developer"),
                                       div(class="about-bio",
                                           "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod ",
                                           "tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim ",
                                           "veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea ",
                                           "commodo consequat."
                                       )
                                   ),
                                   div(class="about-card",
                                       tags$img(class="about-photo", src="tongqiu.png", alt="Dr. Tong Qiu"),
                                       div(class="about-name", "Dr. Tong Qiu"),
                                       div(class="about-role", "Faculty Advisor"),
                                       div(class="about-bio",
                                           "Duis aute irure dolor in reprehenderit in voluptate velit esse cillum ",
                                           "dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non ",
                                           "proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
                                       )
                                   ),
                                   div(class="about-card",
                                       tags$img(class="about-photo", src="jamesclark.png", alt="Dr. James Clark"),
                                       div(class="about-name", "Dr. James Clark"),
                                       div(class="about-role", "Faculty Advisor"),
                                       div(class="about-bio",
                                           "Sed ut perspiciatis unde omnis iste natus error sit voluptatem ",
                                           "accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ",
                                           "ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt."
                                       )
                                   )
                               )
                           )
                       )
              )
              
  )
)

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    click_lng       = NULL,
    click_lat       = NULL,
    time_range      = "2040_2069",
    gjam_raster     = NULL,
    gjam_fallback   = FALSE,
    usfs_region     = NULL
  )
  
  observeEvent(input$time_range, {
    rv$time_range <- input$time_range
  })
  
  composite_scores <- reactive({
    if (is.null(rv$click_lng)) return(NULL)
    all_sp <- tryCatch(get_all_species(rv$click_lng, rv$click_lat), error=function(e) NULL)
    if (is.null(all_sp)) return(NULL)
    keys <- all_sp$key
    
    gjam_result    <- get_gjam_vals(rv$click_lng, rv$click_lat, keys, rv$time_range)
    gjam_raw       <- gjam_result$values
    gjam_n         <- norm_global(gjam_raw, gjam_range)
    gjam_fallback_keys <- gjam_result$fallback_keys
    
    wild_raw  <- wildlife_data$wildlife_value_index[match(keys, wildlife_data$ba_code)]
    wild_n    <- norm_global(wild_raw, wildlife_range)
    
    timber_prices <- sapply(keys, function(k) get_timber_price(rv$click_lng, rv$click_lat, k))
    timber_n      <- norm_global(timber_prices, timber_range)
    
    composite <- weighted_composite_full(gjam_n, wild_n, timber_n, input$w1, input$w2, input$w3)
    
    final_score <- apply_riparian_boost(composite$score, keys, input$riparian_priority)
    
    df <- data.frame(
      key            = keys,
      display_name   = all_sp$display_name,
      gjam           = round(gjam_n, 4),
      wildlife       = round(wild_n, 4),
      timber         = round(timber_n, 4),
      score          = round(final_score, 4),
      n_components   = composite$n_components,
      timber_raw     = round(timber_prices, 2),
      wildlife_raw   = round(wild_raw, 1),
      gjam_is_genus_avg = keys %in% gjam_fallback_keys,
      stringsAsFactors = FALSE
    )
    df[order(df$score, decreasing=TRUE, na.last=TRUE), ]
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
              .setContent('<span style=\\'font-family:Arial;font-size:12px;color:#555;\\'>Select a location within the contiguous United States.</span>')
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
    leafletProxy("gjam_map") %>% clearMarkers() %>%
      addMarkers(lng=lng, lat=lat, icon=thumbtack_icon)
    
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
      div(style="color:var(--pb-sage); font-family:var(--font-body); font-size:11px;",
          "Click the map to select a site. The top 10 species by composite score will appear in the Decision-Making tab.")
    )
    df <- top10_scored()
    if (is.null(df)) return(NULL)
    div(
      tags$span(class="sidebar-section-label", style="margin-top:0;", "Top species at site"),
      tagList(lapply(seq_len(min(5,nrow(df))), function(i) {
        row <- df[i,]
        wd  <- wildlife_data[wildlife_data$ba_code == row$key, ]
        is_riparian <- nrow(wd) > 0 && isTRUE(wd$riparian_recommended[1])
        div(style="display:flex;align-items:center;gap:6px;padding:4px 0;border-bottom:1px solid rgba(255,255,255,0.08);",
            div(style="font-family:var(--font-body);font-size:9px;font-weight:700;color:var(--pb-sage);width:14px;", i),
            div(style="font-family:var(--font-body);font-style:italic;font-size:11px;color:var(--pb-mist);flex:1;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;",
                row$display_name,
                if (is_riparian) tags$span(class="flag-tag flag-riparian", "RP") else NULL),
            div(style="font-family:var(--font-body);font-size:10px;font-weight:700;color:var(--pb-accent);", sprintf("%.3f",row$score))
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
  
  output$gjam_map <- renderLeaflet({
    leaflet(options=leafletOptions(zoomControl=TRUE)) %>%
      addProviderTiles("CartoDB.Positron", options=providerTileOptions(maxZoom=12)) %>%
      addPolygons(data=state_sf, fillOpacity=0, color="#2D5A48", weight=0.8, opacity=0.65) %>%
      setView(lng=-96, lat=38, zoom=4)
  })
  
  output$gjam_species_list_ui <- renderUI({
    keys_available <- names(all_gjam_display)
    sel <- input$gjam_selected_species
    tagList(lapply(seq_along(all_gjam_display), function(i) {
      key  <- keys_available[i]
      nm   <- all_gjam_display[i]
      cls  <- if (!is.null(sel) && sel == key) "gjam-sp-item active" else "gjam-sp-item"
      div(class=cls, `data-key`=key, nm)
    }))
  })
  
  output$gjam_species_title_ui <- renderUI({
    key <- input$gjam_selected_species
    if (is.null(key) || !key %in% names(species_display))
      return(div(class="gjam-species-title",
                 style="font-style:normal;color:var(--pb-ink-muted);font-size:11px;",
                 "Select a species to view its climate projection"))
    fb_note <- if (isTRUE(rv$gjam_fallback)) {
      div(class="gjam-fallback-note",
          "\u26a0 Genus average shown. No direct raster for this species.")
    } else NULL
    tagList(
      div(class="gjam-species-title", species_display[key]),
      fb_note
    )
  })
  
  observeEvent(list(input$gjam_selected_species, input$time_range), {
    key    <- input$gjam_selected_species
    period <- if (!is.null(input$time_range)) input$time_range else "2040_2069"
    if (is.null(key) || key == "") return()
    result <- load_gjam_raster_for_key(key, period, agg_fact=3)
    rv$gjam_raster  <- result$raster
    rv$gjam_fallback <- isTRUE(result$fallback)
  }, ignoreNULL=TRUE)
  
  render_gjam_layer <- function() {
    lng    <- isolate(rv$click_lng)
    lat    <- isolate(rv$click_lat)
    raster <- isolate(rv$gjam_raster)
    
    proxy <- leafletProxy("gjam_map") %>% clearImages() %>% clearMarkers()
    if (!is.null(lng))
      proxy <- proxy %>% addMarkers(lng=lng, lat=lat, icon=thumbtack_icon)
    if (!is.null(raster)) {
      vals <- safe_vals(raster)
      if (length(vals) > 0) {
        rng <- range(vals)
        pal <- colorNumeric(
          palette  = c("#FFF8E1","#FFB300","#E65100","#8B0000"),
          domain   = rng,
          na.color = "transparent"
        )
        proxy <- proxy %>%
          addRasterImage(raster, colors=pal, opacity=0.8, project=TRUE)
      }
    }
  }
  
  observeEvent(input$main_tabs, {
    if (identical(input$main_tabs, "GJAM")) {
      session$onFlushed(function() { render_gjam_layer() }, once=TRUE)
    }
  })
  observeEvent(rv$gjam_raster, {
    if (identical(input$main_tabs, "GJAM")) render_gjam_layer()
  }, ignoreNULL=TRUE)
  
  output$gjam_legend_ui <- renderUI({
    r <- rv$gjam_raster
    if (is.null(r)) return(NULL)
    vals <- safe_vals(r); if (length(vals)==0) return(NULL)
    rng <- range(vals)
    tagList(
      tags$span(class="legend-title", "Predicted basal area"),
      div(class="legend-gradient",
          style="background:linear-gradient(to right,#FFF8E1,#FFB300,#E65100,#8B0000);"),
      div(class="legend-labels",
          tags$span(round(rng[1],2)),
          tags$span(round(rng[2],2)))
    )
  })
  
  output$dm_loc_label <- renderUI({
    if (is.null(rv$click_lat))
      div(class="dm-loc-text", "No site selected. Click the map on the Location tab.")
    else
      div(class="dm-loc-text", sprintf("%.4f\u00b0N  %.4f\u00b0W", rv$click_lat, abs(rv$click_lng)))
  })
  
  output$dm_region_label <- renderUI({
    if (!is.null(rv$usfs_region))
      div(class="dm-loc-text", paste0("USFS ", rv$usfs_region))
    else
      NULL
  })
  
  output$dm_species_list <- renderUI({
    df  <- top10_scored()
    sel <- input$dm_selected_species
    if (is.null(df))
      return(div(class="dm-no-data", "Select a site on the Location tab to see species rankings here."))
    
    items <- lapply(seq_len(nrow(df)), function(i) {
      row <- df[i,]
      cls <- if (!is.null(sel) && identical(row$key, sel)) "dm-sp-row selected" else "dm-sp-row"
      wd  <- wildlife_data[wildlife_data$ba_code == row$key, ]
      pest_flag     <- nrow(wd) > 0 && isTRUE(wd$has_pest_flag[1])
      riparian_flag <- nrow(wd) > 0 && isTRUE(wd$riparian_recommended[1])
      div(class=cls, `data-key`=row$key,
          div(class="dm-sp-rank", i),
          div(class="dm-sp-name", row$display_name,
              if (isTRUE(input$riparian_priority) && riparian_flag)
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
    if (is.null(key) || !key %in% names(species_display)) return(div(class="dm-panel-subtitle", "-"))
    div(class="dm-panel-subtitle", species_display[key])
  })
  
  output$dm_metric_cards <- renderUI({
    df  <- top10_scored()
    key <- input$dm_selected_species
    if (is.null(df) || is.null(key)) {
      return(div(class="dm-metric-row",
                 div(class="dm-metric-card gjam-card",  div(class="dm-metric-label","GJAM"), div(class="dm-metric-val","-")),
                 div(class="dm-metric-card wild-card",  div(class="dm-metric-label","Wildlife"), div(class="dm-metric-val","-")),
                 div(class="dm-metric-card timb-card",  div(class="dm-metric-label","Timber"), div(class="dm-metric-val","-")),
                 div(class="dm-metric-card total-card", div(class="dm-metric-label","Total"), div(class="dm-metric-val","-"))
      ))
    }
    row <- df[df$key == key, ]
    if (nrow(row)==0) return(NULL)
    fmt <- function(x) if (is.na(x)) "-" else sprintf("%.3f", x)
    riparian_note <- if (isTRUE(input$riparian_priority)) {
      wd <- wildlife_data[wildlife_data$ba_code == key, ]
      if (nrow(wd) > 0 && isTRUE(wd$riparian_recommended[1]))
        paste0("+", RIPARIAN_BOOST, " riparian boost applied")
      else "no riparian boost (not recommended)"
    } else NULL
    div(class="dm-metric-row",
        div(class="dm-metric-card gjam-card",
            div(class="dm-metric-label", "GJAM"),
            div(class="dm-metric-val", fmt(row$gjam[1])),
            div(class="dm-metric-raw",
                if (isTRUE(row$gjam_is_genus_avg[1]))
                  "genus average (no direct data for this species)"
                else
                  "climate suitability (normalized)")
        ),
        div(class="dm-metric-card wild-card",
            div(class="dm-metric-label", "Wildlife"),
            div(class="dm-metric-val", fmt(row$wildlife[1])),
            div(class="dm-metric-raw",
                if (!is.na(row$wildlife_raw[1])) paste0("raw: ", row$wildlife_raw[1], " / 100") else "no data")
        ),
        div(class="dm-metric-card timb-card",
            div(class="dm-metric-label", "Timber"),
            div(class="dm-metric-val", fmt(row$timber[1])),
            div(class="dm-metric-raw",
                if (!is.na(row$timber_raw[1])) paste0("$", formatC(row$timber_raw[1], format="f", digits=0), "/MBF") else "no price data")
        ),
        div(class="dm-metric-card total-card",
            div(class="dm-metric-label", "Composite"),
            div(class="dm-metric-val", fmt(row$score[1])),
            div(class="dm-metric-raw",
                if (!is.null(riparian_note)) riparian_note
                else if (!is.na(row$n_components[1]) && row$n_components[1] < 3)
                  paste0("based on ", row$n_components[1], " of 3 components")
                else
                  "based on 3 of 3 components")
        )
    )
  })
  
  output$dm_composite_bars <- renderUI({
    df  <- top10_scored()
    sel <- input$dm_selected_species
    if (is.null(df)) return(div(class="no-data-msg", "-"))
    tagList(lapply(seq_len(nrow(df)), function(i) {
      row    <- df[i,]
      is_sel <- !is.null(sel) && identical(row$key, sel)
      pct    <- if (!is.na(row$score)) row$score * 100 else 0
      col    <- score_color_css(row$score, "composite")
      div(class=if(is_sel) "bar-row hl" else "bar-row",
          div(class="bar-sp-name", row$display_name),
          div(class="bar-track",
              div(class="bar-fill", style=sprintf("width:%.1f%%;background:%s;", pct, col))),
          div(class="bar-val", if (!is.na(row$score)) sprintf("%.3f",row$score) else "-")
      )
    }))
  })
  
  output$dm_component_bars <- renderUI({
    df  <- top10_scored()
    sel <- input$dm_selected_species
    if (is.null(df)) return(div(class="no-data-msg", "-"))
    
    header <- div(style="display:flex;gap:12px;margin-bottom:6px;",
                  div(style="width:110px;flex-shrink:0;"),
                  div(style="flex:1;display:flex;gap:4px;",
                      div(style="flex:1;font-family:var(--font-display);font-size:8px;font-weight:700;letter-spacing:0.06em;text-transform:uppercase;color:var(--pb-gjam);text-align:center;","GJAM"),
                      div(style="flex:1;font-family:var(--font-display);font-size:8px;font-weight:700;letter-spacing:0.06em;text-transform:uppercase;color:var(--pb-wildlife);text-align:center;","Wildlife"),
                      div(style="flex:1;font-family:var(--font-display);font-size:8px;font-weight:700;letter-spacing:0.06em;text-transform:uppercase;color:var(--pb-timber);text-align:center;","Timber")
                  )
    )
    rows <- lapply(seq_len(nrow(df)), function(i) {
      row    <- df[i,]
      is_sel <- !is.null(sel) && identical(row$key, sel)
      g_pct  <- if (!is.na(row$gjam))     row$gjam*100     else 0
      w_pct  <- if (!is.na(row$wildlife))  row$wildlife*100 else 0
      t_pct  <- if (!is.na(row$timber))   row$timber*100   else 0
      div(class=if(is_sel) "bar-row hl" else "bar-row",
          style="align-items:center;",
          div(class="bar-sp-name", row$display_name),
          div(style="flex:1;display:flex;gap:4px;",
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--pb-gjam);opacity:0.7;", g_pct))),
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--pb-wildlife);opacity:0.7;", w_pct))),
              div(class="bar-track", style="flex:1;",
                  div(class="bar-fill", style=sprintf("width:%.1f%%;background:var(--pb-timber);opacity:0.7;", t_pct)))
          )
      )
    })
    tagList(header, rows)
  })
  
  output$dm_wildlife_subscores <- renderUI({
    key <- input$dm_selected_species
    if (is.null(key)) return(div(class="no-data-msg", "Select a species"))
    wd <- wildlife_data[wildlife_data$ba_code == key, ]
    if (nrow(wd)==0) return(div(class="no-data-msg", "No wildlife data"))
    
    fmt_sub <- function(x) if (is.na(x)) "-" else sprintf("%.1f", x)
    browse_desc <- switch(
      tolower(replace_na(wd$palatable_browse_animal[1],"unknown")),
      "high"   = "High palatability",
      "medium" = "Moderate value",
      "low"    = "Low value",
      "Unknown"
    )
    mast_desc <- paste0(
      replace_na(wd$fruit_seed_abundance[1],"?"), " abundance",
      if (!is.na(wd$berry_nut_seed_product[1]) && tolower(wd$berry_nut_seed_product[1])=="yes") " - mast producer" else ""
    )
    lep_count <- if (!is.na(wd$lep_spp_supported[1])) paste0(wd$lep_spp_supported[1], " Lep spp.") else "No data"
    bloom_desc <- if (!is.na(wd$bloom_period[1]) && nzchar(wd$bloom_period[1])) wd$bloom_period[1] else "No bloom data"
    
    tagList(
      tags$span(style="font-family:var(--font-display);font-size:9px;font-weight:700;letter-spacing:0.1em;text-transform:uppercase;color:var(--pb-wildlife);display:block;margin-bottom:7px;",
                "Wildlife subscores"),
      div(class="subscore-grid",
          div(class="subscore-card",
              tags$span(class="subscore-name","Mast / Fruit"),
              tags$span(class="subscore-val", fmt_sub(wd$cat_A[1])),
              tags$span(class="subscore-denom"," / 25"),
              tags$span(class="subscore-desc", mast_desc)
          ),
          div(class="subscore-card",
              tags$span(class="subscore-name","Insect Host"),
              tags$span(class="subscore-val", fmt_sub(wd$cat_B[1])),
              tags$span(class="subscore-denom"," / 25"),
              tags$span(class="subscore-desc", lep_count)
          ),
          div(class="subscore-card",
              tags$span(class="subscore-name","Pollinator"),
              tags$span(class="subscore-val", fmt_sub(wd$cat_C[1])),
              tags$span(class="subscore-denom"," / 25"),
              tags$span(class="subscore-desc", bloom_desc)
          ),
          div(class="subscore-card",
              tags$span(class="subscore-name","Browse / Deer"),
              tags$span(class="subscore-val", fmt_sub(wd$cat_D[1])),
              tags$span(class="subscore-denom"," / 25"),
              tags$span(class="subscore-desc", browse_desc)
          )
      )
    )
  })
  
  output$dm_riparian_info <- renderUI({
    key <- input$dm_selected_species
    if (is.null(key)) return(NULL)
    wd <- wildlife_data[wildlife_data$ba_code == key, ]
    if (nrow(wd) == 0) return(NULL)
    
    is_recommended <- isTRUE(wd$riparian_recommended[1])
    category       <- wd$riparian_category[1]
    wettest        <- wd$wetland_wettest[1]
    
    if (is_recommended) {
      div(class="riparian-info",
          tags$b("Riparian status"),
          tags$p(
            if (!is.na(category) && nzchar(category)) category else "Recommended for riparian planting",
            if (!is.na(wettest) && nzchar(wettest)) paste0(" (wetland indicator: ", wettest, ")") else ""
          )
      )
    } else {
      div(class="riparian-info",
          tags$b("Riparian status"),
          tags$p(class="ripar-no",
                 if (!is.na(category) && nzchar(category)) category else "Not typically riparian")
      )
    }
  })
  
  output$dm_pest_info <- renderUI({
    key <- input$dm_selected_species
    if (is.null(key)) return(NULL)
    wd <- wildlife_data[wildlife_data$ba_code == key, ]
    if (nrow(wd)==0 || !isTRUE(wd$has_pest_flag[1])) return(NULL)
    div(class="pest-alert",
        tags$b("\u26a0  Pest / Pathogen Alert"),
        tags$p(replace_na(wd$pest_flag_list[1], "No detail available")),
        if (!is.na(wd$pest_flag_notes[1]) && nzchar(wd$pest_flag_notes[1]))
          tags$p(style="margin-top:4px;", wd$pest_flag_notes[1])
        else NULL
    )
  })
  
}

shinyApp(ui, server)
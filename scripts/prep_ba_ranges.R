# prep_ba_ranges.R -- run ONCE after prep_maps_v4.R (RStudio: Source). Takes a few minutes to ~20 min.
# For each species: cap = the highest 99.5th percentile of its NON-ZERO cells over the 5 maps
# (2000, and 2060/2090 under SSP 4.5 and 8.5). The app then uses ONE fixed scale per species,
# 0..cap, for every period and scenario, both for map colours and for the species ranking.
library(terra)
base <- path.expand("~/Documents/PBGJAM-data-explorer/maps_v4/by_species")
period_tags <- c("baseline_2000", "ssp245_2060", "ssp245_2090", "ssp585_2060", "ssp585_2090")
sp   <- sub("\\.tif$", "", list.files(file.path(base, "baseline_2000"), pattern = "\\.tif$"))
stopifnot(length(sp) == 100)
q995 <- function(f) {
  v <- values(rast(f), na.rm = TRUE); v <- v[is.finite(v) & v > 0]
  if (!length(v)) c(0, 0) else c(quantile(v, 0.995, names = FALSE), max(v))
}
res <- data.frame(species = sp, cap = NA_real_, max = NA_real_)
for (i in seq_along(sp)) {
  qs <- sapply(period_tags, function(t) q995(file.path(base, t, paste0(sp[i], ".tif"))))
  res$cap[i] <- max(qs[1, ]); res$max[i] <- max(qs[2, ])
  if (i %% 10 == 0) message(i, " / ", length(sp), " species")
}
write.csv(res, file.path(base, "ba_ranges.csv"), row.names = FALSE)
message("Wrote ", file.path(base, "ba_ranges.csv"), "  (species with cap = 0: ", sum(res$cap == 0), ")")

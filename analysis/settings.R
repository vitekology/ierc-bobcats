# Shared settings and helper functions

EXTENT     <- "data/iNaturalist/study_extent_v3.shp"
CRS_RASTER <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +ellps=sphere +units=m +no_defs"  # US National Atlas Equal Area
R_MAX      <- 0.7
K_MIN      <- 1      # minimum records per 1-km cell for the 2nd-order background

# Data kept outside the repository; set the paths in analysis/paths_local.R
CAM_DIR     <- "external/camera_data"
GROW_DIR    <- "external/grow_cameras"
GDB_DIR     <- "external/landcover_gdb"
RASTER_ROOT <- "external/rasters"
if (file.exists("analysis/paths_local.R")) source("analysis/paths_local.R")

# Values of the raster cells under points (X, Y in EPSG:3310); one column per layer
raster_values <- function(d, order_dir, layers) {
  pts <- terra::vect(as.data.frame(d[c("X", "Y")]), geom = c("X", "Y"), crs = "EPSG:3310")
  xy  <- terra::crds(terra::project(pts, CRS_RASTER))
  vapply(layers, \(l) {
    r <- terra::rast(file.path(RASTER_ROOT, order_dir, paste0(l, ".tif")))
    terra::extract(r, terra::cellFromXY(r, xy))[[1]]
  }, numeric(nrow(d)))
}

# NLCD release and LandTrendr window (10 years, ending 5 years before) for a year
nlcd_sfx <- function(y) ifelse(y < 2013, "11", ifelse(y < 2016, "13", ifelse(y < 2019, "16", "19")))
lt_win   <- function(y) c("97_06", "98_07", "99_08", "00_09", "01_10", "02_11", "03_12", "04_13")[pmin(pmax(y - 2010, 1), 8)]

# Habitat values with the time-varying layers matched to each point's year
habitat_values <- function(d, order_dir, year, static) {
  stopifnot(length(year) == nrow(d), !anyNA(year))
  nl  <- nlcd_sfx(year)
  lt  <- lt_win(year)
  lyr <- c(static,
           paste0(rep(c("d_edge_", "d_scrub_", "d_forest_"), each = length(unique(nl))), unique(nl)),
           paste0("d_disturb_", unique(lt)))
  ex   <- raster_values(d, order_dir, lyr)
  pick <- function(prefix, sfx) ex[cbind(seq_len(nrow(ex)), match(paste0(prefix, sfx), lyr))]
  cbind(ex[, static, drop = FALSE],
        d_edge = pick("d_edge_", nl), d_scrub = pick("d_scrub_", nl),
        d_forest = pick("d_forest_", nl), d_disturb = pick("d_disturb_", lt))
}

# Of two predictors with |r| > r_max, keep the one with the lower univariate AIC
select_predictors <- function(d, candidates, r_max = R_MAX) {
  aic_null <- AIC(glm(case ~ 1, family = binomial(), data = d, weights = weight))
  uni <- data.frame(term = candidates, aic = sapply(candidates, \(v)
    AIC(glm(reformulate(v, "case"), family = binomial(), data = d, weights = weight))))
  uni$delta_aic_vs_null <- round(uni$aic - aic_null, 1)

  r  <- cor(d[candidates])
  ij <- which(abs(r) > r_max & upper.tri(r), arr.ind = TRUE)
  pairs <- data.frame(var1 = rownames(r)[ij[, 1]], var2 = colnames(r)[ij[, 2]], r = round(r[ij], 3),
                      dropped = rep(NA_character_, nrow(ij)))
  pairs <- pairs[order(-abs(pairs$r)), ]

  keep <- candidates
  for (i in seq_len(nrow(pairs))) {
    v1 <- pairs$var1[i]; v2 <- pairs$var2[i]
    if (!(v1 %in% keep && v2 %in% keep)) next
    drop <- if (uni$aic[uni$term == v1] <= uni$aic[uni$term == v2]) v2 else v1
    keep <- setdiff(keep, drop)
    pairs$dropped[i] <- drop
  }
  list(keep = keep, pairs = pairs, univariate = uni[order(uni$aic), c("term", "delta_aic_vs_null")])
}

# GIS layers: data locations (GeoPackage) and prediction maps (GeoTIFF, raster CRS, MAP_RES cells)

library(tidyverse)
library(sf)
library(terra)
source("analysis/settings.R")

MAP_RES <- 90
OUT     <- "output/data_locations/bobcat_layers.gpkg"
map_dir <- "output/maps"
terraOptions(progress = 0)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)

## ---- data locations ---------------------------------------------------------------------------

put  <- function(x, layer) st_write(x, OUT, layer = layer, append = FALSE, quiet = TRUE)
dmin <- function(x) if (all(is.na(x))) as.Date(NA) else min(x, na.rm = TRUE)
dmax <- function(x) if (all(is.na(x))) as.Date(NA) else max(x, na.rm = TRUE)

cams <- read_csv("data/camera/camera_deployments.csv", show_col_types = FALSE) |>
  mutate(last = coalesce(end, start)) |>
  group_by(study, xr = round(X), yr = round(Y)) |>
  summarise(site = paste(sort(unique(site)), collapse = "; "),
            n_deploy = n(), camera_days = sum(effort),
            n_detected = sum(detected), detected = as.integer(n_detected > 0),
            first_date = dmin(start), last_date = dmax(last),
            X = mean(X), Y = mean(Y), .groups = "drop") |>
  select(-xr, -yr) |>
  st_as_sf(coords = c("X", "Y"), crs = 3310)
put(cams, "cameras")

studies <- cams |>
  group_by(study) |>
  summarise(locations = n(), deployments = sum(n_deploy), camera_days = sum(camera_days),
            loc_detected = sum(detected), dep_detected = sum(n_detected),
            first_date = dmin(first_date), last_date = dmax(last_date)) |>
  st_cast("MULTIPOINT")
put(studies, "camera_studies")

used <- read_csv("data/annotated/RSF_data_2nd_order.csv", show_col_types = FALSE) |>
  filter(case == 1) |> pull(id)
inat <- read_csv("data/iNaturalist/inat_bobcat_filtered.csv", show_col_types = FALSE) |>
  transmute(id, observed_on = as.Date(observed_on), year,
            positional_accuracy = as.numeric(positional_accuracy),
            in_model = as.integer(id %in% used), X, Y) |>
  st_as_sf(coords = c("X", "Y"), crs = 3310)
put(inat, "inat_bobcats")

cells <- st_read("output/2nd_order/targetgroup_cells.gpkg", quiet = TRUE) |>
  filter(n_records >= K_MIN)
put(cells, "targetgroup_cells")

pts2 <- read_csv("data/annotated/RSF_data_2nd_order.csv", show_col_types = FALSE) |>
  transmute(case, id, year, X, Y) |>
  st_as_sf(coords = c("X", "Y"), crs = 3310)
put(pts2, "used_available_2nd_order")

akde <- readRDS("output/homeranges/akde_95.rds")
hr <- do.call(rbind, lapply(names(akde), \(i) {
  st_sf(id = i, geometry = st_geometry(akde[[i]]) |> st_transform(3310) |> st_union() |> st_make_valid())
})) |>
  mutate(study_area = sub("_.*$", "", id), poly_km2 = round(as.numeric(st_area(geometry)) / 1e6, 1)) |>
  left_join(read_csv("output/homeranges/homerange_summary.csv", show_col_types = FALSE), by = "id") |>
  st_cast("MULTIPOLYGON")
put(hr, "home_ranges_95")

## ---- prediction maps -------------------------------------------------------------------------------

b2  <- with(read_csv("output/2nd_order/coef_unscaled_2nd_order.csv", show_col_types = FALSE), setNames(estimate_unscaled, term))
b3  <- with(read_csv("output/3rd_order/coef_unscaled_3rd_order.csv", show_col_types = FALSE), setNames(estimate_unscaled, term))
cap <- read_csv("output/predictions/capping_3rd_order.csv", show_col_types = FALSE)
rnd <- read_csv("output/predictions/predictions_random_points.csv", show_col_types = FALSE)

current <- c(d_edge = "d_edge_19", d_scrub = "d_scrub_19", d_forest = "d_forest_19", d_disturb = "d_disturb_04_13")
layer_file <- function(term, order_dir) {
  p <- sub("^prox_", "d_", term)
  file.path(RASTER_ROOT, order_dir, paste0(if (p %in% names(current)) current[[p]] else p, ".tif"))
}

extent <- st_read(EXTENT, quiet = TRUE) |> st_transform(3310) |> st_union() |> st_make_valid() |>
  vect() |> project(CRS_RASTER)
bb   <- ext(extent)
grid <- rast(ext(floor(bb$xmin / MAP_RES) * MAP_RES, ceiling(bb$xmax / MAP_RES) * MAP_RES,
                 floor(bb$ymin / MAP_RES) * MAP_RES, ceiling(bb$ymax / MAP_RES) * MAP_RES),
             resolution = MAP_RES, crs = CRS_RASTER)

# predictors averaged onto the map grid (proximity = -distance); capping applied before averaging
tmp_file <- file.path(tempdir(), "capped.tif")
predictor_maps <- function(terms, order_dir, cap = NULL) {
  rast(lapply(terms, \(v) {
    s <- if (startsWith(v, "prox_")) -1 else 1
    r <- rast(layer_file(v, order_dir))
    crs(r) <- CRS_RASTER
    if (!is.null(cap)) {
      lim <- sort(s * c(cap$lower[cap$term == v], cap$upper[cap$term == v]))
      r <- clamp(r, lim[1], lim[2], filename = tmp_file, overwrite = TRUE)
    }
    r <- project(r, grid, method = "average") * s
    names(r) <- v
    r
  }))
}
linear_predictor <- function(x, b) {
  out <- x[[names(b)[1]]] * b[[1]]
  for (v in names(b)[-1]) out <- out + x[[v]] * b[[v]]
  out
}

unlink(list.files(map_dir, full.names = TRUE))

x2     <- predictor_maps(names(b2), "order2")
x3     <- predictor_maps(names(b3), "order3", cap)
x3_raw <- predictor_maps(names(b3), "order3")

n_capped <- mask(sum(x3_raw != x3), extent)
names(n_capped) <- "n_capped"
writeRaster(n_capped, file.path(map_dir, "third_order_n_capped.tif"), overwrite = TRUE,
            datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
rm(x3_raw)

maps <- list(second_order = linear_predictor(x2, b2), third_order = linear_predictor(x3, b3))
prank <- function(r, ref) { e <- ecdf(ref); app(r, \(x) e(x)) }
maps$integrated <- prank(maps$second_order, rnd$lp2) * prank(maps$third_order, rnd$lp3)
ref <- list(second_order = rnd$lp2, third_order = rnd$lp3, integrated = rnd$int)

for (nm in names(maps)) {
  r <- mask(maps[[nm]], extent)
  names(r) <- nm
  maps[[nm]] <- r
  writeRaster(r, file.path(map_dir, paste0(nm, ".tif")), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
  e   <- ecdf(ref[[nm]])
  pct <- app(r, \(x) round(100 * e(x)))
  names(pct) <- nm
  writeRaster(pct, file.path(map_dir, paste0(nm, "_percentile.tif")), overwrite = TRUE,
              datatype = "INT1U", gdal = "COMPRESS=DEFLATE")
}
unlink(tmp_file)

# map values against point predictions
xy <- crds(project(vect(as.data.frame(rnd[c("X", "Y")]), geom = c("X", "Y"), crs = "EPSG:3310"), CRS_RASTER))
r_check <- sapply(names(maps), \(nm) cor(extract(maps[[nm]], cellFromXY(maps[[nm]], xy))[[1]], ref[[nm]],
                                          use = "complete.obs"))
print(round(r_check, 2))

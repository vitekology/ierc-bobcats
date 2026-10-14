# Habitat rasters exported from the geodatabases to GeoTIFF (native CRS and grids)
#   order2: 79 m moving-window means (2nd order); order3: unsmoothed (3rd order)

library(terra)
library(sf)
source("analysis/settings.R")

BUFFER_M <- 10000
N_CHECK  <- 2000

set.seed(20261004)
stopifnot(compareVersion(gdal(), "3.7.0") >= 0)
terraOptions(progress = 0)

std_name <- function(x) {
  x <- sub("^mw79_", "", x)
  x <- sub("^dedge_", "d_edge_", x)
  sub("^d_LT_", "d_disturb_", x)
}

jobs <- do.call(rbind, lapply(c("order2", "order3"), function(ord) {
  s <- describe(file.path(GDB_DIR, paste0("analysis_ready_", sub("order", "", ord), "order.gdb")), sds = TRUE)
  data.frame(order = ord, layer = std_name(s$var), sds = s$name)
}))
jobs$out <- file.path(RASTER_ROOT, jobs$order, paste0(jobs$layer, ".tif"))

nlcd <- c(outer(c("d_edge_", "d_scrub_", "d_forest_"), c("11", "13", "16", "19"), paste0))
lt   <- paste0("d_disturb_", c("97_06", "98_07", "99_08", "00_09", "01_10", "02_11", "03_12", "04_13"))
need <- c("d_stream", "d_roads_maj", "d_roads_min", "d_development", "elevation", "vrm", "slope",
          "forest_age", nlcd, lt)
for (ord in c("order2", "order3")) {
  miss <- setdiff(need, jobs$layer[jobs$order == ord])
  if (length(miss) > 0) stop(ord, ": missing layers: ", paste(miss, collapse = ", "))
}

for (d in unique(dirname(jobs$out))) dir.create(d, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(RASTER_ROOT, "_partial\\.tif$", recursive = TRUE, full.names = TRUE))

grids <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
  r <- rast(jobs$sds[i])
  data.frame(order = jobs$order[i], layer = jobs$layer[i], proj = crs(r, proj = TRUE),
             cell_m = round(res(r)[1], 3), origin_x = round(origin(r)[1], 2), origin_y = round(origin(r)[2], 2))
}))
write.csv(grids, file.path(RASTER_ROOT, "source_grids.csv"), row.names = FALSE)
stopifnot(all(grids$proj == CRS_RASTER))

extent   <- st_read(EXTENT, quiet = TRUE) |> st_transform(3310) |> st_union() |> st_make_valid()
buf      <- project(vect(st_buffer(extent, BUFFER_M)), CRS_RASTER)
check_xy <- crds(spatSample(project(vect(extent), CRS_RASTER), N_CHECK, method = "random"))

todo <- jobs[!file.exists(jobs$out), ]
for (i in seq_len(nrow(todo))) {
  src <- rast(todo$sds[i])
  crs(src)   <- CRS_RASTER
  names(src) <- todo$layer[i]
  tmp <- sub("\\.tif$", "_partial.tif", todo$out[i])

  r <- crop(src, buf, snap = "out")
  r <- mask(r, buf, filename = tmp, overwrite = TRUE, datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER"))

  out <- rast(tmp)
  a   <- extract(src, cellFromXY(src, check_xy))[[1]]
  b   <- extract(out, cellFromXY(out, check_xy))[[1]]
  bad <- is.na(a) != is.na(b) | (!is.na(a) & !is.na(b) & abs(a - b) > 1e-6 * pmax(1, abs(a)))
  if (any(bad)) stop("Exported values differ from the source: ", todo$order[i], "/", todo$layer[i])

  rm(r, out); gc()
  if (!file.rename(tmp, todo$out[i])) stop("Could not rename ", tmp)
  tmpFiles(remove = TRUE)
  message(i, "/", nrow(todo), "  ", todo$order[i], "/", todo$layer[i])
}

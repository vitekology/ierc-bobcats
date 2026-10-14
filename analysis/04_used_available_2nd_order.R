# Used (iNaturalist bobcat records) and available points, 2nd-order RSF
# Available points: target-group background (1-km cells with >= K_MIN mammal records)

library(tidyverse)
library(sf)
source("analysis/settings.R")

set.seed(20261001)

inat_dir <- "data/iNaturalist"
out_dir  <- "output/2nd_order"
CELL_M   <- 1000
RATIO    <- 100
dir.create("data/annotated", showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

extent <- st_read(EXTENT, quiet = TRUE) |> st_transform(3310) |> st_union() |> st_make_valid()

read_pts <- function(f) read_csv(file.path(inat_dir, f), show_col_types = FALSE) |>
  st_as_sf(coords = c("X", "Y"), crs = 3310, remove = FALSE)
bobcat <- read_pts("inat_bobcat_filtered.csv") |> mutate(species = "Lynx rufus")
tg     <- bind_rows(select(bobcat, id, species, user_id, year),
                    select(read_pts("inat_targetgroup_filtered.csv"), id, species, user_id, year))

bb    <- st_bbox(extent)
grid  <- st_make_grid(extent, cellsize = CELL_M, offset = floor(bb[c("xmin", "ymin")] / CELL_M) * CELL_M) |>
  st_sf(geometry = _)
cells <- grid[st_within(grid, extent, sparse = FALSE)[, 1], ]
cells$n_records <- lengths(st_intersects(cells, tg))
cells$n_bobcat  <- lengths(st_intersects(cells, bobcat))
st_write(cells, file.path(out_dir, "targetgroup_cells.gpkg"), delete_dsn = TRUE, quiet = TRUE)

domain <- cells |> filter(n_records >= K_MIN) |> st_union() |> st_make_valid()
used   <- st_filter(bobcat, domain)
xy     <- st_coordinates(st_sample(domain, size = RATIO * nrow(used), type = "random", exact = TRUE))

df <- bind_rows(tibble(case = 1L, id = used$id, year = as.integer(used$year), X = used$X, Y = used$Y),
                tibble(case = 0L, id = NA, year = NA_integer_, X = xy[, 1], Y = xy[, 2]))
write_csv(df, "data/annotated/used_available_2nd_order.csv")
print(count(df, case))

png(file.path(out_dir, "availability_map.png"), width = 1000, height = 2000, res = 170)
par(mar = c(1, 1, 1, 1))
plot(extent, col = "grey95", border = "grey40")
plot(domain, col = "steelblue", border = NA, add = TRUE)
points(used$X, used$Y, pch = 16, cex = 0.5, col = "red")
dev.off()

# Habitat values at camera deployments (layers matched to deployment year) and at random points in the extent

library(tidyverse)
library(terra)
library(sf)
source("analysis/settings.R")

set.seed(20261003)
N_RANDOM     <- 10000
CURRENT_YEAR <- 2019

static2 <- c("d_stream", "elevation", "vrm", "slope", "forest_age")
static3 <- c(static2, "d_roads_maj", "d_roads_min", "d_development")

cams <- read_csv("data/camera/camera_deployments.csv", show_col_types = FALSE) |>
  mutate(year = year(start))

extent <- st_read(EXTENT, quiet = TRUE) |> st_transform(3310) |> st_union() |> st_make_valid()
rnd <- st_sample(extent, N_RANDOM, type = "random", exact = TRUE) |>
  st_coordinates() |> as_tibble() |> mutate(pt = row_number(), year = CURRENT_YEAR, .before = 1)

# o2_: smoothed layers (2nd order), o3_: unsmoothed layers (3rd order); distances as proximity
annotate_points <- function(pts) {
  values <- function(order_dir, static, prefix) {
    as_tibble(habitat_values(pts, order_dir, pts$year, static)) |>
      mutate(across(starts_with("d_"), \(x) -x)) |>
      rename_with(\(n) paste0(prefix, sub("^d_", "prox_", n)))
  }
  bind_cols(pts, values("order2", static2, "o2_"), values("order3", static3, "o3_"))
}

write_csv(annotate_points(cams), "data/annotated/predictors_cameras.csv")
write_csv(annotate_points(rnd),  "data/annotated/predictors_random_points.csv")

# Used (GPS fixes) and available points (within each 95% AKDE home range), 3rd-order RSF

library(tidyverse)
library(sf)

set.seed(20261002)
RATIO <- 100

akde <- readRDS("output/homeranges/akde_95.rds")
gps  <- read_csv("data/gps/bobcats_gps_clean.csv", show_col_types = FALSE) |>
  rename(id = individual.local.identifier) |>
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) |>
  st_transform(3310)
ids <- sort(unique(gps$id))
stopifnot(setequal(ids, names(akde)))

df <- map_dfr(ids, \(i) {
  hr    <- st_transform(akde[[i]], 3310) |> st_union() |> st_make_valid()
  used  <- filter(gps, id == i)
  avail <- st_sample(hr, size = RATIO * nrow(used), type = "random", exact = TRUE)
  bind_rows(tibble(id = i, case = 1L, as_tibble(st_coordinates(used))),
            tibble(id = i, case = 0L, as_tibble(st_coordinates(avail))))
})
dir.create("data/annotated", showWarnings = FALSE, recursive = TRUE)
write_csv(df, "data/annotated/used_available_3rd_order.csv.gz")

# available points of the example animal (home-range figure)
df |>
  filter(id == "KR_M03", case == 0) |>
  st_as_sf(coords = c("X", "Y"), crs = 3310) |>
  st_write("data/annotated/available_KR_M03.gpkg", delete_dsn = TRUE, quiet = TRUE)
print(count(df, id, case) |> pivot_wider(names_from = case, values_from = n), n = Inf)

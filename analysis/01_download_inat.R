# iNaturalist research-grade records of bobcats and of the mammal target group

library(tidyverse)
library(httr2)
library(sf)
source("analysis/settings.R")

inat_dir <- "data/iNaturalist"
DOWNLOAD <- FALSE     # FALSE: refilter the latest saved download
D1 <- "2010-01-01"; D2 <- "2020-12-31"
MAX_ACC <- 500

target_species <- c(
  "Canis latrans", "Urocyon cinereoargenteus", "Vulpes vulpes", "Vulpes macrotis",
  "Procyon lotor", "Bassariscus astutus", "Mephitis mephitis", "Spilogale gracilis",
  "Taxidea taxus", "Pekania pennanti", "Martes caurina", "Erethizon dorsatum",
  "Didelphis virginiana", "Puma concolor", "Ursus americanus", "Odocoileus hemionus",
  "Cervus canadensis", "Sus scrofa")

extent <- st_read(EXTENT, quiet = TRUE) |> st_transform(3310) |> st_union() |> st_make_valid()
bb     <- st_bbox(st_transform(extent, 4326))

## ---- download -----------------------------------------------------------------------------

nz <- function(x) if (is.null(x)) NA else x
parse_obs <- function(o) tibble(
  id = o$id, observed_on = nz(o$observed_on), year = nz(o$observed_on_details$year),
  taxon_id = nz(o$taxon$id), scientific_name = nz(o$taxon$name), taxon_rank = nz(o$taxon$rank),
  longitude = nz(o$geojson$coordinates[[1]]), latitude = nz(o$geojson$coordinates[[2]]),
  positional_accuracy = nz(o$positional_accuracy), coordinates_obscured = nz(o$obscured),
  quality_grade = nz(o$quality_grade), captive = nz(o$captive), license = nz(o$license_code),
  user_id = nz(o$user$id), user_login = nz(o$user$login))

inat_pull <- function(taxon_id) {
  out <- list(); last_id <- 0
  repeat {
    resp <- request("https://api.inaturalist.org/v1/observations") |>
      req_url_query(taxon_id = taxon_id, d1 = D1, d2 = D2, quality_grade = "research",
                    swlat = bb[["ymin"]], swlng = bb[["xmin"]],
                    nelat = bb[["ymax"]], nelng = bb[["xmax"]],
                    per_page = 200, order_by = "id", order = "asc", id_above = last_id) |>
      req_user_agent("bobcat-rsf") |>
      req_retry(max_tries = 5, backoff = ~ 10) |>
      req_perform() |>
      resp_body_json()
    if (length(resp$results) == 0) break
    out[[length(out) + 1]] <- map_dfr(resp$results, parse_obs)
    last_id <- max(out[[length(out)]]$id)
    Sys.sleep(1.1)
  }
  bind_rows(out)
}

latest <- function(prefix) {
  f <- sort(list.files(inat_dir, paste0("^", prefix, "_raw_[0-9]{8}\\.csv$"), full.names = TRUE))
  if (length(f) == 0) stop("No saved download for ", prefix)
  read_csv(tail(f, 1), show_col_types = FALSE, guess_max = Inf)
}

if (DOWNLOAD) {
  stamp <- format(Sys.Date(), "%Y%m%d")
  bobcat_raw <- inat_pull(41976)     # Lynx rufus
  write_csv(bobcat_raw, file.path(inat_dir, paste0("inat_bobcat_raw_", stamp, ".csv")))
  mammal_raw <- inat_pull(40151)     # Mammalia
  write_csv(mammal_raw, file.path(inat_dir, paste0("inat_mammal_raw_", stamp, ".csv")))
} else {
  bobcat_raw <- latest("inat_bobcat")
  mammal_raw <- latest("inat_mammal")
}

## ---- filter --------------------------------------------------------------------------------

filter_inat <- function(d) {
  d <- d |>
    filter(!coordinates_obscured,
           is.na(positional_accuracy) | positional_accuracy <= MAX_ACC,
           !is.na(longitude), !is.na(latitude)) |>
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) |>
    st_transform(3310) |>
    st_filter(extent)
  bind_cols(st_drop_geometry(d), as_tibble(st_coordinates(d)))
}

bobcat <- filter_inat(bobcat_raw)
target <- mammal_raw |>
  mutate(species = sub("^(\\S+ \\S+).*$", "\\1", scientific_name)) |>
  filter(species %in% target_species) |>
  filter_inat()

write_csv(bobcat, file.path(inat_dir, "inat_bobcat_filtered.csv"))
write_csv(target, file.path(inat_dir, "inat_targetgroup_filtered.csv"))
message("Bobcat records: ", nrow(bobcat), "; target-group records: ", nrow(target))

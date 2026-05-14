###
# Bobcat Habitat Selection - Used and available point generation (3rd-order RSF)
# Wengert et al. 2026, Ecology and Evolution
#
# Generates systematically-placed available points within each bobcat's 95%
# AKDE home range at 10x the number of used GPS locations, erases available
# points within 30 m of used locations to account for GPS error, then merges
# used (case=1) and available (case=0) into a combined dataset for export.
#
# The exported CSV and shapefile were annotated with habitat raster values
# in ArcGIS Pro 2.9 using the layers described in Table 1 of the manuscript.
# The annotated output is provided directly in data/annotated/ so this script
# does not need to be re-run to reproduce RSF modeling (script 04).
#
# Input:  output/bobcat_AKDE_95.rda      (from script 01)
#         data/gps/bobcats_gps_clean.csv
# Output: data/annotated/bobcats_up_ap.csv   (for ArcGIS; already annotated
#         data/annotated/bobcats_up_ap.shp    version provided in repo)
#
###

library(sf)
library(tidyverse)
library(tictoc)

tic("total")

gps_dir    <- "data/gps"
output_dir <- "output"
annotated_dir <- "data/annotated"
dir.create(annotated_dir, showWarnings = FALSE)

# load AKDE home range polygons
akde <- readRDS(file.path(output_dir, "bobcat_AKDE_95.rda"))

# load used (GPS) points
bobcats.gps <- read.csv(file.path(gps_dir, "bobcats_gps_clean.csv"), as.is = TRUE)
bobcats.gps$individual.local.identifier <- str_replace(
  bobcats.gps$individual.local.identifier, "_", ".")

# sample size per individual
bobcats.n <- bobcats.gps %>%
  group_by(individual.local.identifier) %>%
  summarise(n = n(), .groups = "drop")

varNames <- unique(bobcats.gps$individual.local.identifier)

result_list <- list()

for (i in varNames) {
  j <- match(i, varNames)

  # home range polygon in CA Albers
  temp_hr <- st_transform(akde[[j]], st_crs(3310))

  # systematically sample available points (10x used)
  n_used  <- bobcats.n$n[j]
  temp.ap <- st_sample(temp_hr, size = 10 * n_used, type = "regular")

  # used points
  temp.up <- bobcats.gps %>%
    filter(individual.local.identifier == i) %>%
    st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>%
    st_transform(st_crs(3310))

  # buffer used points by 30 m and erase from available
  temp.up.buff  <- st_union(st_buffer(temp.up, 30))
  temp.ap.ready <- st_difference(temp.ap, temp.up.buff)

  # assemble
  ap_df      <- as.data.frame(st_coordinates(temp.ap.ready))
  ap_df$case <- 0; ap_df$id <- i

  up_df      <- as.data.frame(st_coordinates(temp.up))
  up_df$case <- 1; up_df$id <- i

  result_list[[i]] <- bind_rows(up_df, ap_df)
}

combined <- bind_rows(result_list)

write.csv(combined,
          file.path(annotated_dir, "bobcats_up_ap.csv"),
          row.names = FALSE)

combined.sf <- st_as_sf(combined, coords = c("X", "Y"), crs = st_crs(3310))
st_write(combined.sf,
         file.path(annotated_dir, "bobcats_up_ap.shp"),
         delete_layer = TRUE)

toc()
message("Used/available points exported: n = ", nrow(combined),
        " (", sum(combined$case), " used, ", sum(combined$case == 0), " available)")

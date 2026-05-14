###
# Bobcat Habitat Selection - Used and available point generation (2nd-order RSF)
# Wengert et al. 2026, Ecology and Evolution
#
# Processes iNaturalist bobcat observations (curated from GBIF), filters to
# the study extent and 2010-2020, generates systematically-placed available
# points (10x used) across the study extent, erases available points within
# 79 m of used points (mean reported positional accuracy = 79 m), and exports
# for ArcGIS habitat annotation.
#
# The exported files were annotated with habitat raster values in ArcGIS Pro
# 2.9. The annotated output is provided directly in data/annotated/ so this
# script does not need to be re-run to reproduce RSF modeling (script 06).
#
# Input:  data/iNaturalist/LR.csv, LRF.csv
#         data/iNaturalist/study_extent.shp
# Output: data/annotated/iNat_used_available.csv
#         data/annotated/iNat_used_available.shp
#
###

library(tidyverse)
library(lubridate)
library(sf)

inat_dir      <- "data/iNaturalist"
annotated_dir <- "data/annotated"
dir.create(annotated_dir, showWarnings = FALSE)

## ---- load and filter iNaturalist observations --------------------------------

LR  <- read.csv(file.path(inat_dir, "LR.csv"),  as.is = TRUE)
LRF <- read.csv(file.path(inat_dir, "LRF.csv"), as.is = TRUE)
df  <- rbind(LR, LRF)

# temporal filter: 2010-2020
df$observed_on <- ymd(df$observed_on, tz = "America/Los_Angeles")
df$year        <- as.integer(format(df$observed_on, "%Y"))
df <- subset(df, year > 2009 & year < 2021)

# spatial filter: remove obscured locations and positional accuracy > 500 m
df <- subset(df, coordinates_obscured == "false")
df <- subset(df, positional_accuracy <= 500 | is.na(positional_accuracy))

## ---- spatial operations ------------------------------------------------------

df.sp <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
df.sp <- st_transform(df.sp, st_crs(3310))

# clip to study extent
extent <- st_read(file.path(inat_dir, "study_extent.shp"))
extent <- st_transform(extent, st_crs(3310))
df.sp.extent <- st_intersection(df.sp, extent)

message("iNat observations in study extent: n = ", nrow(df.sp.extent))

# generate available points (10x used, systematic)
inat.ap <- st_sample(extent, size = 10 * nrow(df.sp.extent), type = "regular")

# mean positional accuracy: 79 m; buffer and erase
df.sp.extent.buff <- st_union(st_buffer(df.sp.extent, 79))
inat.ap.ready     <- st_difference(inat.ap, df.sp.extent.buff)

## ---- assemble and export ----------------------------------------------------

ap_df      <- as.data.frame(st_coordinates(inat.ap.ready))
ap_df$year <- NA; ap_df$case <- 0

up_df      <- as.data.frame(st_coordinates(df.sp.extent))
up_df$year <- df.sp.extent$year; up_df$case <- 1

df.ready <- rbind(up_df, ap_df)

write.csv(df.ready,
          file.path(annotated_dir, "iNat_used_available.csv"),
          row.names = FALSE)

df.ready.sf <- st_as_sf(df.ready, coords = c("X", "Y"), crs = st_crs(3310))
st_write(df.ready.sf,
         file.path(annotated_dir, "iNat_used_available.shp"),
         delete_layer = TRUE)

message("2nd-order used/available exported: n = ", nrow(df.ready),
        " (", sum(df.ready$case), " used, ", sum(df.ready$case == 0), " available)")

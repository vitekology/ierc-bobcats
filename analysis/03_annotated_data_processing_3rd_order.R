###
# Bobcat Habitat Selection - Process ArcGIS-annotated data (3rd-order RSF)
# Wengert et al. 2026, Ecology and Evolution
#
# Reads ArcGIS-annotated used/available points and assigns the temporally-
# appropriate version of each dynamic habitat variable to each individual
# based on their tracking period, then exports the formatted RSF input dataset.
#
# Dynamic variables have multiple versions by NLCD/LandTrendr release year:
#   Edge / scrub / forest: NLCD 2011, 2013, 2016, 2019
#   Disturbance:           LandTrendr 5-15 yr windows matched to tracking end
#
# Input:  data/annotated/bobcats_up_ap_annotated.csv
# Output: data/annotated/bobcats_RSF_data.csv
#
###

library(tidyverse)

annotated_dir <- "data/annotated"

df <- read.csv(file.path(annotated_dir, "bobcats_up_ap_annotated.csv"), as.is = TRUE)
df <- df[, -1]  # remove ArcGIS index column

## ---- temporal matching lookup tables ----------------------------------------

# NLCD year suffix per individual (based on tracking period)
nlcd_year <- c(
  Hoopa.F02 = "11", Hoopa.F03 = "11", Hoopa.F04 = "11",
  Hoopa.M01 = "11", Hoopa.M02 = "11", Hoopa.M03 = "11",
  Hoopa.M04 = "11", KR.M01    = "11", KR.M02    = "11",
  SN.M01    = "11",
  KR.F02    = "13", KR.M03    = "13", KR.M04    = "13",
  KR.M05    = "13", Stirling.M01 = "13",
  Hoopa.M07 = "16", KR.M06    = "16",
  Hoopa.M08 = "19"
)

# LandTrendr layer per individual (5-15 yr disturbance window prior to tracking end)
lt_layer <- c(
  Hoopa.F02    = "lt97_06", Hoopa.M01    = "lt97_06",
  Hoopa.F03    = "lt98_07", Hoopa.F04    = "lt98_07",
  Hoopa.M02    = "lt98_07", Hoopa.M03    = "lt98_07",
  Hoopa.M04    = "lt98_07", KR.M01       = "lt98_07",
  KR.M02       = "lt98_07", SN.M01       = "lt98_07",
  KR.F02       = "lt99_08",
  KR.M03       = "lt00_09", KR.M04       = "lt00_09",
  Stirling.M01 = "lt00_09", KR.M05       = "lt01_10",
  KR.M06       = "lt02_11", Hoopa.M07    = "lt03_12",
  Hoopa.M08    = "lt04_13"
)

## ---- assign temporally-matched values per individual ------------------------

assign_by_id <- function(df, prefix, lookup) {
  result <- rep(NA_real_, nrow(df))
  for (id in names(lookup)) {
    col <- paste0(prefix, lookup[[id]])
    idx <- df$id == id
    if (col %in% names(df)) result[idx] <- as.numeric(df[[col]][idx])
  }
  result
}

df$d_edge   <- assign_by_id(df, "edge_",     nlcd_year)
df$d_scrub  <- assign_by_id(df, "d_scrub_",  nlcd_year)
df$d_forest <- assign_by_id(df, "d_forest_", nlcd_year)
df$d_disturb <- assign_by_id(df, "",         lt_layer)

stopifnot(sum(is.na(df$d_edge))    == 0)
stopifnot(sum(is.na(df$d_scrub))   == 0)
stopifnot(sum(is.na(df$d_forest))  == 0)
stopifnot(sum(is.na(df$d_disturb)) == 0)

## ---- select final columns and export ----------------------------------------

df.final <- df %>%
  select(case, id, d_edge, d_stream, d_roads_maj, d_roads_min,
         elevation, vrm, slope, d_development, forest_age,
         d_disturb, d_scrub, d_forest)

write.csv(df.final,
          file.path(annotated_dir, "bobcats_RSF_data.csv"),
          row.names = FALSE)

message("3rd-order RSF data ready: n = ", nrow(df.final),
        " (", sum(df.final$case), " used, ", sum(df.final$case == 0), " available)")

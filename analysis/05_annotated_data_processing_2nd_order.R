###
# Bobcat Habitat Selection - Process ArcGIS-annotated data (2nd-order RSF)
# Wengert et al. 2026, Ecology and Evolution
#
# Reads ArcGIS-annotated iNaturalist used/available points and assigns the
# temporally-appropriate version of dynamic habitat variables based on
# observation year, then exports the formatted 2nd-order RSF input dataset.
#
# Available points (year = NA) are assigned to the 2019 (most recent) layers,
# reflecting the increasing use of iNaturalist in recent years and aligning
# with the majority of observations.
#
# Dynamic variable matching:
#   Edge / scrub / forest: year < 2013  -> NLCD 2011
#                          2013-2015    -> NLCD 2013
#                          2016-2018    -> NLCD 2016
#                          >= 2019      -> NLCD 2019
#   Disturbance (LandTrendr): 5-15 yr window prior to observation year
#
# Input:  data/annotated/iNat_up_ap_annotated.csv
# Output: data/annotated/iNat_RSF_data.csv
#
###

library(tidyverse)

annotated_dir <- "data/annotated"

df <- read.csv(file.path(annotated_dir, "iNat_up_ap_annotated.csv"), as.is = TRUE)
df <- df[, -1]  # remove ArcGIS index column

# available points assigned to 2019 layers
df$year <- ifelse(df$case == 0, 2019, as.integer(df$year))

## ---- temporally-matched variable assignment ----------------------------------

# NLCD suffix by year
nlcd_sfx <- case_when(
  df$year < 2013               ~ "11",
  df$year >= 2013 & df$year < 2016 ~ "13",
  df$year >= 2016 & df$year < 2019 ~ "16",
  df$year >= 2019              ~ "19"
)

# vectorised lookup across rows
lookup_nlcd <- function(df, prefix, sfx) {
  mapply(function(s, i) {
    col <- paste0(prefix, s)
    if (col %in% names(df)) as.numeric(df[[col]][i]) else NA_real_
  }, sfx, seq_len(nrow(df)))
}

df$d_edge   <- lookup_nlcd(df, "d_edge_",   nlcd_sfx)
df$d_scrub  <- lookup_nlcd(df, "d_scrub_",  nlcd_sfx)
df$d_forest <- lookup_nlcd(df, "d_forest_", nlcd_sfx)

# LandTrendr disturbance (5-15 yr prior window)
df$d_disturb <- case_when(
  df$year >= 2018 ~ df$d_LT_04_13,
  df$year == 2017 ~ df$d_LT_03_12,
  df$year == 2016 ~ df$d_LT_02_11,
  df$year == 2015 ~ df$d_LT_01_10,
  df$year == 2014 ~ df$d_LT_00_09,
  df$year == 2013 ~ df$d_LT_99_08,
  df$year == 2012 ~ df$d_LT_98_07,
  df$year <= 2011 ~ df$d_LT_97_06,
  TRUE            ~ NA_real_
)

# remove ArcGIS NoData rows and coerce to numeric
df <- subset(df, d_edge != "None" & elevation != "None" &
               vrm != "None" & slope != "None" & forest_age != "None")
for (v in c("d_edge", "elevation", "vrm", "slope", "forest_age")) {
  df[[v]] <- as.numeric(df[[v]])
}

stopifnot(sum(is.na(df$d_edge))    == 0)
stopifnot(sum(is.na(df$d_scrub))   == 0)
stopifnot(sum(is.na(df$d_forest))  == 0)
stopifnot(sum(is.na(df$d_disturb)) == 0)

## ---- select final columns and export ----------------------------------------

df.final <- df %>%
  select(case, year, d_edge, d_stream, elevation, vrm, slope,
         forest_age, d_disturb, d_scrub, d_forest)

write.csv(df.final,
          file.path(annotated_dir, "iNat_RSF_data.csv"),
          row.names = FALSE)

message("2nd-order RSF data ready: n = ", nrow(df.final),
        " (", sum(df.final$case), " used, ", sum(df.final$case == 0), " available)")

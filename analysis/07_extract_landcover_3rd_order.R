# Habitat values at 3rd-order points (unsmoothed layers)
# Time-varying layers matched to each animal's tracking period

library(tidyverse)
library(terra)
source("analysis/settings.R")

out_dir <- "output/3rd_order"
W_AVAIL <- 1000
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

static  <- c("d_stream", "d_roads_maj", "d_roads_min", "elevation", "vrm", "slope",
             "d_development", "forest_age")
dynamic <- c("d_edge", "d_scrub", "d_forest", "d_disturb")
preds   <- c(dynamic, static)

nlcd_year <- c(
  Hoopa_F02 = "11", Hoopa_F03 = "11", Hoopa_F04 = "11", Hoopa_M01 = "11", Hoopa_M02 = "11",
  Hoopa_M03 = "11", Hoopa_M04 = "11", KR_M01 = "11", KR_M02 = "11", SN_M01 = "11",
  KR_F02 = "13", KR_M03 = "13", KR_M04 = "13", KR_M05 = "13", Stirling_M01 = "13",
  Hoopa_M07 = "16", KR_M06 = "16", Hoopa_M08 = "19")

lt_window <- c(
  Hoopa_F02 = "97_06", Hoopa_M01 = "97_06",
  Hoopa_F03 = "98_07", Hoopa_F04 = "98_07", Hoopa_M02 = "98_07", Hoopa_M03 = "98_07",
  Hoopa_M04 = "98_07", KR_M01 = "98_07", KR_M02 = "98_07", SN_M01 = "98_07",
  KR_F02 = "99_08", KR_M03 = "00_09", KR_M04 = "00_09", Stirling_M01 = "00_09",
  KR_M05 = "01_10", KR_M06 = "02_11", Hoopa_M07 = "03_12", Hoopa_M08 = "04_13")

df <- read_csv("data/annotated/used_available_3rd_order.csv.gz", show_col_types = FALSE)
stopifnot(all(unique(df$id) %in% names(nlcd_year)), all(unique(df$id) %in% names(lt_window)))

out <- map_dfr(unique(df$id), \(i) {
  d   <- filter(df, id == i)
  lyr <- c(static, paste0(c("d_edge_", "d_scrub_", "d_forest_"), nlcd_year[[i]]),
           paste0("d_disturb_", lt_window[[i]]))
  ex  <- raster_values(d, "order3", lyr)
  colnames(ex) <- c(static, dynamic)
  bind_cols(d, as_tibble(ex)) |> mutate(nlcd = nlcd_year[[i]], lt = lt_window[[i]])
})

dist_vars <- grep("^d_", preds, value = TRUE)
out <- out |>
  filter(complete.cases(out[preds])) |>
  mutate(weight = if_else(case == 1, 1, W_AVAIL),
         across(all_of(dist_vars), \(x) -x)) |>
  rename_with(\(n) sub("^d_", "prox_", n), all_of(dist_vars)) |>
  select(id, case, weight, nlcd, lt, X, Y, all_of(sub("^d_", "prox_", preds)))
preds <- sub("^d_", "prox_", preds)

write_csv(out, "data/annotated/RSF_data_3rd_order.csv.gz")
print(count(out, case))

cors <- as.data.frame(as.table(cor(out[preds])), responseName = "r") |>
  rename(var1 = Var1, var2 = Var2) |>
  filter(as.character(var1) < as.character(var2)) |>
  mutate(r = round(r, 2)) |>
  arrange(desc(abs(r)))
write_csv(cors, file.path(out_dir, "predictor_correlations.csv"))

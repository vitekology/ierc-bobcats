# Habitat values at 2nd-order points (79 m moving-window layers)
# Time-varying layers matched to the record year; available points get years drawn from the used records

library(tidyverse)
library(terra)
source("analysis/settings.R")

set.seed(20261005)

out_dir <- "output/2nd_order"
W_AVAIL <- 1000
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

static  <- c("d_stream", "elevation", "vrm", "slope", "forest_age")
dynamic <- c("d_edge", "d_scrub", "d_forest", "d_disturb")
preds   <- c(dynamic, static)

acc <- read_csv("data/iNaturalist/inat_bobcat_filtered.csv", show_col_types = FALSE)$positional_accuracy
message(sprintf("Positional accuracy: mean %.1f m (n = %d of %d)", mean(acc, na.rm = TRUE), sum(!is.na(acc)), length(acc)))

df <- read_csv("data/annotated/used_available_2nd_order.csv", show_col_types = FALSE) |>
  mutate(year = as.integer(year),
         layer_year = if_else(case == 1, year, sample(year[case == 1], n(), replace = TRUE)))

df <- bind_cols(df, as_tibble(habitat_values(df, "order2", df$layer_year, static)))

dist_vars <- grep("^d_", preds, value = TRUE)
df <- df |>
  filter(complete.cases(df[preds])) |>
  mutate(weight = if_else(case == 1, 1, W_AVAIL),
         across(all_of(dist_vars), \(x) -x)) |>
  rename_with(\(n) sub("^d_", "prox_", n), all_of(dist_vars)) |>
  select(case, weight, id, year, layer_year, X, Y, all_of(sub("^d_", "prox_", preds)))
preds <- sub("^d_", "prox_", preds)

write_csv(df, "data/annotated/RSF_data_2nd_order.csv")
print(count(df, case))

cors <- as.data.frame(as.table(cor(df[preds])), responseName = "r") |>
  rename(var1 = Var1, var2 = Var2) |>
  filter(as.character(var1) < as.character(var2)) |>
  mutate(r = round(r, 2)) |>
  arrange(desc(abs(r)))
write_csv(cors, file.path(out_dir, "predictor_correlations.csv"))

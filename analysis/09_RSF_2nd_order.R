# 2nd-order RSF: weighted logistic regression of iNaturalist records against the target-group background

library(tidyverse)
library(glmmTMB)
library(broom.mixed)
source("analysis/settings.R")

out_dir    <- "output/2nd_order"
candidates <- c("prox_edge", "prox_forest", "prox_scrub", "prox_disturb", "prox_stream",
                "elevation", "slope", "vrm", "forest_age")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv("data/annotated/RSF_data_2nd_order.csv", show_col_types = FALSE)

sel   <- select_predictors(df, candidates)
preds <- sel$keep
print(sel$pairs)
write_csv(sel$univariate |> mutate(kept = term %in% preds), file.path(out_dir, "predictor_selection_2nd_order.csv"))

sds <- map_dbl(df[preds], sd)
ds  <- df |> mutate(across(all_of(preds), \(x) as.numeric(scale(x))))

m <- glmmTMB(reformulate(preds, response = "case"), family = binomial(), data = ds, weights = weight)
if (m$fit$convergence != 0 || !m$sdr$pdHess) warning("Convergence problem")
saveRDS(m, file.path(out_dir, "RSF_2nd_order.rds"))

coefs <- tidy(m, effects = "fixed", conf.int = TRUE) |>
  filter(term != "(Intercept)") |>
  transmute(n_used = sum(df$case == 1), n_avail = sum(df$case == 0),
            term, estimate, std.error, p.value, conf.low, conf.high,
            estimate_unscaled = estimate / unname(sds[term]))
write_csv(coefs |> select(-estimate_unscaled) |> mutate(across(estimate:conf.high, \(x) signif(x, 3))),
          file.path(out_dir, "coefficients_2nd_order.csv"))
write_csv(select(coefs, term, estimate_unscaled), file.path(out_dir, "coef_unscaled_2nd_order.csv"))
print(select(coefs, term, estimate, conf.low, conf.high, p.value), n = Inf)

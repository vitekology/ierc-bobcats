# 3rd-order RSF: weighted logistic mixed model with random slopes per animal (Muff et al. 2020)

library(tidyverse)
library(glmmTMB)
library(broom.mixed)
source("analysis/settings.R")

out_dir   <- "output/3rd_order"
N_THREADS <- 8
candidates <- c("prox_edge", "prox_forest", "prox_scrub", "prox_disturb", "prox_stream",
                "prox_roads_maj", "prox_roads_min", "prox_development",
                "elevation", "slope", "vrm", "forest_age")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_csv("data/annotated/RSF_data_3rd_order.csv.gz", show_col_types = FALSE)

sel   <- select_predictors(df, candidates)
preds <- sel$keep
print(sel$pairs)
write_csv(sel$univariate |> mutate(kept = term %in% preds), file.path(out_dir, "predictor_selection_3rd_order.csv"))

vif <- tibble(variable = preds, vif = round(diag(solve(cor(df[preds]))), 2)) |> arrange(desc(vif))
write_csv(vif, file.path(out_dir, "vif_3rd_order.csv"))

sds <- map_dbl(df[preds], sd)
df  <- df |> mutate(across(all_of(preds), \(x) as.numeric(scale(x))), id = factor(id))

f <- as.formula(paste("case ~", paste(preds, collapse = " + "),
                      "+ (1 | id) +", paste0("(0 + ", preds, " | id)", collapse = " + ")))
n_re <- 1 + length(preds)

# random intercept SD fixed at 1000
m <- glmmTMB(f, family = binomial(), data = df, weights = weight,
             map   = list(theta = factor(c(NA, seq_len(n_re - 1)))),
             start = list(theta = c(log(1e3), rep(0, n_re - 1))),
             control = glmmTMBControl(parallel = N_THREADS))
if (m$fit$convergence != 0 || !m$sdr$pdHess) warning("Convergence problem")
saveRDS(m, file.path(out_dir, "RSF_3rd_order.rds"))

coefs <- tidy(m, effects = "fixed", conf.int = TRUE) |>
  filter(term != "(Intercept)") |>
  transmute(term, estimate, std.error, p.value, conf.low, conf.high,
            estimate_unscaled = estimate / unname(sds[term]))
write_csv(coefs |> select(-estimate_unscaled) |> mutate(across(where(is.double), \(x) signif(x, 3))),
          file.path(out_dir, "coefficients_3rd_order.csv"))
write_csv(select(coefs, term, estimate_unscaled), file.path(out_dir, "coef_unscaled_3rd_order.csv"))
print(coefs, n = Inf)

slope_sd <- unlist(lapply(VarCorr(m)$cond, \(x) attr(x, "stddev")))
tibble(term = sub(".*\\.", "", names(slope_sd)), sd = signif(slope_sd, 3)) |>
  filter(term %in% preds) |>
  write_csv(file.path(out_dir, "random_slope_sd_3rd_order.csv"))

coef(m)$cond$id |> rownames_to_column("id") |> select(id, all_of(preds)) |>
  write_csv(file.path(out_dir, "individual_coefs_3rd_order.csv"))

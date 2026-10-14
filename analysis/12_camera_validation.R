# Camera validation: detection ~ prediction + study + log(effort) + (1 | site), binomial cloglog,
# offset log(effort) (sdmTMB). Structure compared on the model without a prediction surface.

library(tidyverse)
library(sdmTMB)
library(DHARMa)

out_dir <- "output/validation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261005)

EFFORT  <- "estimated"   # "offset" or "estimated" (exponent of effort)
STUDY   <- "fixed"       # "fixed" or "random"
SPATIAL <- "off"         # Matern spatial field, "on" or "off"
MESH_CUTOFF_KM <- 5
N_SIM   <- 250
MIN_DET <- 20

surfaces <- c(`2nd order` = "lp2", `3rd order` = "lp3", `Integrated RSF` = "int")

cams <- read_csv("output/predictions/predictions_cameras.csv", show_col_types = FALSE) |>
  mutate(study = factor(study), site = factor(paste(study, site)),
         loc = factor(paste(round(X), round(Y))),
         log_effort = log(effort), X_km = X / 1000, Y_km = Y / 1000) |>
  as.data.frame()
rnd <- read_csv("output/predictions/predictions_random_points.csv", show_col_types = FALSE)

message(nrow(cams), " deployments, ", sum(cams$detected), " detections, ", nlevels(cams$site), " sites, ",
        nlevels(cams$study), " studies")

fit_val <- function(d, terms = NULL, effort = EFFORT, study = STUDY, spatial = SPATIAL) {
  d   <- droplevels(d)
  rhs <- c(terms,
           if (effort == "estimated") "log_effort",
           if (nlevels(d$study) > 1) switch(study, fixed = "study", random = "(1 | study)"),
           "(1 | site)")
  mesh <- make_mesh(d, xy_cols = c("X_km", "Y_km"), cutoff = MESH_CUTOFF_KM)
  sdmTMB(reformulate(rhs, response = "detected"), data = d, mesh = mesh,
         family = binomial(link = "cloglog"), spatial = spatial, offset = "log_effort", silent = TRUE)
}
converged <- function(m) {
  s <- tryCatch(sanity(m, silent = TRUE), error = \(e) NULL)
  isTRUE(if (is.list(s)) s$all_ok else s)
}
fixed <- function(m, v) {
  t <- filter(tidy(m, effects = "fixed", conf.int = TRUE), term == v)
  c(est = t$estimate, se = t$std.error, lo = t$conf.low, hi = t$conf.high,
    p = 2 * pnorm(-abs(t$estimate / t$std.error)))
}
ran_par <- function(m, v) {
  t <- tidy(m, effects = "ran_pars")
  if (v %in% t$term) t$estimate[t$term == v][1] else NA_real_
}
# population-level DHARMa residuals (new site effects)
dharma <- function(m) {
  sims <- simulate(m, nsim = N_SIM, type = "mle-mvn", re_form = NA)
  pred <- predict(m, re_form = NA, re_form_iid = NA, type = "response")$est
  createDHARMa(simulatedResponse = sims, observedResponse = m$data$detected,
               fittedPredictedResponse = pred, integerResponse = TRUE)
}
auc <- function(score, y) {
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---- model structure (no prediction surface) ---------------------------------------------------

structures <- expand_grid(effort = c("offset", "estimated"), study = c("fixed", "random"), spatial = c("off", "on"))
nulls <- pmap(structures, \(effort, study, spatial) fit_val(cams, NULL, effort, study, spatial))
null_tab <- structures |>
  mutate(AIC = map_dbl(nulls, AIC),
         effort_exponent = map2_dbl(nulls, effort, \(m, e) if (e == "estimated") 1 + fixed(m, "log_effort")[["est"]] else 1),
         range_km = map_dbl(nulls, \(m) ran_par(m, "range")),
         sigma_spatial = map_dbl(nulls, \(m) ran_par(m, "sigma_O")),
         converged = map_lgl(nulls, converged)) |>
  mutate(dAIC = AIC - min(AIC), across(where(is.double), \(x) signif(x, 4))) |>
  arrange(AIC)
write_csv(null_tab, file.path(out_dir, "null_structures.csv"))
print(null_tab, width = Inf)

m0   <- nulls[[which(structures$effort == EFFORT & structures$study == STUDY & structures$spatial == SPATIAL)]]
res0 <- dharma(m0)
png(file.path(out_dir, "dharma_null.png"), width = 2000, height = 1000, res = 200)
plot(res0)
dev.off()

loc_xy <- cams |> group_by(loc) |> summarise(X_km = mean(X_km), Y_km = mean(Y_km), .groups = "drop")
loc_xy <- loc_xy[match(levels(cams$loc), as.character(loc_xy$loc)), ]
mt <- testSpatialAutocorrelation(recalculateResiduals(res0, group = cams$loc),
                                 x = loc_xy$X_km, y = loc_xy$Y_km, plot = FALSE)
message(sprintf("Moran's I of residuals: %.4f, P = %.3f", mt$statistic[["observed"]], mt$p.value))

## ---- prediction surfaces -------------------------------------------------------------------------

models <- list(null = m0)
tab <- imap_dfr(surfaces, \(v, nm) {
  d <- cams
  d$x <- as.numeric(scale(d[[v]]))
  d$x_raw <- d[[v]]
  m <- fit_val(d, "x")
  models[[nm]] <<- m
  s <- fixed(m, "x")
  r <- dharma(m)
  tibble(surface = nm, slope_per_sd = s[["est"]], se = s[["se"]], ci_low = s[["lo"]], ci_high = s[["hi"]],
         z = s[["est"]] / s[["se"]], p = s[["p"]],
         AUC = auc(d[[v]], d$detected), dAIC_vs_null = AIC(m) - AIC(m0),
         dharma_p_uniformity = testUniformity(r, plot = FALSE)$p.value,
         dharma_p_dispersion = testDispersion(r, plot = FALSE)$p.value,
         converged = converged(m))
}) |> mutate(across(where(is.double), \(x) signif(x, 3)))
write_csv(tab, file.path(out_dir, "validation.csv"))
saveRDS(models, file.path(out_dir, "models_validation.rds"))
print(tab, width = Inf)

big <- cams |> filter(detected == 1) |> count(study) |> filter(n >= MIN_DET) |> pull(study) |> as.character()
by_study <- map_dfr(big, \(st) imap_dfr(surfaces, \(v, nm) {
  d <- filter(cams, study == st)
  d$x <- as.numeric(scale(d[[v]]))
  m <- tryCatch(fit_val(d, "x"), error = \(e) NULL)
  if (is.null(m)) return(tibble(study = st, surface = nm, converged = FALSE))
  s <- fixed(m, "x")
  tibble(study = st, surface = nm, deployments = nrow(d), detections = sum(d$detected),
         slope_per_sd = s[["est"]], ci_low = s[["lo"]], ci_high = s[["hi"]], p = s[["p"]],
         converged = converged(m))
})) |> mutate(across(where(is.double), \(x) signif(x, 3)))
write_csv(by_study, file.path(out_dir, "validation_by_study.csv"))

coverage <- imap_dfr(surfaces, \(v, nm) {
  pct <- ecdf(rnd[[v]])(cams[[v]])
  tibble(surface = nm, median_percentile = 100 * median(pct),
         pct_in_top20 = 100 * mean(pct >= 0.8), pct_in_bottom20 = 100 * mean(pct <= 0.2))
}) |> mutate(across(where(is.double), \(x) round(x)))
write_csv(coverage, file.path(out_dir, "camera_coverage.csv"))

###
# Bobcat Habitat Selection - 3rd-order RSF (GLMM)
# Wengert et al. 2026, Ecology and Evolution
#
# Fits a resource selection function for bobcat 3rd-order (home range scale)
# habitat selection. Follows Muff et al. (2020): a GLMM with binomial error
# where the random intercept variance is fixed at 10^6 (effectively a fixed
# intercept) and all predictor slopes vary freely by individual.
#
# Two model versions:
#   model.scaled: standardized predictors -- used for coefficient reporting
#                 (Table S2, Figure S2) and profile-likelihood CIs
#   model.raw:    unstandardized predictors -- used for spatial prediction
#                 and mapping (Figure 4)
#
# Available locations are down-weighted (w = 1000) relative to used locations
# (w = 1) per the Muff et al. approach.
#
# Note: d_edge and d_forest approached r = 0.70 in the correlation matrix
# but were retained as both were a priori hypothesized (see manuscript).
#
# Input:  data/annotated/bobcats_RSF_data.csv
# Output: output/RSF_3order_scaled.RDS
#         output/RSF_3order_raw.RDS
#         output/RSF_3order_summary_scaled.doc
#         output/corr_matrix_3order.tiff
#
# References:
#   Muff et al. 2020, J Anim Ecol 89:80-92
#   Fieberg et al. 2021, J Anim Ecol 90:1027-1043
#
###

library(tidyverse)
library(glmmTMB)
library(psych)
library(performance)
library(tictoc)
library(sjPlot)

annotated_dir <- "data/annotated"
output_dir    <- "output"
dir.create(output_dir, showWarnings = FALSE)

df <- read.csv(file.path(annotated_dir, "bobcats_RSF_data.csv"), as.is = TRUE)

df <- df %>% rename(d_rds_maj = d_roads_maj,
                    d_rds_min = d_roads_min,
                    d_develop = d_development,
                    for_age   = forest_age)

# center and scale predictors
df <- df %>% mutate(
  s.d_edge    = scale(d_edge),    s.d_stream  = scale(d_stream),
  s.d_rds_maj = scale(d_rds_maj), s.d_rds_min = scale(d_rds_min),
  s.elevation = scale(elevation), s.vrm       = scale(vrm),
  s.slope     = scale(slope),     s.d_develop = scale(d_develop),
  s.for_age   = scale(for_age),   s.d_disturb = scale(d_disturb),
  s.d_scrub   = scale(d_scrub),   s.d_forest  = scale(d_forest)
)

# weights: available = 1000, used = 1
df$w <- ifelse(df$case == 1, 1, 1e3)

## ---- correlation check -------------------------------------------------------

tiff(file.path(output_dir, "corr_matrix_3order.tiff"),
     units = "px", width = 5000, height = 5000, res = 300, compression = "lzw")
pairs.panels(df[, c("d_edge","d_stream","d_rds_maj","d_rds_min","elevation",
                    "vrm","slope","d_develop","for_age","d_disturb","d_scrub","d_forest")],
             method = "pearson", hist.col = "#00AFBB", density = TRUE, ellipses = TRUE)
dev.off()

## ---- Muff et al. fixed-intercept GLMM helper --------------------------------

fit_muff <- function(formula_str, df) {
  tmp <- glmmTMB(as.formula(formula_str),
                 family = binomial(), data = df, weights = w, doFit = FALSE)
  tmp$parameters$theta[1] <- log(1e3)
  n <- length(tmp$parameters$theta)
  tmp$mapArg <- list(theta = factor(c(NA, seq_len(n - 1))))
  glmmTMB:::fitTMB(tmp)
}

## ---- scaled model -----------------------------------------------------------

tic("3rd-order RSF: scaled")

scaled_f <- "case ~
  s.d_edge + s.d_stream + s.d_rds_maj + s.d_rds_min +
  s.elevation + s.vrm + s.slope + s.d_develop +
  s.for_age + s.d_disturb + s.d_scrub + s.d_forest +
  (1|id) +
  (0+s.d_edge|id) + (0+s.d_stream|id) + (0+s.d_rds_maj|id) + (0+s.d_rds_min|id) +
  (0+s.elevation|id) + (0+s.vrm|id) + (0+s.slope|id) + (0+s.d_develop|id) +
  (0+s.for_age|id) + (0+s.d_disturb|id) + (0+s.d_scrub|id) + (0+s.d_forest|id)"

model.scaled <- fit_muff(scaled_f, df)
summary(model.scaled)
model.scaled$sdr$pdHess  # should be TRUE
diagnose(model.scaled)
toc()

saveRDS(model.scaled, file.path(output_dir, "RSF_3order_scaled.RDS"))
tab_model(model.scaled, show.intercept = FALSE, digits = 2,
          transform = NULL, show.stat = TRUE,
          file = file.path(output_dir, "RSF_3order_summary_scaled.doc"))

pdf(file.path(output_dir, "RSF_3order_coefplot.pdf"), width = 5, height = 5)
plot_model(model.scaled, transform = NULL, sort.est = TRUE)
dev.off()

## ---- raw (unstandardized) model for spatial prediction ----------------------

tic("3rd-order RSF: raw")

raw_f <- "case ~
  d_edge + d_stream + d_rds_maj + d_rds_min +
  elevation + vrm + slope + d_develop +
  for_age + d_disturb + d_scrub + d_forest +
  (1|id) +
  (0+d_edge|id) + (0+d_stream|id) + (0+d_rds_maj|id) + (0+d_rds_min|id) +
  (0+elevation|id) + (0+vrm|id) + (0+slope|id) + (0+d_develop|id) +
  (0+for_age|id) + (0+d_disturb|id) + (0+d_scrub|id) + (0+d_forest|id)"

model.raw <- fit_muff(raw_f, df)
summary(model.raw)
model.raw$sdr$pdHess
diagnose(model.raw)
toc()

saveRDS(model.raw, file.path(output_dir, "RSF_3order_raw.RDS"))

message("3rd-order RSF complete.")

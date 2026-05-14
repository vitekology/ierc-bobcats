###
# Bobcat Habitat Selection - 2nd-order RSF (GLM)
# Wengert et al. 2026, Ecology and Evolution
#
# Fits a resource selection function for bobcat 2nd-order (regional scale)
# habitat selection using iNaturalist citizen science observations as used
# locations. A GLM with binomial error is used (no random effects, as there
# is no individual-level structure in the iNaturalist data).
#
# Roads and development predictors are excluded from this analysis due to
# expected human accessibility bias in citizen science data -- iNaturalist
# observations are more likely near roads and developed areas regardless
# of bobcat habitat preference (see manuscript).
#
# d_forest is also excluded: highly correlated with d_edge in this dataset.
#
# Two model versions:
#   model.scaled:       standardized predictors -- coefficient reporting
#                       (Table S2, Figure S1)
#   model.scaled.unst:  unstandardized -- spatial prediction and mapping
#
# Input:  data/annotated/iNat_RSF_data.csv
# Output: output/RSF_2order_scaled.RDS
#         output/RSF_2order_unst.RDS
#         output/RSF_2order_summary.doc
#
###

library(tidyverse)
library(glmmTMB)
library(sjPlot)

annotated_dir <- "data/annotated"
output_dir    <- "output"
dir.create(output_dir, showWarnings = FALSE)

df <- read.csv(file.path(annotated_dir, "iNat_RSF_data.csv"), as.is = TRUE)

df <- df %>% rename(for_age = forest_age)

# center and scale predictors
# (d_forest, roads, development excluded -- see header)
df <- df %>% mutate(
  s.d_edge    = scale(d_edge),
  s.d_stream  = scale(d_stream),
  s.elevation = scale(elevation),
  s.vrm       = scale(vrm),
  s.slope     = scale(slope),
  s.for_age   = scale(for_age),
  s.d_disturb = scale(d_disturb),
  s.d_scrub   = scale(d_scrub)
)

# weights: available = 1000, used = 1
df$w <- ifelse(df$case == 1, 1, 1e3)

## ---- scaled model (for coefficient reporting) --------------------------------

model.scaled <- glmmTMB(case ~ s.d_edge + s.d_stream + s.elevation + s.vrm +
                           s.slope + s.for_age + s.d_disturb + s.d_scrub,
                         family = binomial(), data = df, weights = w)
summary(model.scaled)
saveRDS(model.scaled, file.path(output_dir, "RSF_2order_scaled.RDS"))

tab_model(model.scaled, transform = NULL, digits = 2, show.intercept = FALSE,
          show.stat = TRUE,
          file = file.path(output_dir, "RSF_2order_summary.doc"))

pdf(file.path(output_dir, "RSF_2order_coefplot.pdf"), width = 5, height = 5)
plot_model(model.scaled, transform = NULL, sort.est = TRUE)
dev.off()

## ---- unstandardized model (for spatial prediction) --------------------------

model.unst <- glmmTMB(case ~ d_edge + d_stream + elevation + vrm +
                        slope + for_age + d_disturb + d_scrub,
                      family = binomial(), data = df, weights = w)
summary(model.unst)
saveRDS(model.unst, file.path(output_dir, "RSF_2order_unst.RDS"))

message("2nd-order RSF complete.")

###
# Bobcat Habitat Selection - RSF coefficient back-transformation and
# 2nd/3rd-order agreement analysis
# Wengert et al. 2026, Ecology and Evolution
#
# Part 1: Back-transforms standardized 3rd-order RSF coefficients to the
# original scale for use in spatial prediction. Scaled coefficients are
# divided by their predictor standard deviations.
#
# Profile-likelihood CIs (via broom.mixed::tidy with conf.method = "profile")
# take ~24 hours to compute and are provided pre-computed as output/coef_st.csv.
# Uncomment the relevant block to recompute from scratch if needed.
#
# Part 2: Tests correlation between 2nd- and 3rd-order RSF relative selection
# probability at 1,000 randomly placed points within the study extent.
# Produces Figure 3.
#
# Input:  data/annotated/bobcats_RSF_data.csv
#         output/RSF_3order_scaled.RDS
#         output/coef_st.csv              (pre-computed; provided in repo)
#         data/annotated/rpts_1000.csv    (random points with RSF values extracted)
# Output: output/coef_unst.csv
#         output/RSF_agreement_plot.pdf   (Figure 3)
#
###

library(tidyverse)
library(glmmTMB)
library(purrr)
library(broom.mixed)
library(DHARMa)
library(ggeffects)

annotated_dir <- "data/annotated"
output_dir    <- "output"


## ---- PART 1: BACK-TRANSFORMATION --------------------------------------------

df <- read.csv(file.path(annotated_dir, "bobcats_RSF_data.csv"), as.is = TRUE)
df <- df %>% rename(d_rds_maj = d_roads_maj, d_rds_min = d_roads_min,
                    d_develop = d_development, for_age = forest_age)

model <- readRDS(file.path(output_dir, "RSF_3order_scaled.RDS"))

# load pre-computed profile CIs
# To recompute (takes ~24 hrs):
# coef_st <- tidy(model, effects = "fixed", conf.int = TRUE, conf.method = "profile")
# write.csv(coef_st, file.path(output_dir, "coef_st.csv"), row.names = FALSE)
coef_st <- read.csv(file.path(output_dir, "coef_st.csv"), as.is = TRUE)

# predictor standard deviations
sd_all <- df %>%
  select(d_edge:d_forest) %>%
  map(sd) %>%
  stack() %>%
  rename(sd = values) %>%
  mutate(ind = paste("s", ind, sep = "."))

# divide scaled coefficients by SD
coef_unst <- coef_st %>%
  inner_join(sd_all, by = c("term" = "ind")) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(. / sd, 7))) %>%
  select(-effect, -(std.error:p.value), -sd)

print(coef_unst)
write.csv(coef_unst, file.path(output_dir, "coef_unst.csv"), row.names = FALSE)


## ---- PART 2: RSF AGREEMENT (Figure 3) --------------------------------------

df.rpts <- read.csv(file.path(annotated_dir, "rpts_1000.csv"), as.is = TRUE)

# Student-t error distribution fit best per DHARMa diagnostics
val.model <- glmmTMB(RSF_2order ~ RSF_3order, data = df.rpts,
                     family = t_family(link = "identity"))
simulateResiduals(val.model, plot = TRUE)
summary(val.model)

# correlation tests (Spearman reported in manuscript)
cor.test(df.rpts$RSF_2order, df.rpts$RSF_3order, method = "pearson")
cor.test(df.rpts$RSF_2order, df.rpts$RSF_3order, method = "spearman")

# Figure 3
agreement_plot <- ggpredict(val.model, "RSF_3order [all]") %>%
  plot() +
  xlab("Third-order relative selection probability") +
  ylab("Second-order relative selection probability") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

pdf(file.path(output_dir, "RSF_agreement_plot.pdf"), width = 5, height = 3.5)
print(agreement_plot)
dev.off()

message("Back-transformation and RSF agreement analysis complete.")

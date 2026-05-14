###
# Bobcat Habitat Selection - Camera trap validation
# Wengert et al. 2026, Ecology and Evolution
#
# Validates RSF-predicted relative selection probability against independent
# camera trap presence/non-detection data from six studies within the study
# extent. Tests whether RSF values are higher at bobcat detection cameras
# than at non-detection cameras.
#
# Two complementary analyses (Table 2, Table S3):
#   1. GLMMs: case ~ RSF + (1|group), binomial, camera study as random effect.
#      Five RSF versions compared: 2nd-order, 3rd-order, product, mean, sum.
#      Models ranked by R2, AIC weight, and performance score (Table 2).
#   2. One-tailed t-tests and Wilcoxon rank sum tests comparing RSF values
#      at presence vs non-detection cameras (Table S3).
#
# Produces Figure 5 (estimated marginal means for 2nd, 3rd, and product models).
#
# Studies with < 5 bobcat presences/non-detections excluded: COOS, LAS, SQ.
#
# Input:  data/camera/cameras_PA_RSF.csv
# Output: output/validation_GLMM_results.doc   (Table 2)
#         output/validation_plot.pdf            (Figure 5)
#
###

library(tidyverse)
library(glmmTMB)
library(MuMIn)
library(performance)
library(DHARMa)
library(cowplot)
library(ggeffects)
library(sjPlot)

camera_dir <- "data/camera"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

df <- read.csv(file.path(camera_dir, "cameras_PA_RSF.csv"), as.is = TRUE)

# exclude studies with < 5 bobcat presences/non-detections
df <- subset(df, !group %in% c("COOS", "LAS", "SQ"))
df$group <- as.factor(df$group)

# standardize RSF values
df$RSF_2order   <- scale(df$RSF_2order)
df$RSF_3order   <- scale(df$RSF_3order)
df$RSF_product  <- scale(df$RSF_product)
df$RSF_addition <- scale(df$RSF_addition)
df$RSF_average  <- scale(df$RSF_average)


## ---- GLMM models (Table 2) --------------------------------------------------

fit_val <- function(rsf_col, df) {
  glmmTMB(as.formula(paste("case ~", rsf_col, "+ (1|group)")),
          family = binomial(), data = df)
}

m.RSF2     <- fit_val("RSF_2order",   df)
m.RSF3     <- fit_val("RSF_3order",   df)
m.RSF_sum  <- fit_val("RSF_addition", df)
m.RSF_mean <- fit_val("RSF_average",  df)
m.RSF_prod <- fit_val("RSF_product",  df)

# model performance ranking (Table 2)
perf <- compare_performance(m.RSF2, m.RSF3, m.RSF_sum, m.RSF_mean, m.RSF_prod,
                            rank = TRUE)
print(perf)

tab_model(m.RSF2, m.RSF3, m.RSF_sum, m.RSF_mean, m.RSF_prod,
          show.intercept = FALSE, digits = 3, transform = NULL,
          file = file.path(output_dir, "validation_GLMM_results.doc"))


## ---- t-tests and Wilcoxon rank sum tests (Table S3) -------------------------

case1 <- subset(df, case == 1)
case0 <- subset(df, case == 0)

run_tests <- function(var, label) {
  cat("\n====", label, "====\n")
  print(t.test(case1[[var]], case0[[var]], alternative = "greater"))
  print(wilcox.test(case1[[var]], case0[[var]], alternative = "greater"))
}

run_tests("RSF_2order",   "Second-order")
run_tests("RSF_3order",   "Third-order")
run_tests("RSF_average",  "Mean")
run_tests("RSF_addition", "Sum")
run_tests("RSF_product",  "Product")


## ---- Figure 5 ---------------------------------------------------------------

make_plot <- function(model, rsf_col, title_str) {
  ggpredict(model, paste0(rsf_col, " [all]")) %>%
    plot() +
    ggtitle(title_str) +
    xlab("Relative selection probability") +
    ylab("Predicted probability of presence") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
}

p.second  <- make_plot(m.RSF2,     "RSF_2order",  "Second-order")
p.third   <- make_plot(m.RSF3,     "RSF_3order",  "Third-order")
p.product <- make_plot(m.RSF_prod, "RSF_product", "Product")

pdf(file.path(output_dir, "validation_plot.pdf"), width = 12, height = 4)
plot_grid(p.second, p.third, p.product, ncol = 3, align = "none")
dev.off()

message("Validation analysis complete.")

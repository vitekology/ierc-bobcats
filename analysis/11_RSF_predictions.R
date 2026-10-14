# RSF predictions at cameras and random points: 2nd order (lp2), 3rd order with predictors
# capped at the 1st-99th percentile of available values (lp3), and scale-integrated RSF
# (int, product of the percentile ranks of the two orders)

library(tidyverse)
library(mgcv)
library(DHARMa)

rd <- function(f, ...) as_tibble(read.csv(f, ...))
wr <- function(d, f) write.csv(d, f, row.names = FALSE)

out_dir <- "output/predictions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261009)

b2 <- with(rd("output/2nd_order/coef_unscaled_2nd_order.csv"), setNames(estimate_unscaled, term))
b3 <- with(rd("output/3rd_order/coef_unscaled_3rd_order.csv"), setNames(estimate_unscaled, term))

hdr    <- names(read.csv("data/annotated/RSF_data_3rd_order.csv.gz", nrows = 1))
avail3 <- read.csv("data/annotated/RSF_data_3rd_order.csv.gz",
                   colClasses = ifelse(hdr %in% c("case", names(b3)), NA, "NULL"))
avail3 <- avail3[avail3$case == 0, names(b3)]
cap_lo <- map_dbl(avail3, \(x) quantile(x, 0.01))
cap_hi <- map_dbl(avail3, \(x) quantile(x, 0.99))
rm(avail3, hdr); invisible(gc())

predict_rsf <- function(d) {
  x2  <- as.matrix(d[paste0("o2_", names(b2))])
  x3  <- as.matrix(d[paste0("o3_", names(b3))])
  colnames(x3) <- names(b3)
  x3c <- x3
  for (v in names(b3)) x3c[, v] <- pmin(pmax(x3[, v], cap_lo[[v]]), cap_hi[[v]])
  d |>
    select(-matches("^o[23]_")) |>
    mutate(lp2 = drop(x2 %*% b2), lp3 = drop(x3c %*% b3), n_capped = rowSums(x3 != x3c)) |>
    bind_cols(as_tibble(x3 != x3c) |> rename_with(\(n) paste0("capped_", n))) |>
    drop_na(lp2, lp3)
}
rnd  <- predict_rsf(rd("data/annotated/predictors_random_points.csv"))
cams <- predict_rsf(rd("data/annotated/predictors_cameras.csv"))

capping <- tibble(term = names(b3), lower = cap_lo, upper = cap_hi,
                  pct_extent_capped = map_dbl(names(b3), \(v) 100 * mean(rnd[[paste0("capped_", v)]])),
                  pct_cameras_capped = map_dbl(names(b3), \(v) 100 * mean(cams[[paste0("capped_", v)]]))) |>
  mutate(across(starts_with("pct"), \(x) round(x, 1)))
wr(capping, file.path(out_dir, "capping_3rd_order.csv"))
message(sprintf("At least one predictor capped: %.1f%% of the extent, %.1f%% of deployments",
                100 * mean(rnd$n_capped > 0), 100 * mean(cams$n_capped > 0)))
rnd  <- select(rnd,  -starts_with("capped_"))
cams <- select(cams, -starts_with("capped_"))

prank <- function(x, ref) ecdf(ref)(x)
integrate <- function(d) mutate(d, p2 = prank(lp2, rnd$lp2), p3 = prank(lp3, rnd$lp3), int = p2 * p3)
rnd  <- integrate(rnd)
cams <- integrate(cams)
wr(rnd,  file.path(out_dir, "predictions_random_points.csv"))
wr(cams, file.path(out_dir, "predictions_cameras.csv"))

## ---- agreement between orders at 1,000 random points --------------------------------------------

agr  <- filter(rnd, pt <= 1000)
m    <- gam(lp2 ~ s(lp3, bs = "tp"), data = agr, method = "REML")
m_lm <- gam(lp2 ~ lp3, data = agr, method = "ML")
m_ml <- gam(lp2 ~ s(lp3, bs = "tp"), data = agr, method = "ML")

pdf(file.path(out_dir, "agreement_gam_check.pdf"), width = 8, height = 8)
gam.check(m)
dev.off()
kc <- k.check(m)

res <- simulateResiduals(m, n = 1000)
pdf(file.path(out_dir, "agreement_dharma.pdf"), width = 9, height = 4.5)
plot(res)
dev.off()
sac <- testSpatialAutocorrelation(res, x = agr$X / 1000, y = agr$Y / 1000, plot = FALSE)

st <- summary(m)
st
sp <- cor.test(agr$lp2, agr$lp3, method = "spearman", exact = FALSE)
agreement <- tibble(n = nrow(agr),
                    edf = st$s.table[1, "edf"], F = st$s.table[1, "F"], p = st$s.table[1, "p-value"],
                    r2_adj = st$r.sq, dev_explained = st$dev.expl,
                    k_index = kc[1, "k-index"], k_index_p = kc[1, "p-value"],
                    dAIC_vs_linear = AIC(m_ml) - AIC(m_lm),
                    spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
                    dharma_p_uniformity = testUniformity(res, plot = FALSE)$p.value,
                    dharma_p_dispersion = testDispersion(res, plot = FALSE)$p.value,
                    dharma_p_outliers = testOutliers(res, plot = FALSE)$p.value,
                    moran_I = unname(sac$statistic[["observed"]]), moran_p = sac$p.value) |>
  mutate(across(where(is.double), \(x) signif(x, 3)))
wr(agreement, file.path(out_dir, "agreement.csv"))
print(t(agreement))

nd   <- tibble(lp3 = seq(min(agr$lp3), max(agr$lp3), length.out = 200))
pr   <- predict(m, newdata = nd, se.fit = TRUE)
line <- mutate(nd, fit = as.numeric(pr$fit), lo = fit - 1.96 * as.numeric(pr$se.fit),
               hi = fit + 1.96 * as.numeric(pr$se.fit))

p <- ggplot() +
  geom_point(data = agr, aes(lp3, lp2), alpha = 0.25, size = 0.8) +
  geom_ribbon(data = line, aes(lp3, ymin = lo, ymax = hi), fill = "grey60", alpha = 0.4) +
  geom_line(data = line, aes(lp3, fit), colour = "black", linewidth = 0.7) +
  labs(x = "Third-order log relative selection", y = "Second-order log relative selection") +
  theme_bw()
ggsave(file.path(out_dir, "agreement_plot.pdf"), p, width = 4.5, height = 4)

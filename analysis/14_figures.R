# Manuscript figures

library(tidyverse)
library(sf)
library(terra)
library(sdmTMB)
source("analysis/settings.R")

GPKG    <- "output/data_locations/bobcat_layers.gpkg"
EXAMPLE <- "KR_M03"
N_SHOW  <- 5000
fig_dir <- "output/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

blue <- "#2a78d6"; orange <- "#eb6834"; ink <- "#1a1a1a"; grey <- "#8a8984"
theme_set(theme_bw(base_size = 10) +
            theme(panel.grid.minor = element_blank(), legend.title = element_blank(),
                  legend.position = "bottom"))

save_fig <- function(p, name, width, height) ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height)
read_map <- function(files) { r <- rast(files); crs(r) <- CRS_RASTER; terra::project(r, "EPSG:3310", method = "near") }

## ---- Study area ---------------------------------------------------------------------

extent <- st_read(EXTENT, quiet = TRUE) |> st_transform(3310)
cams   <- st_read(GPKG, "cameras", quiet = TRUE)
inat   <- st_read(GPKG, "inat_bobcats", quiet = TRUE)
hr     <- st_read(GPKG, "home_ranges_95", quiet = TRUE)

keys <- c("Camera locations", "iNaturalist bobcat records", "GPS-collared bobcats (home ranges)")
p1 <- ggplot() +
  geom_sf(data = extent, fill = "grey95", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = cams, aes(colour = keys[1]), size = 0.3) +
  geom_sf(data = inat, aes(colour = keys[2]), size = 0.7) +
  geom_sf(data = hr,   aes(colour = keys[3]), fill = NA, linewidth = 0.6) +
  scale_colour_manual(values = setNames(c(grey, orange, blue), keys), breaks = keys) +
  guides(colour = guide_legend(ncol = 1))
save_fig(p1, "study_area", 5, 8)

## ---- Home ranges ---------------------------------------------------------------------

hrs <- read_csv("output/homeranges/homerange_summary.csv", show_col_types = FALSE) |>
  mutate(sex = if_else(sex == "F", "Female", "Male"))
p2a <- ggplot(hrs, aes(hr95_km2, fct_reorder(id, hr95_km2), colour = sex, shape = sex)) +
  geom_linerange(aes(xmin = hr95_km2_low, xmax = hr95_km2_high)) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c(Female = orange, Male = blue)) +
  labs(x = expression("95% home range (" * km^2 * ", log scale) with 95% CI"), y = NULL)
save_fig(p2a, "homerange_size", 4.5, 4.5)

ex <- read_csv("data/annotated/used_available_3rd_order.csv.gz", show_col_types = FALSE) |>
  filter(id == EXAMPLE) |>
  mutate(type = if_else(case == 1, "Used (GPS fixes)", "Available"))
ex <- bind_rows(filter(ex, case == 0) |> slice_sample(n = N_SHOW), filter(ex, case == 1))
p2b <- ggplot() +
  geom_sf(data = filter(hr, id == EXAMPLE), fill = NA, colour = ink, linewidth = 0.4) +
  geom_point(data = ex, aes(X, Y, colour = type, size = type)) +
  scale_colour_manual(values = c(Available = "grey75", `Used (GPS fixes)` = blue)) +
  scale_size_manual(values = c(Available = 0.1, `Used (GPS fixes)` = 0.9)) +
  labs(x = NULL, y = NULL, title = EXAMPLE) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
save_fig(p2b, "example_homerange", 4.5, 4.5)

## ---- Habitat selection coefficients -----------------------------------------------------

term_labels <- c(prox_scrub = "Proximity to scrub", prox_edge = "Proximity to forest edge",
                 prox_forest = "Proximity to forest", prox_disturb = "Proximity to disturbance",
                 prox_stream = "Proximity to streams", prox_roads_min = "Proximity to minor roads",
                 prox_roads_maj = "Proximity to major roads", prox_development = "Proximity to development",
                 elevation = "Elevation", slope = "Slope", vrm = "Ruggedness (VRM)", forest_age = "Forest age")
series_levels <- c("2nd order", "3rd order")

coefs <- bind_rows(
  read_csv("output/2nd_order/coefficients_2nd_order.csv", show_col_types = FALSE) |> mutate(series = series_levels[1]),
  read_csv("output/3rd_order/coefficients_3rd_order.csv", show_col_types = FALSE) |> mutate(series = series_levels[2]))
used_terms <- intersect(names(term_labels), coefs$term)
coefs <- coefs |>
  mutate(term   = factor(term_labels[term], levels = rev(unname(term_labels[used_terms]))),
         series = factor(series, levels = series_levels),
         y      = as.numeric(term) + c(0.18, -0.18)[as.integer(series)])

p3 <- ggplot(coefs, aes(estimate, y, colour = series, shape = series)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = 0.3) +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_y_continuous(breaks = seq_along(levels(coefs$term)), labels = levels(coefs$term)) +
  scale_colour_manual(values = setNames(c(blue, ink), series_levels)) +
  scale_shape_manual(values = setNames(c(16, 15), series_levels)) +
  labs(x = "Standardized coefficient (95% CI)", y = NULL) +
  theme(panel.grid.major.y = element_blank())
save_fig(p3, "coefficients", 6, 4.5)

## ---- Prediction maps (percentile of the study extent) ---------------------------------------

pal <- hcl.colors(10, "Blues 3", rev = TRUE)
r <- read_map(file.path("output/maps", paste0(c("second_order", "third_order", "integrated"), "_percentile.tif")))
names(r) <- c("2nd order", "3rd order", "Integrated")
png(file.path(fig_dir, "prediction_maps.png"), width = 9, height = 5.5, units = "in", res = 300)
plot(r, nc = 3, col = pal, breaks = seq(0, 100, 10), axes = FALSE, mar = c(1, 1, 2, 5))
dev.off()

## ---- Predicted probability of detection ----------------------------------------------------
# deployment of median length, largest study, average site

models  <- readRDS("output/validation/models_validation.rds")
val     <- read_csv("output/validation/validation.csv", show_col_types = FALSE)
rnd     <- read_csv("output/predictions/predictions_random_points.csv", show_col_types = FALSE)
surf    <- c(`2nd order` = "lp2", `3rd order` = "lp3", `Integrated RSF` = "int")
cll     <- function(eta) 1 - exp(-exp(eta))

pred_curve <- function(nm) {
  m  <- models[[nm]]
  d  <- m$data
  nd <- tibble(x_raw = seq(quantile(d$x_raw, 0.01), quantile(d$x_raw, 0.99), length.out = 100)) |>
    mutate(x = (x_raw - mean(d$x_raw)) / sd(d$x_raw),
           log_effort = median(d$log_effort),
           study = factor(names(which.max(table(d$study))), levels = levels(d$study)))
  fe  <- tidy(m, effects = "fixed")
  X   <- model.matrix(~ x + log_effort + study, nd)[, fe$term, drop = FALSE]
  V   <- m$sd_report$cov.fixed
  V   <- V[rownames(V) == "b_j", colnames(V) == "b_j"]
  eta <- drop(X %*% fe$estimate) + nd$log_effort
  se  <- sqrt(rowSums((X %*% V) * X))
  tibble(surface = nm, pct = 100 * ecdf(rnd[[surf[[nm]]]])(nd$x_raw),
         p = cll(eta), lo = cll(eta - 1.96 * se), hi = cll(eta + 1.96 * se))
}
curves <- map_dfr(names(surf), pred_curve) |> mutate(surface = factor(surface, levels = names(surf)))

lab <- val |>
  filter(surface %in% names(surf)) |>
  mutate(surface = factor(surface, levels = names(surf)),
         l1 = sprintf("italic(beta) == %.2f*','~~SE == %.2f", slope_per_sd, se),
         l2 = sprintf("italic(z) == %.2f*','~~italic(P) == %.3f", z, p))

p5 <- ggplot(curves, aes(pct, p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
  geom_line(linewidth = 0.6) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = l1), parse = TRUE, hjust = -0.05, vjust = 1.5, size = 3) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = l2), parse = TRUE, hjust = -0.05, vjust = 3.2, size = 3) +
  facet_wrap(~ surface, nrow = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Predicted relative selection (percentile of study extent)", y = "Predicted probability of detection")
save_fig(p5, "validation_detection", 9, 3.2)

## ---- Validation slopes ----------------------------------------------------------------------

p5s <- val |>
  mutate(surface = factor(surface, levels = rev(names(surf)))) |>
  ggplot(aes(slope_per_sd, surface)) +
  geom_vline(xintercept = 0, colour = grey, linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = blue, linewidth = 0.5) +
  geom_point(colour = blue, size = 2) +
  labs(x = "Change in log detection rate per SD of predicted selection (95% CI)", y = NULL) +
  theme(panel.grid.major.y = element_blank())
save_fig(p5s, "validation_slopes", 5.5, 2.2)

## ---- Capped 3rd-order predictors --------------------------------------------------------------

n_capped <- clamp(read_map("output/maps/third_order_n_capped.tif"), upper = 3)
levels(n_capped) <- data.frame(value = 0:3, capped = c("0", "1", "2", "3 or more"))
coltab(n_capped) <- data.frame(value = 0:3, col = c("grey90", "#fdae6b", "#e6550d", "#7f2704"))
png(file.path(fig_dir, "capped_predictors.png"), width = 4.5, height = 6.5, units = "in", res = 300)
plot(n_capped, axes = FALSE, mar = c(1, 1, 2, 7), main = "3rd-order predictors capped")
dev.off()

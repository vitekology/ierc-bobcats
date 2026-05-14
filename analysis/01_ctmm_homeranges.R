###
# Bobcat Habitat Selection - ctmm home range estimation
# Wengert et al. 2026, Ecology and Evolution
#
# Fits continuous-time movement models and computes 95% AKDE home ranges
# for 18 GPS-collared bobcats. Also performs population-level home range
# meta-analysis overall and by sex.
#
# Outlier removal was performed prior to archiving; the cleaned dataset
# provided in data/gps/ reflects post-outlier-removal locations only.
#
# Movement model selected per individual (see Table S1):
#   OU anisotropic:  Hoopa_F02
#   OUF anisotropic: Hoopa_F03, Hoopa_F04, Hoopa_M01, Hoopa_M02, Hoopa_M03,
#                    Hoopa_M08, KR_F02, KR_M01, KR_M02, KR_M05, KR_M06,
#                    Stirling_M01
#   OUF isotropic:   Hoopa_M04, Hoopa_M07, KR_M03, SN_M01
#   Stirling_M01:    OUF anisotropic (fits[[2]]; within 2 AIC of top model)
#
# Input:  data/gps/bobcats_gps_clean.csv
# Output: output/bobcat_AKDE_95.rda  (list of 95% AKDE sf polygons, by individual)
#
# References:
#   Calabrese et al. 2016, Methods Ecol Evol 7:1124-1132
#   Fleming et al. 2015, Ecology 96:1182-1188
#   Noonan et al. 2019, Ecol Monogr 89:e01344
#
###

library(ctmm)
library(sf)

gps_dir    <- "data/gps"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

# load cleaned telemetry data
bobcats.raw <- read.csv(file.path(gps_dir, "bobcats_gps_clean.csv"), as.is = TRUE)
bobcats.gps <- as.telemetry(bobcats.raw)


## ---- helper functions -------------------------------------------------------

# fit ctmm models and compute AKDE for one individual
# model_index: 1 = top AIC model; use 2 for Stirling_M01
fit_akde <- function(tel, model_index = 1) {
  guess   <- ctmm.guess(tel, CTMM = ctmm(error = FALSE), interactive = FALSE)
  fits    <- ctmm.select(tel, guess, verbose = TRUE, cores = 0, trace = 1)
  tm      <- fits[[model_index]]
  akde_ud <- akde(tel, tm, weights = TRUE)
  list(tel = tel, fits = fits, top_model = tm, akde = akde_ud)
}

# extract 95% AKDE polygon as sf
akde_to_sf95 <- function(result, id) {
  poly <- as.sf(result$akde, error = FALSE)
  subset(poly, name == paste(id, "95% est"))
}


## ---- fit models (may take several hours) ------------------------------------

message("Fitting ctmm models -- this may take several hours...")

r.Hoopa_F02    <- fit_akde(bobcats.gps$Hoopa_F02)
r.Hoopa_F03    <- fit_akde(bobcats.gps$Hoopa_F03)
r.Hoopa_F04    <- fit_akde(bobcats.gps$Hoopa_F04)
r.Hoopa_M01    <- fit_akde(bobcats.gps$Hoopa_M01)
r.Hoopa_M02    <- fit_akde(bobcats.gps$Hoopa_M02)
r.Hoopa_M03    <- fit_akde(bobcats.gps$Hoopa_M03)
r.Hoopa_M04    <- fit_akde(bobcats.gps$Hoopa_M04)
r.Hoopa_M07    <- fit_akde(bobcats.gps$Hoopa_M07)
r.Hoopa_M08    <- fit_akde(bobcats.gps$Hoopa_M08)
r.KR_F02       <- fit_akde(bobcats.gps$KR_F02)
r.KR_M01       <- fit_akde(bobcats.gps$KR_M01)
r.KR_M02       <- fit_akde(bobcats.gps$KR_M02)
r.KR_M03       <- fit_akde(bobcats.gps$KR_M03)
r.KR_M04       <- fit_akde(bobcats.gps$KR_M04)
r.KR_M05       <- fit_akde(bobcats.gps$KR_M05)
r.KR_M06       <- fit_akde(bobcats.gps$KR_M06)
r.SN_M01       <- fit_akde(bobcats.gps$SN_M01)
# Stirling_M01: OUF anisotropic within 2 AIC of top -- use fits[[2]]
r.Stirling_M01 <- fit_akde(bobcats.gps$Stirling_M01, model_index = 2)


## ---- population home range meta-analysis ------------------------------------

all.AKDE <- list(
  r.Hoopa_F02$akde,  r.Hoopa_F03$akde,  r.Hoopa_F04$akde,
  r.Hoopa_M01$akde,  r.Hoopa_M02$akde,  r.Hoopa_M03$akde,
  r.Hoopa_M04$akde,  r.Hoopa_M07$akde,  r.Hoopa_M08$akde,
  r.KR_F02$akde,     r.KR_M01$akde,     r.KR_M02$akde,
  r.KR_M03$akde,     r.KR_M04$akde,     r.KR_M05$akde,
  r.KR_M06$akde,     r.SN_M01$akde,     r.Stirling_M01$akde
)

all.AKDE.m <- list(
  r.Hoopa_M01$akde, r.Hoopa_M02$akde, r.Hoopa_M03$akde,
  r.Hoopa_M04$akde, r.Hoopa_M07$akde, r.Hoopa_M08$akde,
  r.KR_M01$akde,   r.KR_M02$akde,   r.KR_M03$akde,
  r.KR_M04$akde,   r.KR_M05$akde,   r.KR_M06$akde,
  r.SN_M01$akde,   r.Stirling_M01$akde
)

all.AKDE.f <- list(
  r.Hoopa_F02$akde, r.Hoopa_F03$akde, r.Hoopa_F04$akde, r.KR_F02$akde
)

bobcats.meta   <- meta(all.AKDE,   sort = TRUE, boot = TRUE)   # all
bobcats.meta.m <- meta(all.AKDE.m, sort = TRUE, boot = FALSE)  # males
bobcats.meta.f <- meta(all.AKDE.f, sort = TRUE, boot = FALSE)  # females

message("Population HR estimates (km2):")
print(bobcats.meta)
print(bobcats.meta.m)
print(bobcats.meta.f)


## ---- export AKDE polygons ---------------------------------------------------

ids <- c("Hoopa_F02","Hoopa_F03","Hoopa_F04","Hoopa_M01","Hoopa_M02","Hoopa_M03",
         "Hoopa_M04","Hoopa_M07","Hoopa_M08","KR_F02","KR_M01","KR_M02","KR_M03",
         "KR_M04","KR_M05","KR_M06","SN_M01","Stirling_M01")

results <- list(r.Hoopa_F02, r.Hoopa_F03, r.Hoopa_F04, r.Hoopa_M01, r.Hoopa_M02,
                r.Hoopa_M03, r.Hoopa_M04, r.Hoopa_M07, r.Hoopa_M08, r.KR_F02,
                r.KR_M01,    r.KR_M02,    r.KR_M03,    r.KR_M04,    r.KR_M05,
                r.KR_M06,    r.SN_M01,    r.Stirling_M01)

akde_95_list <- mapply(akde_to_sf95, results, ids, SIMPLIFY = FALSE)
names(akde_95_list) <- ids
saveRDS(akde_95_list, file.path(output_dir, "bobcat_AKDE_95.rda"))

message("AKDE polygons saved to output/bobcat_AKDE_95.rda")

# Bobcat Habitat Selection in California and Southern Oregon Forests

Code and data repository for:

> Wengert GM, Jirinec V, Higley JM, Gabriel MW, Clifford DL, Thompson CM, Purcell KL, Moriarty KM, Sacks BN, Sweitzer RA. Habitat selection of bobcats in California and southern Oregon forests. *Ecology and Evolution*, 2026.

---

## Overview

This repository contains the R code and data used to estimate bobcat home ranges, build resource selection functions (RSFs) at two spatial scales, and validate those RSFs against independent camera trap data.

Three data sources underpin the analysis:

- **GPS collar data** from 18 bobcats across three study areas in northern California (Hoopa, Kings River/Sierra Nevada, Stirling), 2011–2018
- **Citizen science observations** from iNaturalist/GBIF (*Lynx rufus* and *L. r. fasciatus*, 2010–2020)
- **Camera trap data** from six independent mesocarnivore monitoring studies within the study region (validation only)

**Habitat annotation** of used/available points was performed in ArcGIS Pro 2.9 using raster layers derived from publicly available sources (see Table 1 of the manuscript: NLCD, LandTrendr, NHD, TIGER/Line, NED, GNN). The annotated point datasets — which are the direct inputs to RSF modeling — are provided in `data/annotated/`. The raster layers themselves are not archived here due to file size, but can be reproduced from the sources listed in Table 1 following the procedures described in the manuscript.

---

## Data availability

**GPS collar locations** (cleaned, post-outlier-removal): provided in `data/gps/`. Raw collar files are available upon request from the corresponding author (gwengert@iercecology.org).

**iNaturalist observations** were downloaded from GBIF (accessed 2021): https://www.gbif.org

**Camera trap data** were collected as part of separate fisher and Pacific marten monitoring studies. Contact the corresponding author for access.

---

## Repository structure

```
bobcat-habitat-selection/
├── README.md
├── data/
│   ├── gps/
│   │   └── bobcats_gps_clean.csv         # Cleaned GPS locations (18 bobcats)
│   ├── iNaturalist/
│   │   ├── LR.csv                         # Lynx rufus iNat observations (GBIF)
│   │   ├── LRF.csv                        # Lynx rufus fasciatus iNat observations (GBIF)
│   │   └── study_extent.shp               # Study region boundary
│   ├── annotated/
│   │   ├── bobcats_up_ap_annotated.csv    # 3rd-order used/available + habitat values
│   │   ├── iNat_up_ap_annotated.csv       # 2nd-order used/available + habitat values
│   │   └── rpts_1000.csv                  # Random points with 2nd/3rd-order RSF values
│   └── camera/
│       └── cameras_PA_RSF.csv             # Camera trap presence/non-detection + RSF values
├── analysis/                              # R scripts (numbered in execution order)
└── output/
    └── coef_st.csv                        # Pre-computed profile-likelihood CIs (script 08)
```

---

## Analysis workflow

Scripts are numbered in order of execution. Each reads from `data/` and writes to `data/annotated/` or `output/`.

| Script | Description |
|--------|-------------|
| `01_ctmm_homeranges.R` | Fit ctmm movement models; compute 95% AKDE home ranges; population meta-analysis |
| `02_used_available_3rd_order.R` | Generate used/available points for 3rd-order RSF within AKDE home ranges |
| *(ArcGIS — see manuscript Table 1)* | Extract habitat raster values at used/available points |
| `03_annotated_data_processing_3rd_order.R` | Assign temporally-matched habitat values per individual; format RSF input |
| `04_used_available_2nd_order.R` | Process iNaturalist observations; generate used/available for 2nd-order RSF |
| *(ArcGIS — see manuscript Table 1)* | Extract habitat raster values at iNaturalist used/available points |
| `05_annotated_data_processing_2nd_order.R` | Assign temporally-matched habitat values by year; format RSF input |
| `06_RSF_3rd_order.R` | Fit 3rd-order RSF (GLMM; Muff et al. 2020); scaled and unstandardized versions |
| `07_RSF_2nd_order.R` | Fit 2nd-order RSF (GLM; binomial); scaled and unstandardized versions |
| `08_RSF_backtransform_agreement.R` | Back-transform scaled coefficients; test 2nd/3rd-order RSF correlation (Figure 3) |
| `09_camera_validation.R` | Validate RSFs against camera trap data; GLMMs, t-tests, Wilcoxon tests (Tables 2, S3; Figure 5) |

Scripts 02 and 04 generate the used/available point files that were annotated in ArcGIS. The annotated outputs are already provided in `data/annotated/`, so scripts 02 and 04 do not need to be re-run to reproduce the RSF models and downstream analyses.

---

## R package dependencies

```r
install.packages(c(
  "tidyverse", "lubridate", "sf",          # data handling and spatial
  "ctmm",                                   # home range estimation
  "glmmTMB", "MuMIn", "broom.mixed",       # RSF modeling
  "performance", "DHARMa", "psych",         # model diagnostics
  "sjPlot", "ggeffects", "cowplot",         # visualization
  "tictoc"                                  # timing
))
```

Note: `ggeffects::ggpredict()` is used for marginal predictions (scripts 08 and 09). The `effects` package is not required.

**Note on deprecated packages:** The original analysis used `rgdal` and `rgeos`, which were retired from CRAN in October 2023. All spatial operations in these scripts use `sf` equivalents.

---

## Computational notes

- `ctmm.select()` in script 01 is computationally intensive. Fitting all 18 individuals sequentially may take several hours. Population-level home range results from `meta()` are printed directly to console.
- Profile-likelihood CIs in script 08 take approximately 24 hours to compute. Pre-computed results are provided as `output/coef_st.csv`; see the script header for instructions to recompute.

---

## Citation

Please cite both the manuscript and this repository:

> Wengert GM, Jirinec V, et al. (2026). Code and data for: Habitat selection of bobcats in California and southern Oregon forests. Zenodo. https://doi.org/[DOI pending]

---

## Contact

Greta M. Wengert — gwengert@iercecology.org  
Integral Ecology Research Center, Blue Lake, CA

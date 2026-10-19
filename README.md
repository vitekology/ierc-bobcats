# Multiscale Habitat Selection of Bobcats in Forested Regions of California and Southern Oregon

Code and data for:

> Wengert GM, Jirinec V, Higley JM, Gabriel MW, Clifford DL, Thompson CM, Purcell KL, Moriarty KM, Sacks BN, Sweitzer RA. Multiscale habitat selection of bobcats in forested regions of California and southern Oregon. *Ecology and Evolution*, 2026.

---

## Overview

The R scripts in `analysis/` estimate bobcat home ranges, fit resource selection functions (RSFs) at two scales, combine them into a scale-integrated RSF, map all three across the study extent, and compare them with independent camera-trap detections.

- **Second order:** research-grade iNaturalist bobcat records against a target-group background (1-km cells with at least one record of a bobcat or of 17 other medium- to large-bodied mammals). Weighted logistic regression.
- **Third order:** GPS locations of 18 collared bobcats against random points within each animal's 95% AKDE home range. Weighted logistic mixed model with random slopes per animal (Muff et al. 2020).
- **Scale-integrated RSF:** product of the percentile ranks of the second- and third-order predictions.
- **Camera comparison:** detection during each camera deployment modeled as a function of predicted relative selection (binomial GLMM, complementary log-log link, log effort as an offset plus a log-effort covariate so that its exponent is estimated, camera study as a fixed effect, site as a random intercept; sdmTMB).

The study extent (`data/iNaturalist/study_extent_v3.shp`) consists of level III ecoregions in northern California and southern Oregon.

---

## Data

| Data | Location | Notes |
|---|---|---|
| GPS locations | `data/gps/bobcats_gps_clean.csv` | 18 bobcats (Hoopa, Sierra/Kings River, Stirling), 2011–2018; 11,943 locations after outlier removal. Raw collar files are available from the corresponding author. |
| iNaturalist records | `data/iNaturalist/` | Downloaded from the iNaturalist API on 28 Sep 2026 (raw pulls `*_raw_20260928.csv`) and filtered by script 01 (`inat_bobcat_filtered.csv`, `inat_targetgroup_filtered.csv`). Research grade, 2010–2020. Each record keeps its observer's licence (`license` column). |
| Camera deployments | `data/camera/` | 2,695 deployments at 1,492 sites from seven studies, 2007–2020 (`camera_deployments.csv`; summary in `camera_studies.csv`). Site names are anonymised. The raw camera files belong to the contributing studies and are not included; contact the corresponding author. |
| Annotated points | `data/annotated/` | Used and available points with habitat values (the direct inputs to the RSFs), habitat values at camera deployments and at 10,000 random points, and the available points of the example home range (`available_KR_M03.gpkg`). |
| Habitat rasters | not included | Derived from NLCD, LandTrendr, NHD, TIGER/Line, NED and GNN (manuscript Table 1). Not archived because of their size; scripts 02, 06–08 and the maps in 13 need them. |

---

## Repository structure

```
ierc-bobcats/
├── README.md
├── ierc-bobcats.Rproj
├── analysis/
│   ├── 100_master_run.R              # runs the scripts in order
│   ├── settings.R                    # shared settings and helper functions
│   ├── 01_download_inat.R … 14_figures.R
│   └── paths_local.R                 # local data paths (not tracked; see below)
└── data/
    ├── annotated/
    │   ├── used_available_2nd_order.csv
    │   ├── used_available_3rd_order.csv.gz
    │   ├── RSF_data_2nd_order.csv
    │   ├── RSF_data_3rd_order.csv.gz
    │   ├── predictors_cameras.csv
    │   ├── predictors_random_points.csv
    │   └── available_KR_M03.gpkg
    ├── camera/
    │   ├── camera_deployments.csv
    │   └── camera_studies.csv
    ├── gps/
    │   └── bobcats_gps_clean.csv
    └── iNaturalist/
        ├── inat_bobcat_raw_20260928.csv
        ├── inat_mammal_raw_20260928.csv
        ├── inat_bobcat_filtered.csv
        ├── inat_targetgroup_filtered.csv
        └── study_extent_v3.*         # study extent (shapefile)
```

All outputs are written to `output/`, which is not tracked.

---

## Analysis workflow

Run the scripts from the repository root, in order, with `analysis/100_master_run.R` (set `RUN` to the scripts to run). Each script runs in a fresh R session.

| Script | Does | Main outputs | Needs data outside the repository |
|---|---|---|---|
| `01_download_inat.R` | Download (`DOWNLOAD <- TRUE`) or re-filter the saved iNaturalist pulls (default) | `data/iNaturalist/*_filtered.csv` | no |
| `02_export_landcover.R` | Export the habitat rasters from the geodatabases to GeoTIFF | rasters in `RASTER_ROOT` | yes (geodatabases) |
| `03_ctmm_homeranges.R` | Movement models and 95% AKDE home ranges (ctmm); population mean | `output/homeranges/` | no |
| `04_used_available_2nd_order.R` | Target-group cells and 2nd-order available points | `data/annotated/used_available_2nd_order.csv`, `output/2nd_order/targetgroup_cells.gpkg` | no |
| `05_used_available_3rd_order.R` | 3rd-order available points within the home ranges | `data/annotated/used_available_3rd_order.csv.gz` | no |
| `06_extract_landcover_2nd_order.R` | Habitat values at 2nd-order points, matched to record year | `data/annotated/RSF_data_2nd_order.csv` | yes (rasters) |
| `07_extract_landcover_3rd_order.R` | Habitat values at 3rd-order points, matched to each animal's tracking period | `data/annotated/RSF_data_3rd_order.csv.gz` | yes (rasters) |
| `08_extract_landcover_predictions.R` | Habitat values at camera deployments and at 10,000 random points | `data/annotated/predictors_*.csv` | yes (rasters) |
| `09_RSF_2nd_order.R` | Collinearity screening and 2nd-order RSF | `output/2nd_order/` | no |
| `10_RSF_3rd_order.R` | Collinearity screening and 3rd-order RSF | `output/3rd_order/` | no |
| `11_RSF_predictions.R` | Predictions at cameras and random points; capping; integrated RSF; agreement GAM | `output/predictions/` | no |
| `12_camera_validation.R` | Choice of model structure, camera comparison, residual checks, camera coverage | `output/validation/` | no |
| `13_export_gis_layers.R` | GeoPackage of data locations; prediction maps | `output/data_locations/`, `output/maps/` | yes for the maps (rasters) |
| `14_figures.R` | Figures | `output/figures/` | needs the maps from 13 |

**Reproducing the results from the archived data.** The habitat values are already in `data/annotated/`, so the models, predictions and camera comparison can be rerun without the rasters:

```r
RUN <- c(1, 3:5, 9:12)   # in analysis/100_master_run.R
```

Scripts 04 and 05 use fixed random seeds and regenerate the archived used/available files.

**Local paths.** Scripts 02, 06–08 and 13 read habitat data kept outside the repository. Either place the data in the default folders under `external/` (set in `analysis/settings.R`) or create `analysis/paths_local.R` (not tracked) with:

```r
GDB_DIR     <- "path/to/landcover/geodatabases"
RASTER_ROOT <- "path/to/exported/rasters"
```

---

## Software

R 4.6.1 with tidyverse, sf (1.1-3), terra (1.9-50), ctmm (1.3.0), future and future.apply, glmmTMB (1.1.15.2), broom.mixed, mgcv (1.9-4), sdmTMB (1.1.0), DHARMa (0.5.0), httr2 (script 01) and callr (master script).

```r
install.packages(c("tidyverse", "sf", "terra", "ctmm", "future", "future.apply", "glmmTMB",
                   "broom.mixed", "mgcv", "sdmTMB", "DHARMa", "httr2", "callr"))
```

---

## Computational notes

- Script 03 fits the movement models in parallel (`N_WORKERS` in the script) and can take several hours.
- Script 10 fits a mixed model to about 1.2 million points and needs several GB of memory.
- Script 12 compares eight model structures, including spatial random fields, and takes a while.

---

## Contact

Greta M. Wengert, gwengert@iercecology.org
Integral Ecology Research Center, Blue Lake, CA

# Runs the analysis scripts in order, each in a fresh R session (from the repository root).
# Scripts 00 and 02 need data kept outside the repository (paths in analysis/settings.R).

RUN <- 0:14

scripts <- c(
  "00_camera_compilation.R",
  "01_download_inat.R",
  "02_export_landcover.R",
  "03_ctmm_homeranges.R",
  "04_used_available_2nd_order.R",
  "05_used_available_3rd_order.R",
  "06_extract_landcover_2nd_order.R",
  "07_extract_landcover_3rd_order.R",
  "08_extract_landcover_predictions.R",
  "09_RSF_2nd_order.R",
  "10_RSF_3rd_order.R",
  "11_RSF_predictions.R",
  "12_camera_validation.R",
  "13_export_gis_layers.R",
  "14_figures.R"
)

stopifnot(dir.exists("analysis"), dir.exists("data"))
scripts <- scripts[as.integer(substr(scripts, 1, 2)) %in% RUN]

timing <- data.frame(script = scripts, minutes = NA_real_)
for (i in seq_along(scripts)) {
  message("\n==== ", scripts[i])
  t0 <- Sys.time()
  callr::rscript(file.path("analysis", scripts[i]), show = TRUE, stderr = "2>&1")
  timing$minutes[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
}
print(timing)

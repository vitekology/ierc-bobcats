# Home ranges of the collared bobcats (ctmm, 95% AKDE) and population mean area

library(tidyverse)
library(ctmm)
library(sf)
library(future)
library(future.apply)

set.seed(2026)

out_dir    <- "output/homeranges"
N_WORKERS  <- 12
MODEL_RANK <- c(Stirling_M01 = 2)
dir.create(file.path(out_dir, "fits"), showWarnings = FALSE, recursive = TRUE)

tel_list <- as.telemetry(read.csv("data/gps/bobcats_gps_clean.csv", as.is = TRUE))
ids      <- names(tel_list)
n_fixes  <- sapply(tel_list, nrow)

jobs <- lapply(ids[order(n_fixes, decreasing = TRUE)], function(i) list(
  id = i, tel = tel_list[[i]],
  rank = if (i %in% names(MODEL_RANK)) MODEL_RANK[[i]] else 1,
  file = file.path(out_dir, "fits", paste0(i, ".rds"))))

plan(multisession, workers = N_WORKERS)
status <- future_lapply(jobs, function(job) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  tryCatch({
    guess <- ctmm.guess(job$tel, CTMM = ctmm(error = FALSE), interactive = FALSE)
    fits  <- ctmm.select(job$tel, guess, verbose = TRUE, cores = 1)
    ud    <- akde(job$tel, fits[[job$rank]], weights = TRUE)
    saveRDS(list(id = job$id, tel = job$tel, fits = fits, akde = ud), job$file)
    list(id = job$id, model = names(fits)[job$rank], error = NA_character_)
  }, error = function(e) list(id = job$id, model = NA_character_, error = conditionMessage(e)))
}, future.seed = TRUE, future.chunk.size = 1) |> bind_rows()
plan(sequential)
if (any(!is.na(status$error))) { print(status); stop("Some animals failed.") }

uds <- setNames(lapply(ids, \(i) readRDS(file.path(out_dir, "fits", paste0(i, ".rds")))$akde), ids)

hr_summary <- map_dfr(ids, \(i) {
  ci <- summary(uds[[i]], units = FALSE)$CI
  tibble(id = i, sex = if_else(grepl("_F", i), "F", "M"), n_fixes = n_fixes[[i]],
         model = status$model[status$id == i],
         hr95_km2_low = ci[1, "low"] / 1e6, hr95_km2 = ci[1, "est"] / 1e6, hr95_km2_high = ci[1, "high"] / 1e6)
}) |> mutate(across(starts_with("hr95"), \(x) round(x, 1)))
write_csv(hr_summary, file.path(out_dir, "homerange_summary.csv"))
print(hr_summary, n = Inf)

is_f <- grepl("_F", ids)
meta_all <- meta(uds,        sort = TRUE, boot = TRUE)
meta_m   <- meta(uds[!is_f], sort = TRUE, boot = FALSE)
meta_f   <- meta(uds[is_f],  sort = TRUE, boot = FALSE)
capture.output(cat("== All animals ==\n"), print(meta_all),
               cat("\n== Males ==\n"),     print(meta_m),
               cat("\n== Females ==\n"),   print(meta_f),
               file = file.path(out_dir, "meta_analysis.txt"))

akde_95 <- setNames(lapply(ids, \(i) {
  poly <- as.sf(uds[[i]], error = FALSE)
  poly[poly$name == paste(i, "95% est"), ]
}), ids)
saveRDS(akde_95, file.path(out_dir, "akde_95.rds"))

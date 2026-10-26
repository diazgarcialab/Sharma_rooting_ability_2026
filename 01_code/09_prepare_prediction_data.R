# ==============================================================================
# 09_prepare_prediction_data.R -- assemble the inputs for genomic and phenomic
# prediction: harmonised IDs, one spectrum per accession, and the marker sets.
#
#   Rscript 01_code/09_prepare_prediction_data.R
#
# Inputs (00_raw_data/)
#   spectra_interpolated.csv    1408 scans x 721 wavelengths (X380..X1100) + Tags
#   spectra_dictionary.csv      Tags (sample_id) -> genotype
#   genotype_matrix_{5,10,20,50}k.raw   PLINK --recode A, accession names as IID
#
# Outputs (02_clean_data/prediction/)
#   spectra_by_accession.csv    one averaged spectrum per accession
#   spectra_scan_qc.csv         per-scan QC, so a dropped scan is traceable
#   prediction_panel.csv        which accession has phenotypes / markers / spectra
#   id_harmonisation.csv        every ID change made, with the rule that made it
#
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/09_prepare_prediction_data.R"
  }))), "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(tibble)
})

ACCEPT_UT12_TYPO <- FALSE

DIR_PRED <- file.path(DIR_CLEAN, "prediction")
dir.create(DIR_PRED, showWarnings = FALSE, recursive = TRUE)

# Underscores to spaces, collapse runs of whitespace, case-fold. Nothing fuzzier.
norm_id <- function(x) {
  x <- str_replace_all(as.character(x), "_", " ")
  tolower(str_squish(x))
}

# ---- 1. the phenotype side ---------------------------------------------------
blups <- read_csv(file.path(DIR_CLEAN, "blups_by_trait.csv"), show_col_types = FALSE)
pheno <- read_csv(file.path(DIR_CLEAN, "pheno_long.csv"), show_col_types = FALSE)
panel <- tibble(genotype = blups$genotype, key = norm_id(blups$genotype))
stopifnot(!any(duplicated(panel$key)))
message("phenotyped accessions: ", nrow(panel))

id_log <- tibble(source = character(), original = character(),
                 resolved = character(), rule = character())
# Only STRUCTURAL changes are logged. Case-folding touches nearly every ID and
# would bury the handful of edits that actually matter under a thousand rows.
log_ids <- function(src, orig, res) {
  structural <- str_squish(str_replace_all(orig, "_", " ")) != str_squish(orig)
  if (!any(structural)) return(invisible(NULL))
  id_log <<- bind_rows(id_log, tibble(
    source = src, original = orig[structural], resolved = res[structural],
    rule = "underscore -> space"))
}

# ---- 2. spectra --------------------------------------------------------------
spec_raw <- read_csv(file.path(DIR_RAW, "spectra_interpolated.csv"),
                     show_col_types = FALSE, progress = FALSE)
dict <- read_csv(file.path(DIR_RAW, "spectra_dictionary.csv"), show_col_types = FALSE)
wl_cols <- grep("^X[0-9]+$", names(spec_raw), value = TRUE)
WL <- as.numeric(sub("^X", "", wl_cols))
message("spectral scans: ", nrow(spec_raw), " over ", length(wl_cols),
        " wavelengths (", min(WL), "-", max(WL), " nm)")

# Two scans are labelled "Example ..." -- instrument test shots, not plants.
n_example <- sum(!spec_raw$Tags %in% as.character(dict$sample_id))
spec_raw <- spec_raw |> filter(Tags %in% as.character(dict$sample_id))
message("dropped ", n_example, " scans whose Tag is not a sample_id (instrument tests)")

dict <- dict |> mutate(key = norm_id(genotype))
log_ids("spectra_dictionary", dict$genotype, dict$key)

unmatched <- dict |> filter(!key %in% panel$key)
if (nrow(unmatched)) {
  message("\nspectral samples with no phenotyped accession (", nrow(unmatched), "):")
  print(as.data.frame(unmatched |> select(sample_id, genotype)), row.names = FALSE)
  if (ACCEPT_UT12_TYPO && any(unmatched$genotype == "UT12-12-084")) {
    dict <- dict |> mutate(key = if_else(genotype == "UT12-12-084", "ut12-084", key))
    id_log <- bind_rows(id_log, tibble(source = "spectra_dictionary",
      original = "UT12-12-084", resolved = "UT12-084",
      rule = "MANUAL: ACCEPT_UT12_TYPO = TRUE"))
    message("  -> UT12-12-084 mapped to UT12-084 (ACCEPT_UT12_TYPO is TRUE)")
  } else {
    message("  -> dropped. See the note at the top of this script.")
  }
}

# Per-scan QC before averaging: a scan whose reflectance is flat, negative or
# wildly off its siblings is an instrument artefact, and averaging it in would
# quietly corrupt that accession's spectrum.
M <- as.matrix(spec_raw[, wl_cols])
scan_qc <- tibble(
  row      = seq_len(nrow(M)),
  Tags     = spec_raw$Tags,
  min_refl = apply(M, 1, min),
  sd_refl  = apply(M, 1, sd),
  n_na     = rowSums(is.na(M)))
scan_qc <- scan_qc |>
  group_by(Tags) |>
  mutate(cor_to_median = {
    med <- apply(M[row, , drop = FALSE], 2, median)
    vapply(row, \(i) suppressWarnings(cor(M[i, ], med)), numeric(1))
  }) |>
  ungroup() |>
  mutate(flag_na      = n_na > 0,
         flag_flat    = sd_refl < 1e-8,
         flag_outlier = !is.na(cor_to_median) & cor_to_median < 0.95,
         keep         = !(flag_na | flag_flat | flag_outlier))
write_csv(scan_qc, file.path(DIR_PRED, "spectra_scan_qc.csv"))
message("\nscan QC: ", sum(!scan_qc$keep), " of ", nrow(scan_qc), " scans flagged",
        "  (NA ", sum(scan_qc$flag_na), ", flat ", sum(scan_qc$flag_flat),
        ", low correlation to the sample median ", sum(scan_qc$flag_outlier), ")")

spec_ok <- spec_raw[scan_qc$keep, ]
# sample_id is numeric in the dictionary and character in the scan table, so it
# is cast once here rather than relying on the join to coerce it.
spec_by_acc <- spec_ok |>
  mutate(Tags = as.character(Tags)) |>
  left_join(dict |> mutate(sample_id = as.character(sample_id)) |>
              select(Tags = sample_id, key), by = "Tags") |>
  filter(key %in% panel$key) |>
  group_by(key) |>
  summarise(n_scans = n(), across(all_of(wl_cols), mean), .groups = "drop")

message("accessions with a spectrum: ", nrow(spec_by_acc),
        "  (median ", median(spec_by_acc$n_scans), " scans each)")
write_csv(spec_by_acc |> left_join(panel, by = "key") |>
            select(genotype, key, n_scans, all_of(wl_cols)),
          file.path(DIR_PRED, "spectra_by_accession.csv"))

# ---- 3. marker matrices ------------------------------------------------------
DENSITIES <- c("5k", "10k", "20k", "50k")
geno_keys <- list()
for (d in DENSITIES) {
  f <- file.path(DIR_RAW, paste0("genotype_matrix_", d, ".raw"))
  if (!file.exists(f)) { message("missing ", basename(f), " -- skipped"); next }
  hdr <- scan(f, what = "", nlines = 1, quiet = TRUE)
  ids <- read_table(f, col_types = cols_only(IID = col_character()), progress = FALSE)$IID
  k <- norm_id(ids)
  log_ids(paste0("genotype_matrix_", d), ids, k)
  geno_keys[[d]] <- k
  message(sprintf("%-4s markers: %d samples x %d SNPs, %d matched to a phenotyped accession",
                  d, length(k), length(hdr) - 6, sum(k %in% panel$key)))
}
stopifnot(length(geno_keys) > 0)
if (length(geno_keys) > 1) {
  same <- vapply(geno_keys[-1], \(x) identical(x, geno_keys[[1]]), logical(1))
  if (!all(same))
    warning("the marker files do not all hold the same samples in the same order")
  else message("all marker files hold the same ", length(geno_keys[[1]]),
               " samples in the same order")
}
gkey <- geno_keys[[1]]

# ---- 4. who has what ---------------------------------------------------------
pheno_years <- pheno |>
  mutate(key = norm_id(genotype)) |>
  distinct(key, year) |>
  mutate(has = TRUE) |>
  pivot_wider(names_from = year, values_from = has, names_prefix = "pheno_",
              values_fill = FALSE)

panel_full <- panel |>
  left_join(pheno_years, by = "key") |>
  mutate(has_markers = key %in% gkey,
         has_spectra = key %in% spec_by_acc$key) |>
  left_join(blups |> transmute(key = norm_id(genotype), species), by = "key")
write_csv(panel_full, file.path(DIR_PRED, "prediction_panel.csv"))
write_csv(id_log, file.path(DIR_PRED, "id_harmonisation.csv"))

n_both <- sum(panel_full$has_markers & panel_full$has_spectra)
cat("\n"); cat(strrep("=", 72), "\n")
cat("PREDICTION PANEL\n"); cat(strrep("=", 72), "\n")
cat(sprintf("  phenotyped accessions            : %d\n", nrow(panel_full)))
cat(sprintf("  with markers                     : %d\n", sum(panel_full$has_markers)))
cat(sprintf("  with spectra                     : %d\n", sum(panel_full$has_spectra)))
cat(sprintf("  with BOTH markers and spectra    : %d\n", n_both))
cat(sprintf("  species represented              : %d\n",
            n_distinct(panel_full$species[panel_full$has_markers])))
cat(sprintf("  IDs needing underscore repair    : %d across %d files (see id_harmonisation.csv)\n",
            n_distinct(id_log$original), n_distinct(id_log$source)))
cat("\n  by year (accessions with markers):\n")
for (y in c("2024", "2025", "2026")) {
  col <- paste0("pheno_", y)
  if (col %in% names(panel_full))
    cat(sprintf("    %s : %d phenotyped, %d also genotyped, %d also with spectra\n", y,
                sum(panel_full[[col]], na.rm = TRUE),
                sum(panel_full[[col]] & panel_full$has_markers, na.rm = TRUE),
                sum(panel_full[[col]] & panel_full$has_markers &
                    panel_full$has_spectra, na.rm = TRUE)))
}
message("\nwritten to ", DIR_PRED)

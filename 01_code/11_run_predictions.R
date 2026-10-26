# ==============================================================================
# 11_run_predictions.R -- genomic, phenomic and taxonomy-based prediction
#
#   Rscript 01_code/11_run_predictions.R
#   NREP=5 Rscript 01_code/11_run_predictions.R      # quicker, noisier
#
# Requires 09_prepare_prediction_data.R to have run. 10_prediction_engine.R holds
# the model code and is sourced; 12_verify_prediction.R checks that engine and
# should pass before any of these numbers are believed.
#
# Outputs (02_clean_data/prediction/)
#   prediction_accuracy.csv     one row per scenario x trait x model x CV design
#   prediction_predictions.csv  the held-out predictions themselves
#   spectral_preprocessing.csv  raw vs SNV vs derivatives
#   learning_curve.csv          accuracy against training-set size
#   kernels_pca.csv             PCA scores for the figure
#
# ==============================================================================

CODE_DIR <- dirname(normalizePath(local({
  for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) f[1] else "01_code/11_run_predictions.R"
})))
source(file.path(CODE_DIR, "00_setup.R"))
source(file.path(CODE_DIR, "10_prediction_engine.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(tibble)
})

DIR_PRED <- file.path(DIR_CLEAN, "prediction")
if (!dir.exists(DIR_PRED))
  stop("run 09_prepare_prediction_data.R first")

NREP      <- as.integer(Sys.getenv("NREP", unset = "20"))
NFOLD     <- 5L
SEED      <- 20260819
DENSITIES <- c("5k", "10k", "20k", "50k")
MODEL_DENSITY <- "10k"        # density used for the model comparison, as in the draft
MIN_SPECIES_N <- 5L           # a species smaller than this cannot be held out alone
# "random" by default; CV_DESIGNS="random,loso" re-enables leave-one-species-out.
CV_DESIGNS <- trimws(strsplit(Sys.getenv("CV_DESIGNS", unset = "random"), ",")[[1]])
stopifnot(all(CV_DESIGNS %in% c("random", "loso")))

TRAITS <- c("Callus-Stage", "Post-transplant", "Root weight",
            "Cutting weight", "Transplant response")

set.seed(SEED)

# ---- 1. load -----------------------------------------------------------------
panel  <- read_csv(file.path(DIR_PRED, "prediction_panel.csv"), show_col_types = FALSE)
pheno  <- read_csv(file.path(DIR_CLEAN, "pheno_long.csv"), show_col_types = FALSE) |>
  mutate(key = tolower(str_squish(str_replace_all(genotype, "_", " "))))
spec   <- read_csv(file.path(DIR_PRED, "spectra_by_accession.csv"), show_col_types = FALSE)
wl_cols <- grep("^X[0-9]+$", names(spec), value = TRUE)

read_geno <- function(d) {
  f <- file.path(DIR_RAW, paste0("genotype_matrix_", d, ".raw"))
  x <- read_table(f, progress = FALSE, col_types = cols(
    FID = col_character(), IID = col_character(), PAT = col_character(),
    MAT = col_character(), SEX = col_character(), PHENOTYPE = col_character(),
    .default = col_double()))
  keys <- tolower(str_squish(str_replace_all(x$IID, "_", " ")))
  Z <- as.matrix(x[, setdiff(names(x), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))])
  rownames(Z) <- keys
  Z
}

message("building kernels ...")
Zlist <- lapply(DENSITIES, read_geno); names(Zlist) <- DENSITIES
Glist <- lapply(Zlist, kernel_G)
gkeys <- rownames(Zlist[[1]])

Xspec <- as.matrix(spec[, wl_cols]); rownames(Xspec) <- spec$key
Sfull <- kernel_S(Xspec)

sp_of <- setNames(panel$species, panel$key)
Tfull <- kernel_T(sp_of[gkeys]); dimnames(Tfull) <- list(gkeys, gkeys)
for (d in DENSITIES) dimnames(Glist[[d]]) <- list(gkeys, gkeys)
message("  G kernels: ", paste(DENSITIES, collapse = ", "),
        " on ", length(gkeys), " accessions")
message("  S kernel : ", nrow(Sfull), " accessions, ", length(wl_cols), " wavelengths")

# ---- 2. phenotype vectors per scenario ---------------------------------------
# Within-year scenarios use that year's observed value; the across-year scenario
# trains on 2024+2025 and validates on accessions seen ONLY in 2026, so the
# validation genotypes are genuinely new.
year_value <- function(trait, yr) {
  pheno |> filter(trait == !!trait, year == !!yr) |>
    group_by(key) |> summarise(v = mean(value, na.rm = TRUE), .groups = "drop")
}

# The across-year design trains on one year and validates on accessions that
# appear ONLY in a later year, so the validation genotypes are genuinely new.
#
# Training on 2024 AND 2025 together does not work: every accession phenotyped in
# 2026 had already been phenotyped in 2024 or 2025, so pooling both training
# years leaves nothing unseen to validate on. Training on 2024 alone leaves
# 95-98 genotyped accessions unique to each later year, which is the design the
# draft used and the one kept here.
scenarios <- list(
  within_2024 = list(type = "within", year = "2024"),
  within_2025 = list(type = "within", year = "2025"),
  within_2026 = list(type = "within", year = "2026"),
  across_24_25 = list(type = "across", train = "2024", test = "2025"),
  across_24_26 = list(type = "across", train = "2024", test = "2026"))

# ---- 3. CV designs -----------------------------------------------------------
folds_random <- function(keys, nfold, rep_seed) {
  set.seed(rep_seed)
  split(seq_along(keys), sample(rep(seq_len(nfold), length.out = length(keys))))
}
folds_species <- function(keys) {
  sp <- sp_of[keys]
  tab <- table(sp)
  big <- names(tab)[tab >= MIN_SPECIES_N]
  out <- lapply(big, \(s) which(sp == s))
  names(out) <- big
  out
}

# ---- 4. the runner -----------------------------------------------------------
run_cv <- function(y, keys, Ks, design, label) {
  # y indexed like keys; Ks a named list of kernels already subset to keys
  folds <- if (design == "random") NULL else folds_species(keys)
  reps  <- if (design == "random") seq_len(NREP) else 1L
  out <- list()
  for (r in reps) {
    fl <- if (design == "random") folds_random(keys, NFOLD, SEED + r * 977) else folds
    pred <- rep(NA_real_, length(y))
    for (f in seq_along(fl)) {
      te <- fl[[f]]; tr <- setdiff(seq_along(y), te)
      tr <- tr[is.finite(y[tr])]
      if (length(tr) < 20 || !length(te)) next
      pred[te] <- if (length(Ks) == 1L) gblup_1k(y, Ks[[1]], tr, te)
                  else                  gblup_mk(y, Ks, tr, te)
    }
    m <- pred_metrics(y, pred)
    # Within-fold correlation, averaged over folds and weighted by fold size.
    # This is the headline for LOSO; see the note at the top of the file.
    wf <- vapply(fl, function(ix) {
      o <- y[ix]; p <- pred[ix]
      keep <- is.finite(o) & is.finite(p)
      if (sum(keep) < 5 || sd(p[keep]) < 1e-10 || sd(o[keep]) < 1e-10) return(c(NA, 0))
      c(cor(o[keep], p[keep]), sum(keep))
    }, numeric(2))
    rw <- wf[1, ]; nw <- wf[2, ]
    ok_f <- is.finite(rw) & nw > 0
    r_within <- if (any(ok_f)) sum(rw[ok_f] * nw[ok_f]) / sum(nw[ok_f]) else NA_real_
    out[[length(out) + 1]] <- tibble(
      rep = r, r_pooled = m[["r"]], r_within = r_within,
      n_folds_scored = sum(ok_f), n_folds = length(fl),
      slope = m[["slope"]], rmse = m[["rmse"]], n = m[["n"]],
      fold_label = if (design == "random") NA_character_
                   else paste(names(fl), collapse = ";"))
  }
  bind_rows(out) |> mutate(label = label, design = design, .before = 1)
}

acc <- list(); pred_store <- list()
t0 <- Sys.time()
job <- 0

# ---- 4a. marker density, genomic model only ----------------------------------
message("\n[1/4] marker density x scenario x trait")
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  for (tr in TRAITS) {
    if (sc$type == "within") {
      v <- year_value(tr, sc$year)
      if (!nrow(v)) next
      keys <- intersect(v$key, gkeys)
      if (length(keys) < 40) next
      y <- setNames(v$v, v$key)[keys]
      te_mask <- NULL
    } else {
      vtr <- bind_rows(lapply(sc$train, \(yy) year_value(tr, yy))) |>
        group_by(key) |> summarise(v = mean(v), .groups = "drop")
      vte <- year_value(tr, sc$test)
      if (!nrow(vtr) || !nrow(vte)) next
      new_keys <- setdiff(vte$key, vtr$key)          # never seen in training years
      keys <- intersect(c(vtr$key, new_keys), gkeys)
      if (length(keys) < 40 || length(intersect(new_keys, gkeys)) < 10) next
      y <- c(setNames(vtr$v, vtr$key), setNames(vte$v, vte$key))[keys]
      te_mask <- keys %in% new_keys
    }
    for (d in DENSITIES) {
      K <- Glist[[d]][keys, keys, drop = FALSE]
      if (sc$type == "within") {
        for (des in CV_DESIGNS) {
          a <- run_cv(y, keys, list(K), des, paste(sc_name, tr, "G", d, sep = "|"))
          if (nrow(a)) acc[[length(acc) + 1]] <- a |>
            mutate(scenario = sc_name, trait = tr, model = "G", density = d)
        }
      } else {
        te <- which(te_mask); trn <- which(!te_mask)
        trn <- trn[is.finite(y[trn])]
        p <- gblup_1k(y, K, trn, te)
        m <- pred_metrics(y[te], p)
        acc[[length(acc) + 1]] <- tibble(
          design = "holdout", rep = 1L, r_pooled = m[["r"]], r_within = NA_real_,
          n_folds_scored = 1L, n_folds = 1L,
          slope = m[["slope"]], rmse = m[["rmse"]], n = m[["n"]],
          fold_label = NA_character_,
          label = paste(sc_name, tr, "G", d, sep = "|"),
          scenario = sc_name, trait = tr, model = "G", density = d)
        if (d == MODEL_DENSITY)
          pred_store[[length(pred_store) + 1]] <- tibble(
            scenario = sc_name, trait = tr, model = "G", density = d,
            key = keys[te], observed = y[te], predicted = p)
      }
      job <- job + 1
    }
    message(sprintf("   %-12s %-20s done (%d accessions, %.0f s elapsed)",
                    sc_name, tr, length(keys),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

# ---- 4b. model comparison on the common set ----------------------------------
message("\n[2/4] genomic vs spectral vs taxonomy vs combined (common accession set)")
common <- intersect(gkeys, rownames(Sfull))
for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  if (sc$type != "within") next
  for (tr in TRAITS) {
    v <- year_value(tr, sc$year); if (!nrow(v)) next
    keys <- intersect(v$key, common)
    if (length(keys) < 40) next
    y <- setNames(v$v, v$key)[keys]
    KG <- Glist[[MODEL_DENSITY]][keys, keys, drop = FALSE]
    KS <- Sfull[keys, keys, drop = FALSE]
    KT <- Tfull[keys, keys, drop = FALSE]
    models <- list(G = list(KG), S = list(KS), T = list(KT), `G+S` = list(KG, KS))
    for (mn in names(models)) {
      for (des in CV_DESIGNS) {
        a <- run_cv(y, keys, models[[mn]], des,
                    paste(sc_name, tr, mn, MODEL_DENSITY, sep = "|"))
        if (nrow(a)) acc[[length(acc) + 1]] <- a |>
          mutate(scenario = sc_name, trait = tr, model = mn, density = MODEL_DENSITY)
      }
    }
    message(sprintf("   %-12s %-20s done (%d accessions, %.0f s elapsed)",
                    sc_name, tr, length(keys),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

accuracy <- bind_rows(acc) |>
  mutate(r = if_else(design == "loso", r_within, r_pooled)) |>
  select(scenario, trait, model, density, design, rep, r, r_pooled, r_within,
         n_folds_scored, n_folds, slope, rmse, n)
write_csv(accuracy, file.path(DIR_PRED, "prediction_accuracy.csv"))
if (length(pred_store))
  write_csv(bind_rows(pred_store), file.path(DIR_PRED, "prediction_predictions.csv"))

# ---- 4c. spectral preprocessing ----------------------------------------------
message("\n[3/4] spectral preprocessing comparison")
prep_res <- list()
for (meth in c("raw", "snv", "snv_d1", "snv_d2")) {
  Xp <- preprocess_spectra(Xspec, meth)
  Sp <- kernel_S(Xp); dimnames(Sp) <- dimnames(Sfull)
  for (tr in TRAITS) {
    v <- year_value(tr, "2025"); if (!nrow(v)) next
    keys <- intersect(v$key, rownames(Sp)); if (length(keys) < 40) next
    y <- setNames(v$v, v$key)[keys]
    a <- run_cv(y, keys, list(Sp[keys, keys, drop = FALSE]), "random",
                paste("prep", meth, tr, sep = "|"))
    if (nrow(a)) prep_res[[length(prep_res) + 1]] <-
      a |> mutate(preprocessing = meth, trait = tr)
  }
  message("   ", meth, " done (", sprintf("%.0f", as.numeric(difftime(Sys.time(), t0, units="secs"))), " s)")
}
if (length(prep_res))
  write_csv(bind_rows(prep_res) |>
              mutate(r = r_pooled) |>
              select(preprocessing, trait, design, rep, r, slope, rmse, n),
            file.path(DIR_PRED, "spectral_preprocessing.csv"))

# ---- 4d. learning curve ------------------------------------------------------
message("\n[4/4] training-set size learning curve")
lc <- list()
for (tr in TRAITS) {
  v <- year_value(tr, "2025"); if (!nrow(v)) next
  keys <- intersect(v$key, gkeys); if (length(keys) < 80) next
  y <- setNames(v$v, v$key)[keys]
  K <- Glist[[MODEL_DENSITY]][keys, keys, drop = FALSE]
  sizes <- unique(pmin(c(50, 100, 150, 200, 250), floor(length(keys) * 0.8)))
  for (ntr in sizes) {
    rs <- replicate(max(5, NREP %/% 2), {
      set.seed(NULL)
      tr_i <- sample(seq_along(y), ntr)
      te_i <- setdiff(seq_along(y), tr_i)
      p <- gblup_1k(y, K, tr_i[is.finite(y[tr_i])], te_i)
      pred_metrics(y[te_i], p)[["r"]]
    })
    lc[[length(lc) + 1]] <- tibble(trait = tr, n_train = ntr,
                                   r_mean = mean(rs, na.rm = TRUE),
                                   r_sd = sd(rs, na.rm = TRUE))
  }
  message("   ", tr, " done")
}
if (length(lc)) write_csv(bind_rows(lc), file.path(DIR_PRED, "learning_curve.csv"))

# ---- 5. PCA for the figure ---------------------------------------------------
# prcomp cannot take missing calls, so dosages are imputed to the marker mean --
# the same convention kernel_G uses, so the PCA and the kernel see one dataset.
Zpc <- Zlist[[MODEL_DENSITY]]
mu_col <- colMeans(Zpc, na.rm = TRUE)
na_ix <- which(is.na(Zpc), arr.ind = TRUE)
if (nrow(na_ix)) Zpc[na_ix] <- mu_col[na_ix[, "col"]]
Zpc <- Zpc[, apply(Zpc, 2, \(v) all(is.finite(v)) && sd(v) > 0), drop = FALSE]
message("PCA on ", ncol(Zpc), " markers (", nrow(na_ix), " calls imputed to the marker mean)")
pc_g <- prcomp(Zpc, scale. = FALSE)
pc_s <- prcomp(Xspec, scale. = TRUE)
kernels_pca <- bind_rows(
  tibble(source = "markers", key = gkeys,
         PC1 = pc_g$x[, 1], PC2 = pc_g$x[, 2],
         var1 = summary(pc_g)$importance[2, 1], var2 = summary(pc_g)$importance[2, 2]),
  tibble(source = "spectra", key = rownames(Xspec),
         PC1 = pc_s$x[, 1], PC2 = pc_s$x[, 2],
         var1 = summary(pc_s)$importance[2, 1], var2 = summary(pc_s)$importance[2, 2])) |>
  left_join(panel |> select(key, genotype, species), by = "key")
write_csv(kernels_pca, file.path(DIR_PRED, "kernels_pca.csv"))

# ---- 6. summary --------------------------------------------------------------
el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat("\n", strrep("=", 74), "\n", sep = "")
cat("PREDICTION RESULTS   (", sprintf("%.1f", el), " min)\n", sep = "")
cat(strrep("=", 74), "\n")
summ <- accuracy |>
  group_by(scenario, trait, model, density, design) |>
  summarise(r = mean(r, na.rm = TRUE), r_sd = sd(r, na.rm = TRUE),
            slope = mean(slope, na.rm = TRUE),
            folds_scored = mean(n_folds_scored), folds = mean(n_folds),
            n = max(n), .groups = "drop") |>
  mutate(r = if_else(is.nan(r), NA_real_, r))
cat("\nGenomic model, random CV, by marker density (mean r):\n")
print(as.data.frame(summ |> filter(model == "G", design == "random") |>
  select(scenario, trait, density, r) |>
  pivot_wider(names_from = density, values_from = r) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))), row.names = FALSE)

cat("\nModel comparison at ", MODEL_DENSITY, " markers (mean r, random CV):\n", sep = "")
cmp <- summ |> filter(density == MODEL_DENSITY, design == "random") |>
  select(scenario, trait, model, r) |>
  pivot_wider(names_from = model, values_from = r) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
print(as.data.frame(cmp), row.names = FALSE)

if ("loso" %in% CV_DESIGNS) {
  cat("\nRandom CV (pooled r) vs leave-one-species-out (mean within-species r):\n")
  gap <- summ |> filter(density == MODEL_DENSITY, design %in% c("random", "loso")) |>
    select(scenario, trait, model, design, r) |>
    pivot_wider(names_from = design, values_from = r) |>
    mutate(drop = random - loso, across(where(is.numeric), \(x) round(x, 3)))
  print(as.data.frame(gap), row.names = FALSE)
}

cat("\nAcross-year (train 2024, validate on accessions unique to the later year):\n")
ay <- summ |> filter(design == "holdout") |>
  select(scenario, trait, density, r, n) |>
  pivot_wider(names_from = density, values_from = r) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
if (nrow(ay)) print(as.data.frame(ay), row.names = FALSE) else
  cat("  (no across-year rows)\n")

write_csv(summ, file.path(DIR_PRED, "prediction_summary.csv"))
message("\nwritten to ", DIR_PRED)

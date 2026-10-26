# ==============================================================================
# 03_phenotype_models.R
#
# Mixed-model analysis of the adventitious-rooting traits: BLUPs, variance
# components, repeatability / heritability, and the transplant-response index.
#
# 00_raw_data/BLUPs_by_year_and_trait.csv holds ONE value per accession per year
# per trait. Despite the file name these are not BLUPs, and they are not measured
# the same way in every year:
#
#   trait            2024                    2025                   2026
#   ---------------- ----------------------- ---------------------- --------------
#   Callus-Stage     one score per bundle    mean of the individual one score per
#                    (integer 0-4)           cuttings (<= 12)       bundle (0-4)
#   Post-transplant  mean of 6 plants        mean of 6 plants       one score per
#                                                                   bundle (0-4)
#   Cutting weight   mean of the cuttings    mean of the cuttings   not measured
#   Root weight      mean of 6 plants        mean of 6 plants       not measured
#
# Outputs (02_clean_data/)
#   pheno_long.csv               accession x year x trait, raw and standardised
#   blups_by_trait.csv           BLUP and BLUE per accession per trait
#   variance_components.csv      sigma2_a, sigma2_e, LRT for sigma2_a > 0
#   repeatability.csv            single-year, mean-basis and generalised, with CIs
#   repeatability_sensitivity.csv same, refitted without 2026
#   year_covariance_structure.csv observed var/cov/cor by year and year pair
#   variance_partition.csv       between- vs within-species accession variance
#   trait_correlations.csv       Pearson and Spearman among trait BLUPs
#   year_stability.csv           across-year correlations per trait
#   coverage_by_trait.csv        accessions evaluated in 1 / 2 / 3 years
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/03_phenotype_models.R"
  }))), "00_setup.R"))

for (p in c("lme4", "purrr")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("install.packages(\"", p, "\")")
}
suppressPackageStartupMessages({library(lme4); library(purrr)})

N_BOOT <- as.integer(Sys.getenv("ROOTING_NBOOT", unset = "1000"))

# ---- 1. Long-format phenotypes ----------------------------------------------
raw <- read_csv(file.path(DIR_RAW, "BLUPs_by_year_and_trait.csv"),
                show_col_types = FALSE)

MEAS_BASIS <- tibble::tribble(
  ~trait,             ~year,  ~basis,                        ~n_sub,
  "Callus-Stage",     "2024", "single score per bundle",      1L,
  "Callus-Stage",     "2025", "mean of individual cuttings", 12L,
  "Callus-Stage",     "2026", "single score per bundle",      1L,
  "Post-transplant",  "2024", "mean of 6 plants",             6L,
  "Post-transplant",  "2025", "mean of 6 plants",             6L,
  "Post-transplant",  "2026", "single score per bundle",      1L,
  "Cutting weight",   "2024", "mean of the cuttings",        12L,
  "Cutting weight",   "2025", "mean of the cuttings",        12L,
  "Root weight",      "2024", "mean of 6 plants",             6L,
  "Root weight",      "2025", "mean of 6 plants",             6L
)

TRAIT_LEVELS <- c("Callus-Stage", "Post-transplant", "Root weight",
                  "Cutting weight", "Transplant response")
TRAIT_LABELS <- c("Callus-Stage"        = "Callus-stage rooting",
                  "Post-transplant"     = "Post-transplant rooting",
                  "Root weight"         = "Root weight",
                  "Cutting weight"      = "Cutting weight",
                  "Transplant response" = "Transplant response")

pheno_raw <- raw |>
  pivot_longer(-c(genotype, species), names_to = "key", values_to = "value") |>
  separate_wider_delim(key, delim = "_", names = c("trait", "year")) |>
  filter(!is.na(value)) |>
  left_join(MEAS_BASIS, by = c("trait", "year"))
stopifnot(!any(is.na(pheno_raw$basis)))

# ---- 2. Observed year covariance structure (the diagnostic behind the model) --
year_cov_structure <- pheno_raw |>
  select(genotype, trait, year, value) |>
  pivot_wider(names_from = year, values_from = value) |>
  group_by(trait) |>
  group_modify(function(d, key) {
    yrs <- setdiff(names(d), "genotype")
    yrs <- yrs[vapply(yrs, \(y) sum(!is.na(d[[y]])) > 0, logical(1))]
    vars <- tibble(term = "variance", year_a = yrs, year_b = NA_character_,
                   n = vapply(yrs, \(y) sum(!is.na(d[[y]])), integer(1)),
                   value = vapply(yrs, \(y) var(d[[y]], na.rm = TRUE), numeric(1)))
    if (length(yrs) < 2) return(vars)
    gr <- t(combn(yrs, 2))
    cv <- purrr::map_dfr(seq_len(nrow(gr)), function(i) {
      a <- gr[i, 1]; b <- gr[i, 2]
      ok <- complete.cases(d[[a]], d[[b]])
      tibble(term = c("covariance", "correlation"), year_a = a, year_b = b,
             n = sum(ok),
             value = c(cov(d[[a]], d[[b]], use = "complete.obs"),
                       cor(d[[a]], d[[b]], use = "complete.obs")))
    })
    bind_rows(vars, cv)
  }) |>
  ungroup()

# Under compound symmetry sigma2_a is the mean covariance and the implied
# residual variance is var(year) - sigma2_a. A negative value proves the raw
# scale is unusable, which is exactly why the analysis is run on z-scores.
cs_check <- year_cov_structure |>
  group_by(trait) |>
  summarise(sigma2_a_implied = mean(value[term == "covariance"]),
            min_sigma2_e_implied = min(value[term == "variance"]) - sigma2_a_implied,
            cs_violated = min_sigma2_e_implied < 0, .groups = "drop")
write_csv(left_join(year_cov_structure, cs_check, by = "trait"),
          file.path(DIR_CLEAN, "year_covariance_structure.csv"))

message("compound-symmetry check on the raw scale:")
print(as.data.frame(cs_check |> mutate(across(where(is.numeric), \(x) round(x, 3)))),
      row.names = FALSE)

# ---- 3. Within-year standardisation -----------------------------------------
pheno <- pheno_raw |>
  group_by(trait, year) |>
  mutate(year_mean = mean(value), year_sd = sd(value),
         z = (value - year_mean) / year_sd) |>
  ungroup()

# Scale used to express standardised BLUPs back in trait units: the mean of the
# within-year means and SDs, i.e. an "average year".
trait_scale <- pheno |>
  distinct(trait, year, year_mean, year_sd) |>
  group_by(trait) |>
  summarise(mean_bar = mean(year_mean), sd_bar = mean(year_sd), .groups = "drop")

# ---- 4. Transplant-response index -------------------------------------------
# "Does this accession root better after transplanting than its callus-stage
# performance predicts?" = residual of post-transplant on callus-stage, fitted
# within year on standardised values. A ratio was rejected: callus-stage is 0 for
# many accessions and the ratio is badly skewed.
wide_z <- pheno |>
  select(genotype, species, year, trait, z) |>
  pivot_wider(names_from = trait, values_from = z)

index_fits <- wide_z |>
  filter(!is.na(`Callus-Stage`), !is.na(`Post-transplant`)) |>
  group_by(year) |>
  group_modify(function(d, key) {
    fit <- lm(`Post-transplant` ~ `Callus-Stage`, data = d)
    d$value <- unname(residuals(fit))
    d$slope <- unname(coef(fit)[2]); d$r2 <- summary(fit)$r.squared
    d
  }) |>
  ungroup()

message(sprintf("\ntransplant-response index, within-year regression:\n%s",
                paste(sprintf("  %s  slope=%.2f  R2=%.2f  n=%d",
                              unique(index_fits$year),
                              index_fits$slope[!duplicated(index_fits$year)],
                              index_fits$r2[!duplicated(index_fits$year)],
                              as.integer(table(index_fits$year))),
                      collapse = "\n")))

pheno <- bind_rows(
  pheno,
  index_fits |>
    transmute(genotype, species, trait = "Transplant response", year, value,
              basis = "residual of post-transplant on callus-stage (within year)",
              n_sub = NA_integer_) |>
    group_by(year) |>
    mutate(year_mean = mean(value), year_sd = sd(value),
           z = (value - year_mean) / year_sd) |>
    ungroup()
) |>
  mutate(trait     = factor(trait, levels = TRAIT_LEVELS),
         accession = factor(genotype),
         year      = factor(year),
         species   = factor(species, levels = SPECIES_LEVELS))

trait_scale <- bind_rows(
  trait_scale,
  tibble(trait = "Transplant response", mean_bar = 0,
         sd_bar = mean(tapply(index_fits$value, index_fits$year, sd)))
)

write_csv(pheno |> select(genotype, species, trait, year, value, z, basis, n_sub),
          file.path(DIR_CLEAN, "pheno_long.csv"))

# ---- 5. Coverage -------------------------------------------------------------
coverage <- pheno |>
  group_by(trait, genotype) |>
  summarise(n_years = dplyr::n_distinct(year),
            years   = paste(sort(unique(as.character(year))), collapse = "+"),
            .groups = "drop")

write_csv(
  bind_rows(
    coverage |> count(trait, n_years) |>
      transmute(trait, type = "n_years", level = as.character(n_years), n),
    coverage |> count(trait, years) |>
      transmute(trait, type = "combination", level = years, n)),
  file.path(DIR_CLEAN, "coverage_by_trait.csv"))

message("\naccessions by number of years evaluated")
print(as.data.frame(coverage |> count(trait, n_years) |>
                      pivot_wider(names_from = n_years, values_from = n,
                                  names_prefix = "yr_")), row.names = FALSE)

# ---- 6. Model fitting --------------------------------------------------------

fit_trait <- function(d, trait_name, n_boot = N_BOOT, label = "all years") {

  d <- d |> filter(!is.na(z)) |>
    mutate(accession = droplevels(accession), year = droplevels(year))
  if (nlevels(d$year) < 2 || nrow(d) < 20) return(NULL)

  # Year effects are already removed by the within-year standardisation, so the
  # fixed part is an intercept only.
  m <- lme4::lmer(z ~ 1 + (1 | accession), data = d, REML = TRUE)

  vc  <- as.data.frame(lme4::VarCorr(m))
  sa2 <- vc$vcov[vc$grp == "accession"]
  se2 <- vc$vcov[vc$grp == "Residual"]

  # Boundary LRT: the null sits on the edge of the parameter space, so the
  # statistic is a 50:50 mixture of chi2_0 and chi2_1 -> halve the p-value.
  lrt <- as.numeric(2 * (logLik(lme4::refitML(m)) - logLik(lm(z ~ 1, data = d))))
  p_a <- 0.5 * pchisq(max(lrt, 0), df = 1, lower.tail = FALSE)

  n_i    <- as.integer(table(d$accession))
  n_harm <- length(n_i) / sum(1 / n_i)

  rep_single <- sa2 / (sa2 + se2)                 # one year
  rep_mean   <- sa2 / (sa2 + se2 / n_harm)        # mean of the years available
  pev        <- as.vector(attr(lme4::ranef(m, condVar = TRUE)$accession, "postVar"))
  rep_gen    <- 1 - mean(pev) / sa2               # Cullis et al. (2006)

  ci <- rep(NA_real_, 6)
  if (n_boot > 0) {
    stat <- function(fit) {
      v  <- as.data.frame(lme4::VarCorr(fit))
      g  <- v$vcov[v$grp == "accession"]; e <- v$vcov[v$grp == "Residual"]
      pv <- as.vector(attr(lme4::ranef(fit, condVar = TRUE)$accession, "postVar"))
      c(g / (g + e), g / (g + e / n_harm), if (g > 0) 1 - mean(pv) / g else 0)
    }
    bb <- try(suppressWarnings(suppressMessages(
      lme4::bootMer(m, stat, nsim = n_boot, type = "parametric",
                    use.u = FALSE, .progress = "none"))), silent = TRUE)
    if (!inherits(bb, "try-error"))
      ci <- as.vector(vapply(1:3, \(k) quantile(bb$t[, k], c(0.025, 0.975),
                                                na.rm = TRUE), numeric(2)))
  }

  sc   <- trait_scale[trait_scale$trait == trait_name, ]
  re   <- lme4::ranef(m)$accession
  blup <- tibble(genotype = rownames(re),
                 trait = trait_name,
                 blup_z = re[, 1] + lme4::fixef(m)[[1]]) |>
    mutate(blup = blup_z * sc$sd_bar + sc$mean_bar)

  fb <- lm(z ~ 0 + accession, data = d)
  cf <- coef(fb)
  blue <- tibble(genotype = sub("^accession", "", names(cf)),
                 trait = trait_name, blue_z = unname(cf)) |>
    mutate(blue = blue_z * sc$sd_bar + sc$mean_bar)

  list(
    trait = trait_name,
    vc = tibble(trait = trait_name, subset = label,
                n_obs = nrow(d), n_accessions = nlevels(d$accession),
                n_years = nlevels(d$year), n_harmonic_years = n_harm,
                sigma2_a = sa2, sigma2_e = se2, lrt = lrt, p_sigma2_a = p_a),
    rep = tibble(trait = trait_name, subset = label,
                 rep_single_year = rep_single, rep_single_lo = ci[1], rep_single_hi = ci[2],
                 rep_entry_mean  = rep_mean,   rep_mean_lo   = ci[3], rep_mean_hi   = ci[4],
                 rep_generalised = rep_gen,    rep_gen_lo    = ci[5], rep_gen_hi    = ci[6],
                 n_boot = n_boot),
    blup = blup, blue = blue)
}

traits <- levels(pheno$trait)

fits <- purrr::compact(lapply(traits, function(tr) {
  message("fitting ", tr, " ...")
  fit_trait(pheno |> filter(trait == tr), tr)
}))

# Sensitivity: drop 2026, the year scored once per bundle.
fits_no26 <- purrr::compact(lapply(traits, function(tr) {
  d <- pheno |> filter(trait == tr, year != "2026")
  if (dplyr::n_distinct(d$year) < 2) return(NULL)
  fit_trait(d, tr, n_boot = 0, label = "2024 + 2025 only")
}))

variance_components <- bind_rows(lapply(fits, `[[`, "vc"))
repeatability       <- bind_rows(lapply(fits, `[[`, "rep"))
sens <- bind_rows(bind_rows(lapply(fits_no26, `[[`, "vc")) |>
                    select(trait, subset, n_obs, n_accessions, sigma2_a, sigma2_e),
                  .id = NULL) |>
  left_join(bind_rows(lapply(fits_no26, `[[`, "rep")) |>
              select(trait, rep_single_year, rep_entry_mean, rep_generalised),
            by = "trait")

blups <- bind_rows(lapply(fits, `[[`, "blup")) |>
  left_join(bind_rows(lapply(fits, `[[`, "blue")), by = c("genotype", "trait"))

write_csv(variance_components, file.path(DIR_CLEAN, "variance_components.csv"))
write_csv(repeatability,       file.path(DIR_CLEAN, "repeatability.csv"))
write_csv(sens,                file.path(DIR_CLEAN, "repeatability_sensitivity.csv"))

blups_wide <- blups |>
  mutate(trait = factor(trait, levels = TRAIT_LEVELS)) |>
  pivot_wider(names_from = trait,
              values_from = c(blup, blup_z, blue, blue_z),
              names_glue = "{trait}_{.value}") |>
  left_join(raw |> select(genotype, species), by = "genotype") |>
  relocate(genotype, species)
write_csv(blups_wide, file.path(DIR_CLEAN, "blups_by_trait.csv"))

# ---- 7. Between- vs within-species accession variance ------------------------
partition <- purrr::compact(lapply(traits, function(tr) {
  d <- pheno |> filter(trait == tr, !is.na(z)) |>
    mutate(accession = droplevels(accession), species = droplevels(species))
  if (nlevels(d$species) < 3) return(NULL)
  m <- try(suppressMessages(
    lme4::lmer(z ~ 1 + (1 | species) + (1 | species:accession), data = d)),
    silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  v <- as.data.frame(lme4::VarCorr(m))
  g <- function(n) { x <- v$vcov[v$grp == n]; if (length(x)) x else NA_real_ }
  tibble(trait = tr, sigma2_species = g("species"),
         sigma2_accession_within = g("species:accession"),
         sigma2_residual = g("Residual")) |>
    mutate(pct_between_species = 100 * sigma2_species /
             (sigma2_species + sigma2_accession_within))
})) |> bind_rows()
write_csv(partition, file.path(DIR_CLEAN, "variance_partition.csv"))

# ---- 8. Trait correlations (on BLUPs) ---------------------------------------
bmat <- blups |> select(genotype, trait, blup_z) |>
  pivot_wider(names_from = trait, values_from = blup_z)
keep <- intersect(TRAIT_LEVELS, names(bmat))
X <- as.matrix(bmat[, keep])

cor_long <- function(method) {
  m <- cor(X, use = "pairwise.complete.obs", method = method)
  n <- crossprod(!is.na(X))
  grid <- expand.grid(trait_x = rownames(m), trait_y = colnames(m),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) |>
    as_tibble() |>
    mutate(r = as.vector(m), n = as.vector(n), method = method)
  # p-values from the pairwise-complete test for each trait pair
  grid$p <- purrr::pmap_dbl(grid[c("trait_x", "trait_y")], function(trait_x, trait_y) {
    if (identical(trait_x, trait_y)) return(0)
    ok <- !is.na(X[, trait_x]) & !is.na(X[, trait_y])
    if (sum(ok) < 4) return(NA_real_)
    suppressWarnings(cor.test(X[ok, trait_x], X[ok, trait_y],
                              method = method)$p.value)
  })
  grid
}
trait_correlations <- bind_rows(cor_long("pearson"), cor_long("spearman")) |>
  mutate(across(c(trait_x, trait_y), \(x) factor(x, levels = TRAIT_LEVELS)))
write_csv(trait_correlations, file.path(DIR_CLEAN, "trait_correlations.csv"))

# ---- 9. Year-to-year stability ----------------------------------------------
year_stability <- pheno |>
  select(genotype, trait, year, z) |>
  pivot_wider(names_from = year, values_from = z) |>
  group_by(trait) |>
  group_modify(function(d, key) {
    yrs <- sort(setdiff(names(d), "genotype"))
    if (length(yrs) < 2) return(tibble())
    gr <- t(combn(yrs, 2))
    purrr::map_dfr(seq_len(nrow(gr)), function(i) {
      a <- gr[i, 1]; b <- gr[i, 2]
      ok <- !is.na(d[[a]]) & !is.na(d[[b]])
      if (sum(ok) < 5)
        return(tibble(year_a = a, year_b = b, n = sum(ok),
                      r = NA_real_, r_lo = NA_real_, r_hi = NA_real_,
                      rho = NA_real_, p = NA_real_))
      ct <- suppressWarnings(cor.test(d[[a]][ok], d[[b]][ok]))
      tibble(year_a = a, year_b = b, n = sum(ok),
             r = unname(ct$estimate),
             r_lo = ct$conf.int[1], r_hi = ct$conf.int[2],
             rho = suppressWarnings(cor(d[[a]][ok], d[[b]][ok], method = "spearman")),
             p = ct$p.value)
    })
  }) |>
  ungroup()
write_csv(year_stability, file.path(DIR_CLEAN, "year_stability.csv"))

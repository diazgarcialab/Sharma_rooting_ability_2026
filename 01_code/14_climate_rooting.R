# ==============================================================================
# 14_climate_rooting.R -- is rooting associated with the climate of the site an
# accession was collected from?
#
#   Rscript 01_code/14_climate_rooting.R
#
# Output (02_clean_data/)
#   climate_rooting_associations.csv   one row per trait x climate variable
#
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/14_climate_rooting.R"
  }))), "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(lme4)
})

MIN_N        <- 50L    # a trait-variable pair needs at least this many accessions
MIN_SPECIES  <- 3L     # ... spread over at least this many species

panel <- read_csv(file.path(DIR_CLEAN, "accessions_panel.csv"), show_col_types = FALSE)
blups <- read_csv(file.path(DIR_CLEAN, "blups_by_trait.csv"), show_col_types = FALSE)

CLIM <- c(grep("Temperature|Precipitation|Diurnal", names(panel), value = TRUE), "altitude")
CLIM <- setdiff(CLIM, c("clim_id", "n_recs", "colour", "shape"))

d <- panel |>
  filter(georeferenced) |>
  select(genotype, species, all_of(CLIM)) |>
  inner_join(blups |> select(genotype, ends_with("_blup")), by = "genotype")

TRAITS <- grep("_blup$", names(d), value = TRUE)
message(sprintf("georeferenced accessions with BLUPs: %d across %d species",
                nrow(d), n_distinct(d$species)))
message(sprintf("climate variables: %d | traits: %d", length(CLIM), length(TRAITS)))

res <- list()
for (tr in TRAITS) for (cv in CLIM) {
  ok <- is.finite(d[[cv]]) & is.finite(d[[tr]])
  dd <- d[ok, ]
  if (nrow(dd) < MIN_N || n_distinct(dd$species) < MIN_SPECIES) next
  y <- dd[[tr]]; x <- as.vector(scale(dd[[cv]]))
  m <- suppressMessages(try(lmer(y ~ x + (1 | species), data = dd, REML = TRUE),
                            silent = TRUE))
  if (inherits(m, "try-error")) next
  bta <- unname(fixef(m)[2]); se <- unname(sqrt(diag(vcov(m)))[2])
  # Species variance retained by the model, for context on how much of the panel
  # structure is being absorbed.
  vc <- as.data.frame(VarCorr(m))
  v_sp <- vc$vcov[vc$grp == "species"]; v_res <- vc$vcov[vc$grp == "Residual"]
  res[[length(res) + 1]] <- tibble(
    trait = sub("_blup$", "", tr), climate = cv,
    r_across = cor(dd[[cv]], y), beta_within = bta, se = se, t = bta / se,
    icc_species = v_sp / (v_sp + v_res), n = nrow(dd),
    n_species = n_distinct(dd$species))
}
out <- bind_rows(res) |>
  mutate(p_within = 2 * pnorm(-abs(t)),
         q_within = p.adjust(p_within, method = "BH")) |>
  arrange(p_within)

write_csv(out, file.path(DIR_CLEAN, "climate_rooting_associations.csv"))

cat("\n", strrep("=", 76), "\n", sep = "")
cat("CLIMATE OF ORIGIN AND ROOTING -- associations, species held as random\n")
cat(strrep("=", 76), "\n")
cat(sprintf("%d tests | %d with BH q < 0.05 | %d with q < 0.10\n\n",
            nrow(out), sum(out$q_within < 0.05), sum(out$q_within < 0.10)))
top <- out |> filter(q_within < 0.10) |>
  mutate(across(c(r_across, beta_within, icc_species), \(x) round(x, 3)),
         p_within = signif(p_within, 2), q_within = signif(q_within, 2)) |>
  select(trait, climate, r_across, beta_within, p_within, q_within, n)
if (nrow(top)) print(as.data.frame(top), row.names = FALSE) else
  cat("  nothing reached q < 0.10\n")

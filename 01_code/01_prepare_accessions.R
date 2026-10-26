# ==============================================================================
# 01_prepare_accessions.R
#
# Builds the canonical accession table used by every downstream figure:
#   * the accessions that have rooting-ability BLUPs (the study panel, n = 308)
#   * their species and geographic group
#   * their collection coordinates and bioclim variables, joined from
#     00_raw_data/climatic_data_allaccessions.csv
#
#
# Outputs (02_clean_data/)
#   accessions_panel.csv            one row per BLUP accession (n = 308)
#   accessions_georeferenced.csv    the subset carrying coordinates
#   id_match_report.csv             how every ID was resolved
#   species_summary.csv             per-species counts, panel vs georeferenced
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/01_prepare_accessions.R"
  }))), "00_setup.R"))

# ---- 1. Read raw inputs ------------------------------------------------------
blups <- read_csv(file.path(DIR_RAW, "BLUPs_by_year_and_trait.csv"),
                  show_col_types = FALSE)
clim  <- read_csv(file.path(DIR_RAW, "climatic_data_allaccessions.csv"),
                  show_col_types = FALSE)

stopifnot(all(c("genotype", "species") %in% names(blups)))
stopifnot(all(c("genotype", "latitude", "longitude") %in% names(clim)))

message("BLUP panel      : ", nrow(blups), " accessions, ",
        dplyr::n_distinct(blups$species), " species")
message("climate table   : ", nrow(clim), " georeferenced records")

# Trait columns are everything except the two ID columns; kept so the same table
# can feed the phenotype figures.
trait_cols <- setdiff(names(blups), c("genotype", "species"))

# ---- 2. ID harmonisation -----------------------------------------------------
norm_id <- function(x) {
  x |>
    stringr::str_trim() |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]", "")
}

clim <- clim |>
  mutate(clim_id = genotype, key = norm_id(genotype)) |>
  select(-genotype)

# Climate records collapsed to one row per normalised key (identical keys share a
# site; take the mean, which is a no-op when the coordinates are identical).
clim_by_key <- clim |>
  group_by(key) |>
  summarise(across(where(is.numeric), \(v) mean(v, na.rm = TRUE)),
            clim_id  = paste(unique(clim_id), collapse = "; "),
            n_recs   = dplyr::n(),
            .groups  = "drop")

panel <- blups |>
  mutate(key = norm_id(genotype))

# -- pass 1 + 2: exact / normalised -------------------------------------------
exact_ids <- clim$clim_id
panel <- panel |>
  mutate(match_type = dplyr::case_when(
    genotype %in% exact_ids ~ "exact",
    key %in% clim_by_key$key ~ "normalised",
    TRUE ~ NA_character_
  ))

# -- pass 3: site match (one BLUP id -> several sub-accessions, one coordinate) -
COORD_TOL <- 0.01   # degrees (~1 km); refuse the join if the site is not tight

site_match <- function(k) {
  hits <- clim[startsWith(clim$key, k) & clim$key != k, , drop = FALSE]
  if (!nrow(hits)) return(NULL)
  if (diff(range(hits$latitude,  na.rm = TRUE)) > COORD_TOL ||
      diff(range(hits$longitude, na.rm = TRUE)) > COORD_TOL) return(NULL)
  hits
}

unresolved <- which(is.na(panel$match_type))
site_rows  <- list()
for (i in unresolved) {
  hits <- site_match(panel$key[i])
  if (is.null(hits)) next
  panel$match_type[i] <- "site"
  site_rows[[panel$key[i]]] <- hits |>
    summarise(across(where(is.numeric), \(v) mean(v, na.rm = TRUE)),
              clim_id = paste0(dplyr::first(clim_id), " ... ", dplyr::last(clim_id)),
              n_recs  = dplyr::n()) |>
    mutate(key = panel$key[i])
}
site_tbl <- if (length(site_rows)) dplyr::bind_rows(site_rows) else clim_by_key[0, ]

panel$match_type[is.na(panel$match_type)] <- "no_coordinates"

# ---- 3. Join climate ---------------------------------------------------------
lookup <- dplyr::bind_rows(clim_by_key, site_tbl) |>
  distinct(key, .keep_all = TRUE)

accessions <- panel |>
  left_join(lookup, by = "key") |>
  left_join(SPECIES_LUT |> select(species, label, group, colour, shape),
            by = "species") |>
  mutate(
    species     = factor(species, levels = SPECIES_LEVELS),
    georeferenced = !is.na(latitude) & !is.na(longitude)
  ) |>
  relocate(genotype, species, label, group, match_type, georeferenced,
           latitude, longitude, altitude)

# Fail loudly if the BLUP file gains a species that is not in the look-up table.
unknown_sp <- setdiff(as.character(blups$species), SPECIES_LEVELS)
if (length(unknown_sp)) {
  stop("Species missing from SPECIES_LUT in 00_setup.R: ",
       paste(unknown_sp, collapse = ", "))
}

# ---- 4. Report ---------------------------------------------------------------
tab <- table(accessions$match_type)
message("\nID matching")
for (nm in names(tab)) message(sprintf("  %-16s %4d", nm, tab[[nm]]))
message(sprintf("  %-16s %4d / %d (%.0f%%)", "georeferenced",
                sum(accessions$georeferenced), nrow(accessions),
                100 * mean(accessions$georeferenced)))

match_report <- accessions |>
  select(genotype, species, match_type, matched_climate_id = clim_id,
         n_climate_records = n_recs, latitude, longitude) |>
  arrange(match_type, genotype)

species_summary <- accessions |>
  count(species, group, name = "n_panel") |>
  left_join(accessions |> filter(georeferenced) |> count(species, name = "n_mapped"),
            by = "species") |>
  mutate(n_mapped = tidyr::replace_na(n_mapped, 0)) |>
  arrange(group, species)

# ---- 5. Write ----------------------------------------------------------------
write_csv(accessions, file.path(DIR_CLEAN, "accessions_panel.csv"))
write_csv(filter(accessions, georeferenced),
          file.path(DIR_CLEAN, "accessions_georeferenced.csv"))
write_csv(match_report,    file.path(DIR_CLEAN, "id_match_report.csv"))
write_csv(species_summary, file.path(DIR_CLEAN, "species_summary.csv"))

message("\nwrote 4 files to ", DIR_CLEAN)
print(as.data.frame(species_summary), row.names = FALSE)

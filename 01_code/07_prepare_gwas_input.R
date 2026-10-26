# ==============================================================================
# 07_prepare_gwas_input.R
#
# Builds the PLINK / GEMMA input files for the rooting GWAS from the sequencing
# manifest 00_raw_data/Wild_Vitis_sequenced_accessions_May22.xlsx.
#
# Outputs (02_clean_data/gwas/)
#   accession_libraries.csv        every accession x library, with plate details
#   fam_phenotypes_all.txt         every library, input to the resolve step
#   individuals_blups.txt          FID IID, provisional (first library)
#   sample_ids.txt                 one ID per line, provisional
#   filtered_phenotypes_blups.txt  FID IID + one column per trait, provisional
#   fam_phenotypes.txt             FID IID dummy + traits, provisional
#   pheno_<trait>.txt              FID IID value, one file per trait
#   trait_index.tsv                which GEMMA -n index is which trait
#   gwas_id_map.csv                accession -> library audit
#   manifest_conflicts.csv         written only if the manifest is ambiguous
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/07_prepare_gwas_input.R"
  }))), "00_setup.R"))

if (!requireNamespace("readxl", quietly = TRUE))
  stop('install.packages("readxl")   # needed to read the sequencing manifest')

PHENO_TYPE <- Sys.getenv("ROOTING_PHENO_TYPE", unset = "blup")
stopifnot(PHENO_TYPE %in% c("blup", "blue", "blup_z", "blue_z"))

MANIFEST <- Sys.getenv(
  "ROOTING_MANIFEST",
  unset = file.path(DIR_RAW, "Wild_Vitis_sequenced_accessions_May22.xlsx"))
if (!file.exists(MANIFEST)) stop("sequencing manifest not found: ", MANIFEST)

DIR_GWAS <- file.path(DIR_CLEAN, "gwas")
dir.create(DIR_GWAS, showWarnings = FALSE, recursive = TRUE)

TRAIT_LEVELS <- c("Callus-Stage", "Post-transplant", "Root weight",
                  "Cutting weight", "Transplant response")
TRAIT_SLUG   <- c("Callus-Stage"        = "callus_stage",
                  "Post-transplant"     = "post_transplant",
                  "Root weight"         = "root_weight",
                  "Cutting weight"      = "cutting_weight",
                  "Transplant response" = "transplant_response")

norm_id <- function(x) {
  x |> as.character() |> stringr::str_trim() |> stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]", "")
}

# ---- 1. Inputs ---------------------------------------------------------------
blup_file <- file.path(DIR_CLEAN, "blups_by_trait.csv")
if (!file.exists(blup_file)) source(file.path(DIR_CODE, "03_phenotype_models.R"))
blups <- read_csv(blup_file, show_col_types = FALSE)

raw_man <- readxl::read_excel(MANIFEST, sheet = 1, .name_repair = "minimal")
names(raw_man) <- stringr::str_trim(names(raw_man))

need <- c("Genotype", "Filename-READ1")
if (length(setdiff(need, names(raw_man))))
  stop("manifest is missing columns: ",
       paste(setdiff(need, names(raw_man)), collapse = ", "))

pick <- function(d, col) if (col %in% names(d)) as.character(d[[col]]) else NA_character_

man <- tibble(
  accession        = stringr::str_trim(as.character(raw_man$Genotype)),
  read1            = stringr::str_trim(as.character(raw_man$`Filename-READ1`)),
  manifest_species = pick(raw_man, "Species1"),
  plate            = pick(raw_man, "Plate"),
  well             = pick(raw_man, "Well #"),
  duplicated_flag  = pick(raw_man, "Duplicated"),
  sequenced        = pick(raw_man, "Sequenced")) |>
  filter(!is.na(accession), accession != "", !is.na(read1), read1 != "") |>
  mutate(geno_id = sub("-REA.*$", "", read1),
         key     = norm_id(accession)) |>
  select(-read1)

message(sprintf("manifest: %d libraries, %d distinct accessions",
                nrow(man), dplyr::n_distinct(man$key)))

# ---- 2. Manifest conflicts ---------------------------------------------------
# A library ID must belong to exactly one accession. Report any that do not, and
# keep the first row so the run can proceed deterministically.
conflicts <- man |>
  group_by(geno_id) |>
  filter(dplyr::n_distinct(accession) > 1) |>
  ungroup() |>
  select(geno_id, accession, manifest_species, plate, well, sequenced) |>
  arrange(geno_id, accession)

if (nrow(conflicts)) {
  write_csv(conflicts, file.path(DIR_GWAS, "manifest_conflicts.csv"))
  message("\nlibrary IDs claimed by more than one accession (first row kept):")
  print(as.data.frame(conflicts), row.names = FALSE)
  warning(sprintf("%d ambiguous library ID(s); see manifest_conflicts.csv",
                  dplyr::n_distinct(conflicts$geno_id)), call. = FALSE)
  man <- man |> group_by(geno_id) |> slice(1) |> ungroup()
} else {
  unlink(file.path(DIR_GWAS, "manifest_conflicts.csv"))
}

# ---- 3. Accession -> libraries ----------------------------------------------
libs <- blups |>
  select(genotype, species) |>
  mutate(key = norm_id(genotype)) |>
  inner_join(man, by = "key", relationship = "many-to-many") |>
  group_by(genotype) |>
  mutate(lib_index = dplyr::row_number(), n_libraries = dplyr::n()) |>
  ungroup() |>
  arrange(genotype, lib_index)

matched_acc <- dplyr::n_distinct(libs$genotype)
message(sprintf("\nphenotyped accessions found in the manifest: %d / %d",
                matched_acc, nrow(blups)))

unmatched <- setdiff(blups$genotype, libs$genotype)
if (length(unmatched)) {
  message("  NOT found: ", paste(unmatched, collapse = ", "))
  write_csv(tibble(genotype = unmatched),
            file.path(DIR_GWAS, "unmatched_for_review.csv"))
} else {
  unlink(file.path(DIR_GWAS, "unmatched_for_review.csv"))
  message("  every phenotyped accession has at least one sequencing library")
}

message("\nlibraries per accession:")
print(as.data.frame(libs |> distinct(genotype, n_libraries) |> count(n_libraries)),
      row.names = FALSE)


sp_check <- libs |>
  # strip the "V. " prefix BEFORE lowercasing, and drop the hybrid marker so
  # "V. x doaniana" and "doaniana" are treated as the same name
  mutate(blup_sp = tolower(sub("^V\\.\\s*(x|\u00d7)?\\s*", "",
                              stringr::str_trim(as.character(species)))),
         man_sp  = tolower(sub("^(x|\u00d7)\\s*", "",
                              stringr::str_trim(manifest_species)))) |>
  filter(!is.na(man_sp), man_sp != "", blup_sp != man_sp) |>
  distinct(genotype, species, manifest_species)
if (nrow(sp_check)) {
  message(sprintf("\nspecies disagreements between the BLUP file and the manifest: %d",
                  nrow(sp_check)))
  print(as.data.frame(head(sp_check, 20)), row.names = FALSE)
  write_csv(sp_check, file.path(DIR_GWAS, "species_disagreements.csv"))
} else {
  unlink(file.path(DIR_GWAS, "species_disagreements.csv"))
}

write_csv(libs |> select(genotype, species, geno_id, lib_index, n_libraries,
                         manifest_species, plate, well, duplicated_flag, sequenced),
          file.path(DIR_GWAS, "accession_libraries.csv"))
write_csv(libs |> transmute(genotype, species, geno_id, lib_index, n_libraries),
          file.path(DIR_GWAS, "gwas_id_map.csv"))


writeLines(sort(unique(man$geno_id[nzchar(man$geno_id)])),
           file.path(DIR_GWAS, "manifest_all_libraries.txt"))

# ---- 4. Phenotype values -----------------------------------------------------
value_cols <- paste0(TRAIT_LEVELS, "_", PHENO_TYPE)
if (length(setdiff(value_cols, names(blups))))
  stop("blups_by_trait.csv is missing: ",
       paste(setdiff(value_cols, names(blups)), collapse = ", "))

pheno <- blups |> select(genotype, all_of(value_cols)) |>
  rename_with(\(x) sub(paste0("_", PHENO_TYPE, "$"), "", x), all_of(value_cols))

fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 6))

all_tbl <- libs |>
  select(genotype, geno_id, lib_index) |>
  left_join(pheno, by = "genotype") |>
  arrange(genotype, lib_index)



if (any(grepl("\t", all_tbl$genotype))) stop("accession name contains a tab")
write.table(
  data.frame(all_tbl$geno_id, all_tbl$geno_id, all_tbl$genotype, dummy = "1",
             lapply(all_tbl[TRAIT_LEVELS], fmt), check.names = FALSE),
  file.path(DIR_GWAS, "fam_phenotypes_all.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- 5. Provisional single-library files ------------------------------------
one <- all_tbl |> filter(lib_index == 1)
stopifnot(!any(duplicated(one$genotype)), !any(duplicated(one$geno_id)))

write.table(data.frame(one$geno_id, one$geno_id),
            file.path(DIR_GWAS, "individuals_blups.txt"),
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(one$geno_id, file.path(DIR_GWAS, "sample_ids.txt"))
write.table(data.frame(one$geno_id, one$geno_id,
                       lapply(one[TRAIT_LEVELS], fmt), check.names = FALSE),
            file.path(DIR_GWAS, "filtered_phenotypes_blups.txt"),
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(one$geno_id, one$geno_id, dummy = "1",
                       lapply(one[TRAIT_LEVELS], fmt), check.names = FALSE),
            file.path(DIR_GWAS, "fam_phenotypes.txt"),
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

for (tr in TRAIT_LEVELS) {
  d <- data.frame(one$geno_id, one$geno_id, fmt(one[[tr]]))
  d <- d[!is.na(one[[tr]]), , drop = FALSE]
  write.table(d, file.path(DIR_GWAS, paste0("pheno_", TRAIT_SLUG[tr], ".txt")),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

trait_index <- tibble(
  n_index    = seq_along(TRAIT_LEVELS) + 1L,   # +1 for the dummy column
  trait      = TRAIT_LEVELS,
  slug       = unname(TRAIT_SLUG[TRAIT_LEVELS]),
  pheno_file = paste0("pheno_", TRAIT_SLUG[TRAIT_LEVELS], ".txt"),
  n_samples  = vapply(TRAIT_LEVELS, \(t) sum(!is.na(one[[t]])), integer(1)),
  value_type = PHENO_TYPE)
write.table(trait_index, file.path(DIR_GWAS, "trait_index.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write_csv(trait_index |> mutate(n_accessions = nrow(one)),
          file.path(DIR_GWAS, "gwas_coverage.csv"))

# ---- 5b. The same pipeline as plain commands --------------------------------

cmd <- c(
  sprintf("## Rooting-ability GWAS. Generated by 07_prepare_gwas_input.R on %s",
          format(Sys.Date())),
  sprintf("## Phenotype values: %s   |   traits: %d   |   accessions: %d",
          PHENO_TYPE, length(TRAIT_LEVELS), nrow(one)),
  "## Equivalent to: bash gemma_rooting.sh",
  "",
  "## 0. copy the inputs up from the laptop",
  "##   rsync -av 02_clean_data/gwas/ <user>@<server>:/data6/diazgarcia/rooting/",
  "cd /data6/diazgarcia/rooting",
  "",
  "## 1. RESOLVE: one sequencing library per accession.",
  sprintf("##    %d accessions hold %d libraries; more than one library of the same",
          nrow(one), nrow(all_tbl)),
  "##    plant would enter the GWAS as two individuals. Keep the first library",
  "##    of each accession that is actually present in the VCF.",
  "##    fam_phenotypes_all.txt is TAB separated: 38 accession names contain a space.",
  "bcftools query -l $VCF > vcf_samples.txt",
  paste0("awk -F'\\t' -v missing=accessions_not_in_vcf.txt '",
         "NR==FNR{inv[$1]=1; next} ",
         "{acc=$3; seen[acc]=1; ",
         "if (!(acc in taken) && ($1 in inv)) {taken[acc]=1; ",
         "printf \"%s %s\", $1, $2; for(i=4;i<=NF;i++) printf \" %s\", $i; printf \"\\n\"}} ",
         "END{for (a in seen) if (!(a in taken)) print a > missing}' ",
         "vcf_samples.txt fam_phenotypes_all.txt > fam_phenotypes.txt"),
  "awk '{print $1, $2}' fam_phenotypes.txt > individuals_blups.txt",
  "wc -l individuals_blups.txt   # <- the n to report in the paper",
  "",
  "## 2. genotype file, restricted to those samples",
  sprintf("plink1.9 --vcf $VCF --keep individuals_blups.txt --make-bed --out geno --threads 28 --allow-extra-chr --maf %s",
          "0.05"),
  "",
  "## 3. write every trait into the .fam as phenotype columns, keeping row order.",
  "##    GEMMA pairs the .fam, the kinship and the covariates by ROW, not by ID.",
  "cp geno.fam geno.fam.orig",
  sprintf(paste0("awk -v nt=%d 'NR==FNR{s=\"\"; for(i=3;i<=NF;i++) s=s\" \"$i; ph[$2]=s; next} ",
                 "{printf \"%%s %%s %%s %%s %%s\", $1,$2,$3,$4,$5; if ($2 in ph) printf \"%%s\", ph[$2]; ",
                 "else for(i=1;i<=nt;i++) printf \" NA\"; printf \"\\n\"}' ",
                 "fam_phenotypes.txt geno.fam.orig > geno.fam"),
          length(TRAIT_LEVELS) + 1L),
  "",
  "## 4. kinship. -n 1 is the complete dummy column, so every sample contributes",
  "gemma -bfile geno -gk 2 -n 1 -o kinship",
  "",
  "## 5. principal components",
  "plink1.9 --bfile geno --pca 10 --out geno_pca --allow-extra-chr",
  "",
  "## 6. covariates: a column of 1s for the intercept, then the 10 PCs, in .fam order.",
  "##    NOTE: GEMMA does not add an intercept itself. The chloride script passed",
  "##    the ten PCs alone; the leading 1 below is the fix.",
  paste0("awk 'NR==FNR{s=\"\"; for(i=3;i<=NF;i++) s=s\" \"$i; pc[$2]=s; next} ",
         "{print \"1\"pc[$2]}' geno_pca.eigenvec geno.fam > covariates_numeric.txt"),
  "",
  "## 7. association tests, one per trait",
  unlist(lapply(seq_len(nrow(trait_index)), \(i) c(
    sprintf("#  %s   (n = %d before the VCF intersection)",
            trait_index$trait[i], trait_index$n_samples[i]),
    sprintf("gemma -bfile geno -k output/kinship.sXX.txt -lmm 1 -c covariates_numeric.txt -n %d -o gwas_%s_%s_PCA10",
            trait_index$n_index[i], trait_index$slug[i], PHENO_TYPE),
    ""))),
  "## results:",
  sprintf("##   /data6/diazgarcia/rooting/output/gwas_<trait>_%s_PCA10.assoc.txt", PHENO_TYPE))
writeLines(cmd, file.path(DIR_GWAS, "gemma_commands_rooting.txt"))

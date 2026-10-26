# ==============================================================================
# 16_revision_figures.R -- supplementary figures added during revision of
# TPG-2026-08-0296 (Sharma et al., The Plant Genome).
#
#   Rscript 01_code/16_revision_figures.R
#
# Produces, in 03_figures/
#   FigureS1.{pdf,png,tiff}              experimental design matrix   (R1)
#   FigureS3.{...}       per-species climate slopes   (R1 line 216)
#   FigureS4.{...}                     LD decay curve               (R1 255-257, R2 L255)
#   FigureS5.{...}       root weight BLUPs by species (R1 318-323)
#   FigureS6.{...}                trait scatter plots          (R1 line 335)
# and, in 02_clean_data/
#   ld_decay.csv                                 the binned decay table
#
# ==============================================================================

source(file.path(dirname(normalizePath(local({
  for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) f[1] else "01_code/16_supplemental_figures.R"
}))), "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(stringr); library(patchwork)
})

ital <- function(x) SPECIES_LABELS[x]

# ------------------------------------------------------------------------------
# S1. Experimental design: accessions (rows, grouped by species) x season
# ------------------------------------------------------------------------------
pheno <- read_csv(file.path(DIR_CLEAN, "pheno_long.csv"), show_col_types = FALSE)

TRAIT_KEEP <- c("Callus-Stage", "Post-transplant", "Root weight", "Cutting weight")

design <- pheno |>
  filter(trait %in% TRAIT_KEEP) |>
  distinct(genotype, species, trait, year) |>
  mutate(trait = factor(trait, levels = TRAIT_KEEP),
         species = factor(species, levels = SPECIES_LEVELS))

# order accessions by species, then by how many seasons they contribute
ord <- design |>
  count(genotype, species, name = "n_cells") |>
  arrange(species, desc(n_cells), genotype)
design <- design |> mutate(genotype = factor(genotype, levels = ord$genotype))

# species band labels, placed at the midpoint of each species block
bands <- ord |>
  mutate(idx = row_number()) |>
  group_by(species) |>
  summarise(mid = mean(idx), lo = min(idx) - 0.5, hi = max(idx) + 0.5, n = n(),
            .groups = "drop")

# Label each species block at its midpoint. Matching labels to breaks BY NAME is
# what makes this robust: positional pairing silently falls back to the
# accession IDs if the two orderings ever diverge.
y_brk <- as.character(ord$genotype[round(bands$mid)])
y_lab <- SPECIES_LABELS[as.character(bands$species)]
y_lab[is.na(y_lab)] <- as.character(bands$species)[is.na(y_lab)]
y_lab <- setNames(as.character(y_lab), y_brk)

p_s1 <- ggplot(design, aes(x = factor(year), y = genotype, fill = species)) +
  geom_tile(width = 0.92, height = 0.92) +
  facet_wrap(~ trait, nrow = 1) +
  geom_hline(data = bands, aes(yintercept = hi), colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = SPECIES_COLS, labels = SPECIES_LABELS, name = NULL,
                    guide = guide_legend(ncol = 3, override.aes = list(size = 3))) +
  scale_y_discrete(breaks = y_brk, labels = y_lab) +
  labs(x = "Season", y = NULL) +
  theme_paper() +
  theme(axis.text.y = element_text(size = rel(0.8), face = "italic"),
        panel.grid  = element_blank(),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

save_figure(p_s1, "FigureS1", width = 7.0, height = 8.6)

# ------------------------------------------------------------------------------
# S3. Within-species slopes of each trait on the climate of the collection site
# ------------------------------------------------------------------------------
panel <- read_csv(file.path(DIR_CLEAN, "accessions_panel.csv"), show_col_types = FALSE)
blups <- read_csv(file.path(DIR_CLEAN, "blups_by_trait.csv"),  show_col_types = FALSE)
assoc <- read_csv(file.path(DIR_CLEAN, "climate_rooting_associations.csv"),
                  show_col_types = FALSE)

MIN_SP_N <- 10L                      # species need this many georeferenced accessions
sig <- assoc |> filter(q_within < 0.05) |> arrange(p_within) |> slice_head(n = 6)

d <- panel |>
  filter(georeferenced) |>
  select(genotype, species, any_of(unique(sig$climate))) |>
  inner_join(blups |> select(genotype, ends_with("_blup")), by = "genotype")

per_species <- purrr::map_dfr(seq_len(nrow(sig)), function(k) {
  tr <- paste0(sig$trait[k], "_blup"); cv <- sig$climate[k]
  dd <- d |> filter(is.finite(.data[[cv]]), is.finite(.data[[tr]]))
  dd$x <- as.vector(scale(dd[[cv]]))
  dd |>
    group_by(species) |>
    filter(n() >= MIN_SP_N) |>
    group_modify(~ {
      m <- stats::lm(.x[[tr]] ~ x, data = .x)
      s <- summary(m)$coefficients
      tibble(beta = s[2, 1], se = s[2, 2], n = nrow(.x))
    }) |>
    ungroup() |>
    mutate(trait = sig$trait[k], climate = cv,
           pooled = sig$beta_within[k], pooled_se = sig$se[k])
})

per_species <- per_species |>
  mutate(pair = paste0(trait, "\n", str_replace_all(climate, "_", " ")),
         species = factor(species, levels = rev(SPECIES_LEVELS)))

write_csv(per_species, file.path(DIR_CLEAN, "within_species_climate_slopes.csv"))

p_s3 <- ggplot(per_species, aes(x = beta, y = species)) +
  geom_vline(xintercept = 0, colour = INK[["secondary"]], linewidth = 0.3) +
  geom_rect(aes(xmin = pooled - 1.96 * pooled_se, xmax = pooled + 1.96 * pooled_se,
                ymin = -Inf, ymax = Inf), fill = PAL[["blue"]], alpha = 0.08,
            inherit.aes = FALSE,
            data = per_species |> distinct(pair, pooled, pooled_se)) +
  geom_vline(aes(xintercept = pooled), colour = PAL[["blue"]], linewidth = 0.4,
             linetype = "22",
             data = per_species |> distinct(pair, pooled)) +
  geom_errorbarh(aes(xmin = beta - 1.96 * se, xmax = beta + 1.96 * se),
                 height = 0, linewidth = 0.35, colour = INK[["secondary"]]) +
  geom_point(aes(size = n), colour = INK[["primary"]]) +
  facet_wrap(~ pair, scales = "free_x", nrow = 2) +
  scale_size_area(max_size = 2.6, name = "n") +
  scale_y_discrete(labels = SPECIES_LABELS) +
  labs(x = "Change in trait BLUP per SD of the climate variable", y = NULL) +
  theme_paper() +
  theme(axis.text.y = element_text(face = "italic"),
        strip.text = element_text(face = "bold", size = rel(0.9)),
        legend.position = "bottom")

save_figure(p_s3, "FigureS3", width = 7.2, height = 5.4)

# ------------------------------------------------------------------------------
# S4. Linkage disequilibrium decay
# ------------------------------------------------------------------------------
raw_path <- file.path(DIR_RAW, "genotype_matrix_50k.raw")
hdr  <- scan(raw_path, what = "", nlines = 1, quiet = TRUE)
snps <- hdr[-seq_len(6)]
G <- as.matrix(data.table::fread(raw_path, drop = 1:6, na.strings = "NA"))
storage.mode(G) <- "double"
# mean-impute, then standardise so that r is a simple cross-product
cm <- colMeans(G, na.rm = TRUE)
idx <- which(is.na(G), arr.ind = TRUE)
if (nrow(idx)) G[idx] <- cm[idx[, 2]]
sdv  <- apply(G, 2, stats::sd)
keep <- sdv > 0
G <- scale(G[, keep, drop = FALSE])
snps <- snps[keep]
chrom <- as.integer(sub("^chr(\\d+)_.*$", "\\1", snps))
pos   <- as.integer(sub("^chr\\d+_(\\d+)_.*$", "\\1", snps))
n_ind <- nrow(G)

BRK <- c(0, 1, 2, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200,
         300, 400, 500, 1000) * 1000
tot <- numeric(length(BRK) - 1); cnt <- numeric(length(BRK) - 1)
MAXLAG <- 300L          # markers ahead to consider; 1 Mb cap applied below

for (cc in sort(unique(chrom))) {
  m <- which(chrom == cc); o <- order(pos[m]); m <- m[o]
  P <- pos[m]; X <- G[, m, drop = FALSE]; k <- length(m)
  for (i in seq_len(k - 1L)) {
    j2 <- min(k, i + MAXLAG)
    if (j2 <= i) next
    dd <- P[(i + 1L):j2] - P[i]
    sel <- which(dd <= 1e6)
    if (!length(sel)) next
    r  <- as.vector(crossprod(X[, i], X[, (i + 1L):j2, drop = FALSE][, sel, drop = FALSE])) / n_ind
    b  <- .bincode(dd[sel], BRK, right = TRUE, include.lowest = TRUE)
    ok <- !is.na(b)
    tot <- tot + as.numeric(tapply(r[ok]^2, factor(b[ok], levels = seq_along(tot)),
                                   sum, default = 0))
    cnt <- cnt + tabulate(b[ok], nbins = length(cnt))
  }
}

set.seed(1)
ii <- sample.int(ncol(G), 60000, replace = TRUE)
jj <- sample.int(ncol(G), 60000, replace = TRUE)
ok <- chrom[ii] != chrom[jj]          # chrom/pos are already subset to `keep`
r_bg <- colSums(G[, ii[ok], drop = FALSE] * G[, jj[ok], drop = FALSE]) / n_ind
BG <- mean(r_bg^2)

ld <- tibble(lo = head(BRK, -1) / 1000, hi = tail(BRK, -1) / 1000,
             n_pairs = cnt, mean_r2 = tot / pmax(cnt, 1)) |>
  filter(n_pairs > 0) |>
  mutate(mid = sqrt(pmax(lo, 0.5) * hi))
write_csv(bind_rows(ld, tibble(lo = NA, hi = NA, n_pairs = sum(ok),
                               mean_r2 = BG, mid = NA)),
          file.path(DIR_CLEAN, "ld_decay.csv"))

p_s4 <- ggplot(ld, aes(mid, mean_r2)) +
  geom_hline(yintercept = BG, linetype = "22", colour = PAL[["vermilion"]],
             linewidth = 0.4) +
  annotate("text", x = max(ld$mid), y = BG, vjust = -0.6, hjust = 1,
           size = 2.4, colour = PAL[["vermilion"]],
           label = sprintf("inter-chromosomal baseline, r² = %.3f", BG)) +
  geom_line(colour = PAL[["blue"]], linewidth = 0.5) +
  geom_point(aes(size = n_pairs), colour = PAL[["blue"]]) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000),
                labels = c("1", "2", "5", "10", "25", "50", "100", "250", "500", "1000")) +
  scale_size_area(max_size = 2.8, name = "marker pairs", labels = scales::comma) +
  labs(x = "Distance between markers (kb, log scale)",
       y = expression("Mean "*r^2)) +
  theme_paper() + theme(legend.position = "bottom")

save_figure(p_s4, "FigureS4", width = 5.2, height = 3.4)

# ------------------------------------------------------------------------------
# S5. Root weight BLUPs by species (companion to Figure 2C, 2D)
# ------------------------------------------------------------------------------
rw <- blups |>
  select(genotype, species, value = `Root weight_blup`) |>
  filter(is.finite(value)) |>
  mutate(species = factor(species, levels = SPECIES_LEVELS))

ordr <- rw |> group_by(species) |> summarise(m = median(value), n = n(), .groups = "drop") |>
  arrange(m)
rw <- rw |> mutate(species = factor(species, levels = ordr$species))

p_s5 <- ggplot(rw, aes(x = value, y = species, fill = species)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.3, alpha = 0.75) +
  geom_point(position = position_jitter(height = 0.14, width = 0),
             size = 0.5, alpha = 0.6, colour = INK[["secondary"]]) +
  scale_fill_manual(values = SPECIES_COLS, guide = "none") +
  scale_y_discrete(labels = function(s) sprintf("%s (%d)", SPECIES_LABELS[s],
                                                ordr$n[match(s, ordr$species)])) +
  labs(x = "Root weight BLUP (g)", y = NULL) +
  theme_paper() + theme(axis.text.y = element_text(face = "italic"))

save_figure(p_s5, "FigureS5", width = 4.6, height = 4.2)

# ------------------------------------------------------------------------------
# S6. Pairwise scatter plots behind the correlation matrix of Figure 2E
# ------------------------------------------------------------------------------
TR <- c("Callus-Stage", "Post-transplant", "Root weight", "Cutting weight",
        "Transplant response")
wide <- blups |> select(genotype, species, all_of(paste0(TR, "_blup")))
names(wide) <- c("genotype", "species", TR)

pairs_df <- purrr::map_dfr(seq_along(TR), function(a)
  purrr::map_dfr(seq_along(TR), function(b) {
    if (b <= a) return(NULL)
    tibble(x = wide[[TR[a]]], y = wide[[TR[b]]], species = wide$species,
           xn = TR[a], yn = TR[b])
  })) |>
  filter(is.finite(x), is.finite(y)) |>
  mutate(pair = factor(paste0(yn, "  vs  ", xn)),
         species = factor(species, levels = SPECIES_LEVELS))

labs_r <- pairs_df |> group_by(pair) |>
  summarise(r = stats::cor(x, y), n = n(),
            x = min(x), y = max(y), .groups = "drop")

p_s6 <- ggplot(pairs_df, aes(x, y)) +
  geom_point(aes(colour = species), size = 0.55, alpha = 0.75) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = INK[["primary"]], linewidth = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = PAL[["vermilion"]], linetype = "22", linewidth = 0.4) +
  geom_text(data = labs_r, aes(x = x, y = y, label = sprintf("r = %.2f (n = %d)", r, n)),
            hjust = 0, vjust = 1, size = 2.2, colour = INK[["secondary"]]) +
  facet_wrap(~ pair, scales = "free", ncol = 4) +
  scale_colour_manual(values = SPECIES_COLS, labels = SPECIES_LABELS, name = NULL,
                      guide = guide_legend(nrow = 3, override.aes = list(size = 1.8))) +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom", legend.text = element_text(face = "italic"),
        strip.text = element_text(size = rel(0.85)))

save_figure(p_s6, "FigureS6", width = 7.2, height = 6.4)
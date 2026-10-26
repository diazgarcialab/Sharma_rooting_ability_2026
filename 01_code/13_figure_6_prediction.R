# ==============================================================================
# 13_Figure6.R -- Figure 6: the old Figures 6 and 7, fused.
#
#   Rscript 01_code/13_Figure6.R
#
# Requires 09_prepare_prediction_data.R and 11_run_predictions.R to have run.
#
# Panels
#   A  wavelength x wavelength correlation matrix        (old 6A)
#   B  trait-wavelength correlations                     (old 6B)
#   C  variance explained across wavelengths             (old 6C + 6D, one panel)
#   D  PCA on markers                                    (old 7A)
#   E  PCA on spectra                                    (old 7B)
#   F  genomic relationship matrix                       (old 7C)
#   G  spectral relationship matrix                      (old 7D)
#   H  taxonomy relationship matrix                      (old 7E)
#   I  accuracy by marker density and scenario           (old 7F)
#   J  model comparison: genomic / spectral / both / taxonomy   (old 7G)
#
# ==============================================================================

CODE_DIR <- dirname(normalizePath(local({
  for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) f[1] else "01_code/13_figure_6_prediction.R"
})))
source(file.path(CODE_DIR, "00_setup.R"))
source(file.path(CODE_DIR, "10_prediction_engine.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(patchwork); library(tibble)
})

DIR_PRED <- file.path(DIR_CLEAN, "prediction")
for (f in c("prediction_summary.csv", "kernels_pca.csv", "spectra_by_accession.csv"))
  if (!file.exists(file.path(DIR_PRED, f)))
    stop("missing ", f, " -- run 09 and 11 first")

# ---- styling, matching Figure 5 ---------------------------------------------
PT     <- 0.85
SIG    <- "#C1272D"
DEEP   <- "#2A4365"
PALE   <- "#90B4CE"
RAMP_CORR <- c("#1F3B70", "#4A93BF", "#AFD7DA", "#F4F2EF", "#FAC38C", "#DC6E4A", "#7F2733")
MODEL_LAB <- c(G = "Genomic", S = "Spectral", `G+S` = "Genomic + spectral",
               `T` = "Taxonomy")
MODEL_COL <- c(Genomic = DEEP, Spectral = "#E69F00",
               `Genomic + spectral` = SIG, Taxonomy = "#7A9CB5")
SCEN_LAB <- c(within_2024 = "Within 2024", within_2025 = "Within 2025",
              within_2026 = "Within 2026",
              across_24_25 = "Train 2024, test 2025",
              across_24_26 = "Train 2024, test 2026")
# Short forms for the crowded legends in B; the caption carries the full names.
TRAIT_TINY <- c("Callus-Stage" = "Callus-stage", "Post-transplant" = "Post-transplant",
                "Root weight" = "Root weight", "Cutting weight" = "Cutting weight",
                "Transplant response" = "Transplant response index")
TRAIT_SHORT <- c("Callus-Stage" = "Callus-stage rooting",
                 "Post-transplant" = "Post-transplant rooting",
                 "Root weight" = "Root weight",
                 "Cutting weight" = "Cutting weight",
                 "Transplant response" = "Transplant response index")

# ---- data --------------------------------------------------------------------
spec  <- read_csv(file.path(DIR_PRED, "spectra_by_accession.csv"), show_col_types = FALSE)
panel <- read_csv(file.path(DIR_PRED, "prediction_panel.csv"), show_col_types = FALSE)
summ  <- read_csv(file.path(DIR_PRED, "prediction_summary.csv"), show_col_types = FALSE)
pca   <- read_csv(file.path(DIR_PRED, "kernels_pca.csv"), show_col_types = FALSE)
blups <- read_csv(file.path(DIR_CLEAN, "blups_by_trait.csv"), show_col_types = FALSE)

wl_cols <- grep("^X[0-9]+$", names(spec), value = TRUE)
WL <- as.numeric(sub("^X", "", wl_cols))
X  <- as.matrix(spec[, wl_cols]); rownames(X) <- spec$key

# ---- A. wavelength x wavelength correlation ----------------------------------
STEP <- 5                                  # 721 x 721 cells is more than a panel can show
idx  <- seq(1, length(WL), by = STEP)
Cw   <- cor(X[, idx], use = "pairwise.complete.obs")
dA <- as_tibble(as.table(Cw), .name_repair = "minimal")
names(dA) <- c("i", "j", "r")
dA <- dA |> mutate(wi = as.numeric(sub("^X", "", i)), wj = as.numeric(sub("^X", "", j)))

lo_A <- floor(min(Cw, na.rm = TRUE) * 10) / 10
pA <- ggplot(dA, aes(wi, wj, fill = r)) +
  geom_raster() +
  scale_fill_gradientn(colours = RAMP_PRECIP, limits = c(lo_A, 1),
                       name = "Pearson r", breaks = c(lo_A, (lo_A + 1) / 2, 1),
                       labels = \(b) sprintf("%.1f", b),
                       guide = guide_colourbar(title.position = "top")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(400, 1000, 200)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(400, 1000, 200)) +
  coord_equal() +
  labs(x = "Wavelength (nm)", y = "Wavelength (nm)") +
  theme_paper() +
  theme(legend.position = "bottom", legend.key.width = unit(20, "pt"),
        legend.key.height = unit(4, "pt"), panel.grid = element_blank(),
        legend.title = element_text(size = rel(0.8)),
        legend.margin = margin(t = -2))

# ---- B. trait x wavelength correlation ---------------------------------------
tr_cols <- paste0(names(TRAIT_SHORT), "_blup")
tr_cols <- tr_cols[tr_cols %in% names(blups)]
bl <- blups |> mutate(key = tolower(str_squish(str_replace_all(genotype, "_", " "))))
common <- intersect(rownames(X), bl$key)
Bm <- as.matrix(bl[match(common, bl$key), tr_cols, drop = FALSE])
Xc <- X[common, , drop = FALSE]
corr_tw <- vapply(seq_len(ncol(Bm)), function(k)
  suppressWarnings(cor(Xc, Bm[, k], use = "pairwise.complete.obs"))[, 1], numeric(ncol(Xc)))
colnames(corr_tw) <- sub("_blup$", "", tr_cols)
dB <- as_tibble(corr_tw) |> mutate(wl = WL) |>
  pivot_longer(-wl, names_to = "trait", values_to = "r") |>
  mutate(trait = factor(TRAIT_TINY[trait], levels = unname(TRAIT_TINY)))

pB <- ggplot(dB, aes(wl, r, colour = trait)) +
  geom_hline(yintercept = 0, colour = INK[["muted"]], linewidth = 0.25) +
  geom_line(linewidth = 0.45) +
  scale_colour_manual(values = setNames(
    c(DEEP, SIG, "#E69F00", "#009E73", "#7A9CB5"), unname(TRAIT_TINY)), name = NULL) +
  scale_x_continuous(breaks = seq(400, 1000, 200)) +
  labs(x = "Wavelength (nm)", y = "Correlation with BLUP") +
  theme_paper() +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.key.height = unit(7, "pt"),
        legend.key.width = unit(10, "pt"), legend.margin = margin(t = -4),
        legend.text = element_text(size = rel(0.72)))

# ---- C. variance explained across wavelengths --------------------------------
# One-way ANOVA per wavelength, vectorised: for grouping factor g,
#   between-group SS / total SS, adjusted to a variance-component ratio.
scan_raw <- read_csv(file.path(DIR_RAW, "spectra_interpolated.csv"),
                     show_col_types = FALSE, progress = FALSE)
dict <- read_csv(file.path(DIR_RAW, "spectra_dictionary.csv"), show_col_types = FALSE) |>
  mutate(key = tolower(str_squish(str_replace_all(genotype, "_", " "))),
         Tags = as.character(sample_id))
scan_raw <- scan_raw |> mutate(Tags = as.character(Tags)) |>
  inner_join(dict |> select(Tags, key), by = "Tags") |>
  inner_join(panel |> select(key, species), by = "key")
Xs <- as.matrix(scan_raw[, wl_cols])

vc_ratio <- function(Xs, g) {
  g <- factor(g); n <- nrow(Xs); k <- nlevels(g)
  ni <- as.vector(table(g))
  grand <- colMeans(Xs)
  gm <- rowsum(Xs, g) / ni
  ssb <- colSums(ni * sweep(gm, 2, grand)^2)
  sst <- colSums(sweep(Xs, 2, grand)^2)
  ssw <- sst - ssb
  msb <- ssb / (k - 1); msw <- ssw / (n - k)
  n0 <- (n - sum(ni^2) / n) / (k - 1)
  s2a <- pmax((msb - msw) / n0, 0)
  s2a / (s2a + msw)
}
dC <- tibble(wl = WL,
             Accession = vc_ratio(Xs, scan_raw$key),
             Species   = vc_ratio(Xs, scan_raw$species)) |>
  pivot_longer(-wl, names_to = "level", values_to = "p")

pC <- ggplot(dC, aes(wl, p, colour = level)) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c(Accession = DEEP, Species = SIG), name = NULL) +
  scale_x_continuous(breaks = seq(400, 1100, 200)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Wavelength (nm)", y = "Variance explained") +
  theme_paper() +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.key.height = unit(7, "pt"),
        legend.key.width = unit(10, "pt"), legend.margin = margin(t = -4),
        legend.text = element_text(size = rel(0.72)))

# ---- D, E. PCA ---------------------------------------------------------------
# All 18 species, using the colour x shape scheme from 00_setup.R that Figure 1
# already uses. Eighteen categories is far past the point where hue alone can
# separate them, which is exactly why that scheme pairs a reused colour with a
# shape that encodes geographic group -- and reusing it means the reader learns
# one species legend for the whole paper rather than a new one per figure.
sp_present <- SPECIES_LEVELS[SPECIES_LEVELS %in% unique(pca$species)]
pca <- pca |> filter(species %in% sp_present) |>
  mutate(sp_lab = factor(species, levels = sp_present))

# Species names are set in italics. Built as plotmath expressions rather than
# with ggtext so the figure needs no extra package; the hybrid marker in
# "V. x doaniana" is rendered as a separate italic token because a multiplication
# sign inside a plotmath string does not survive every device.
italic_label <- function(sp) {
  lab <- SPECIES_LABELS[sp]
  vapply(lab, function(l) {
    if (grepl("\u00d7", l)) {
      # plotmath renders %*% as a multiplication sign; a literal U+00D7 inside a
      # quoted plotmath string comes out as "<U+00D7>" on some devices.
      parts <- trimws(strsplit(l, "\u00d7")[[1]])
      sprintf('italic("%s")%%*%%italic("%s")', parts[1], parts[2])
    } else sprintf('italic("%s")', l)
  }, character(1), USE.NAMES = FALSE)
}
SP_LAB_EXPR <- parse(text = italic_label(sp_present))

pca_panel <- function(src, ttl) {
  d <- pca |> filter(source == src)
  v1 <- unique(d$var1)[1] * 100; v2 <- unique(d$var2)[1] * 100
  ggplot(d, aes(PC1, PC2, colour = sp_lab, shape = sp_lab)) +
    geom_point(size = PT + 0.45, stroke = 0.3, alpha = 0.9) +
    scale_colour_manual(values = SPECIES_COLS[sp_present], name = NULL,
                        labels = SP_LAB_EXPR, drop = FALSE) +
    scale_shape_manual(values = SPECIES_SHAPES[sp_present], name = NULL,
                       labels = SP_LAB_EXPR, drop = FALSE) +
    labs(x = sprintf("PC1 (%.1f%%)", v1), y = sprintf("PC2 (%.1f%%)", v2), title = ttl) +
    theme_paper() +
    theme(plot.title = element_text(size = rel(0.95), hjust = 0.5),
          legend.position = "none")
}
pD <- pca_panel("markers", "Markers (10k)")
pE <- pca_panel("spectra", "Spectra")

# pD carries the legend; plot_layout(guides = "collect") lifts it to the row so
# it appears once. patchwork's `&` operator is avoided throughout -- it fails to
# dispatch against ggplot2 4.0's S7 themes.
pD <- pD + theme(legend.position = "bottom",
                 legend.key.height = unit(7, "pt"),
                 legend.key.width = unit(9, "pt"),
                 legend.text = element_text(size = rel(0.62)),
                 legend.margin = margin(t = -2)) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE,
                               override.aes = list(size = 1.5, stroke = 0.3)),
         shape  = guide_legend(nrow = 3, byrow = TRUE,
                               override.aes = list(size = 1.5, stroke = 0.3)))

# ---- F, G, H. relationship matrices ------------------------------------------
ord_keys <- panel |> filter(has_markers) |> arrange(species, genotype) |> pull(key)
Zg <- local({
  x <- read_table(file.path(DIR_RAW, "genotype_matrix_10k.raw"), progress = FALSE,
    col_types = cols(FID = col_character(), IID = col_character(), PAT = col_character(),
                     MAT = col_character(), SEX = col_character(),
                     PHENOTYPE = col_character(), .default = col_double()))
  Z <- as.matrix(x[, setdiff(names(x), c("FID","IID","PAT","MAT","SEX","PHENOTYPE"))])
  rownames(Z) <- tolower(str_squish(str_replace_all(x$IID, "_", " ")))
  Z
})
Kg <- kernel_G(Zg); dimnames(Kg) <- list(rownames(Zg), rownames(Zg))
Ks <- kernel_S(X);  dimnames(Ks) <- list(rownames(X), rownames(X))
sp_of <- setNames(panel$species, panel$key)
Kt <- kernel_T(sp_of[ord_keys]); dimnames(Kt) <- list(ord_keys, ord_keys)

mat_panel <- function(K, keys, ttl, thin = 2, cols = RAMP_CORR) {
  keys <- keys[keys %in% rownames(K)]
  keys <- keys[seq(1, length(keys), by = thin)]
  M <- K[keys, keys, drop = FALSE]
  lim <- quantile(abs(M[upper.tri(M)]), 0.98, na.rm = TRUE)
  d <- as_tibble(as.table(M), .name_repair = "minimal")
  names(d) <- c("i", "j", "v")
  d <- d |> mutate(i = match(i, keys), j = match(j, keys),
                   v = pmax(pmin(v, lim), -lim))
  ggplot(d, aes(i, j, fill = v)) +
    geom_raster() +
    scale_fill_gradientn(colours = cols, name = NULL) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = ttl) +
    theme_paper() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "none",
          plot.title = element_text(size = rel(0.9), hjust = 0.5))
}
pF <- mat_panel(Kg, ord_keys, "Genomic")
pG <- mat_panel(Ks, ord_keys[ord_keys %in% rownames(Ks)], "Spectral")
# The taxonomy kernel is binary, so it gets a two-colour scale of its own rather
# than the diverging ramp used for the two continuous kernels -- a diverging
# scale on 0/1 data implies a meaningful midpoint that does not exist.
pH <- mat_panel(Kt, ord_keys, "Taxonomy", cols = c(INK[["secondary"]], "#E8863C"))

# ---- I. accuracy by marker density and scenario ------------------------------
dI <- summ |>
  filter(model == "G", design %in% c("random", "holdout")) |>
  mutate(scen = factor(SCEN_LAB[scenario], levels = unname(SCEN_LAB)),
         trait = factor(TRAIT_SHORT[trait], levels = unname(TRAIT_SHORT)),
         density = factor(density, levels = c("5k", "10k", "20k", "50k")))

pI <- ggplot(dI, aes(density, r, colour = trait, group = trait)) +
  geom_hline(yintercept = 0, colour = INK[["muted"]], linewidth = 0.25) +
  geom_line(linewidth = 0.4) +
  geom_point(size = PT + 0.35, stroke = 0) +
  facet_wrap(~ scen, nrow = 1) +
  scale_colour_manual(values = setNames(
    c(DEEP, SIG, "#E69F00", "#009E73", "#7A9CB5"), unname(TRAIT_SHORT)), name = NULL) +
  labs(x = "Markers", y = "Predictive ability (r)") +
  theme_paper() +
  theme(legend.position = "bottom", legend.key.height = unit(8, "pt"),
        legend.text = element_text(size = rel(0.78)),
        axis.text.x = element_text(size = rel(0.75)),
        panel.spacing = unit(5, "pt"))

# ---- J. model comparison -----------------------------------------------------
dJ <- summ |>
  filter(design == "random", density == "10k", scenario == "within_2025") |>
  mutate(model = factor(MODEL_LAB[model], levels = unname(MODEL_LAB)),
         trait = factor(TRAIT_SHORT[trait], levels = unname(TRAIT_SHORT)))

pJ <- ggplot(dJ, aes(trait, r, fill = model)) +
  geom_hline(yintercept = 0, colour = INK[["muted"]], linewidth = 0.25) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  scale_fill_manual(values = MODEL_COL, name = NULL) +
  scale_x_discrete(labels = \(x) str_wrap(x, 12)) +
  labs(x = NULL, y = "Predictive ability (r)") +
  theme_paper() +
  theme(legend.position = "bottom", legend.key.height = unit(8, "pt"),
        legend.key.width = unit(9, "pt"),
        legend.text = element_text(size = rel(0.78)),
        axis.text.x = element_text(size = rel(0.8)))

# ---- assemble ----------------------------------------------------------------
tag_it <- function(p, t) p + labs(tag = t)
row1 <- tag_it(pA, "A") | tag_it(pB, "B") | tag_it(pC, "C")
row2 <- (tag_it(pD, "D") | tag_it(pE, "E") | tag_it(pF, "F") |
         tag_it(pG, "G") | tag_it(pH, "H")) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))
row3 <- tag_it(pI, "I")
row4 <- tag_it(pJ, "J")

fig <- wrap_elements(full = row1) / wrap_elements(full = row2) /
       wrap_elements(full = row3) / wrap_elements(full = row4) +
  # Panel I is squeezed by about 28% and the height goes to row 2, where the
  # PCA and matrix panels were cramped.
  plot_layout(heights = c(1.15, 1.20, 0.90, 1.05))

save_figure(fig, "Figure6", width = 7.28, height = 9.6)


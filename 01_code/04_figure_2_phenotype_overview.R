# ==============================================================================
# 04_Figure2.R
#
# FIGURE 3 -- Phenotypic variation, genetic structure and year-to-year
# consistency of adventitious-rooting traits. Merges the old Figures 3 and 4 and
# adds the evaluation-coverage and repeatability panels.
#
#   A  evaluation coverage: which years each accession was phenotyped in
#   B  raw trait distributions by year
#   C  callus-stage rooting BLUPs by species
#   D  post-transplant rooting BLUPs by species
#   E  correlations among trait BLUPs
#   F  repeatability (entry-mean basis)
#
# Season-to-season correlations are no longer plotted; they remain in
# 02_clean_data/year_stability.csv and drive the sensitivity fit shown in F.
#
# Input   02_clean_data/*  (from 03_phenotype_models.R)
# Output  03_figures/Figure2.{pdf,png,tiff}
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/04_figure_2_phenotype_overview.R"
  }))), "00_setup.R"))

need <- c("pheno_long.csv", "blups_by_trait.csv", "repeatability.csv",
          "repeatability_sensitivity.csv", "trait_correlations.csv",
          "year_stability.csv")
if (!all(file.exists(file.path(DIR_CLEAN, need)))) {
  source(file.path(DIR_CODE, "03_phenotype_models.R"))
}

TRAIT_LEVELS <- c("Callus-Stage", "Post-transplant", "Root weight",
                  "Cutting weight", "Transplant response")
TRAIT_SHORT  <- c("Callus-Stage" = "Callus-stage\nrooting",
                  "Post-transplant" = "Post-transplant\nrooting",
                  "Root weight" = "Root weight\n(g)",
                  "Cutting weight" = "Cutting weight\n(g)",
                  "Transplant response" = "Transplant response\nindex")
TRAIT_CODE   <- c("Callus-Stage" = "CS", "Post-transplant" = "PT",
                  "Root weight" = "RW", "Cutting weight" = "CW",
                  "Transplant response" = "TR")
TRAIT_STRIP  <- c("Callus-Stage" = "Callus-stage score",
                  "Post-transplant" = "Post-transplant score",
                  "Root weight" = "Root weight (g)",
                  "Cutting weight" = "Cutting weight (g)")

TRAIT_LABELS_ALL <- c("Callus-Stage"        = "Callus-stage rooting",
                      "Post-transplant"     = "Post-transplant rooting",
                      "Root weight"         = "Root weight",
                      "Cutting weight"      = "Cutting weight",
                      # named "index" on the axes of panels E and F so the
                      # figure matches the term used throughout the text
                      "Transplant response" = "Transplant response index")
TRAIT_LABELS_LONG <- TRAIT_LABELS_ALL

rd <- function(f) read_csv(file.path(DIR_CLEAN, f), show_col_types = FALSE)
pheno  <- rd("pheno_long.csv")   |> mutate(trait = factor(trait, TRAIT_LEVELS))
blupsw <- rd("blups_by_trait.csv")
repeat_all  <- rd("repeatability.csv")
repeat_sens <- rd("repeatability_sensitivity.csv")
trait_cor   <- rd("trait_correlations.csv")
year_stab   <- rd("year_stability.csv")

blups <- blupsw |>
  select(genotype, species, ends_with("_blup_z")) |>
  pivot_longer(ends_with("_blup_z"), names_to = "trait", values_to = "blup_z") |>
  mutate(trait = factor(sub("_blup_z$", "", trait), TRAIT_LEVELS),
         species = factor(species, SPECIES_LEVELS)) |>
  left_join(SPECIES_LUT |> select(species, group), by = "species")

# ==============================================================================
# A -- evaluation coverage
# ==============================================================================
# UpSet-style: bars are the number of accessions with each combination of years,
# the dot matrix underneath says which years those are.
cov_combo <- pheno |>
  distinct(genotype, year) |>
  group_by(genotype) |>
  summarise(combo = paste(sort(unique(year)), collapse = "+"),
            n_years = dplyr::n(), .groups = "drop") |>
  count(combo, n_years, name = "n_acc") |>
  arrange(desc(n_acc)) |>
  mutate(combo = factor(combo, levels = combo))

all_years <- sort(unique(as.character(pheno$year)))
cov_dots <- tidyr::expand_grid(combo = cov_combo$combo, year = all_years) |>
  mutate(present = purrr::map2_lgl(combo, year,
                                   \(c, y) grepl(y, as.character(c), fixed = TRUE))) |>
  left_join(cov_combo |> select(combo, n_acc), by = "combo")

pA_bar <- ggplot(cov_combo, aes(combo, n_acc)) +
  geom_col(fill = INK[["secondary"]], width = 0.68) +
  geom_text(aes(label = n_acc), vjust = -0.45, size = BASE_SIZE / .pt * 0.8,
            colour = INK[["primary"]]) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.18)),
                     name = "Accessions") +
  theme_paper() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        # the dot matrix below carries wider y labels; without this the shared
        # axis gutter pushes "Accessions" away from its own axis
        axis.title.y = element_text(margin = margin(r = 1)),
        plot.margin = margin(4, 4, 6, 4))

pA_dot <- ggplot(cov_dots, aes(combo, year)) +
  geom_point(aes(colour = present), size = 2.1) +
  geom_line(data = filter(cov_dots, present), aes(group = combo),
            colour = INK[["primary"]], linewidth = 0.5) +
  scale_colour_manual(values = c("TRUE" = INK[["primary"]], `FALSE` = INK[["faint"]]),
                      guide = "none") +
  scale_y_discrete(limits = rev(all_years)) +
  theme_paper() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(4, 4, 4, 4))

pA <- (pA_bar + labs(tag = "A") +
         theme(plot.tag = element_text(size = rel(1.7), face = "bold"))) /
      pA_dot + plot_layout(heights = c(1, 0.44))

# ==============================================================================
# B -- raw trait distributions by year
# ==============================================================================
pB <- pheno |>
  filter(trait != "Transplant response") |>
  mutate(trait = factor(TRAIT_STRIP[as.character(trait)],
                        levels = unname(TRAIT_STRIP[TRAIT_LEVELS[1:4]])),
         year = factor(year)) |>
  ggplot(aes(year, value, fill = year)) +
  geom_violin(colour = NA, alpha = 0.5, width = 0.9, scale = "width") +
  geom_boxplot(width = 0.20, outlier.shape = NA, linewidth = 0.28,
               colour = INK[["primary"]], fill = "white") +
  scale_fill_manual(values = YEAR_COLS, guide = "none") +
  facet_wrap(~ trait, scales = "free", nrow = 2) +
  scale_x_discrete(drop = TRUE) +
  labs(y = "Observed value", x = NULL) +
  theme_paper() +
  theme(panel.grid.major.x = element_blank(),
        strip.text = element_text(size = rel(0.85), face = "bold"),
        panel.spacing.x = unit(6, "pt"))

# ==============================================================================
# C, D -- BLUPs by species
# ==============================================================================
species_panel <- function(tr, show_y_title = TRUE) {
  d <- blups |> filter(trait == tr, !is.na(blup_z))
  ord <- d |> group_by(species) |>
    summarise(m = median(blup_z), n = dplyr::n(), .groups = "drop") |>
    arrange(m)
  d <- d |> mutate(species = factor(species, levels = ord$species))
  labs_n <- setNames(sprintf("%s (%d)", SPECIES_LABELS[as.character(ord$species)],
                             ord$n), ord$species)

  ggplot(d, aes(species, blup_z)) +
    geom_hline(yintercept = 0, colour = INK[["faint"]], linewidth = 0.3) +
    geom_boxplot(aes(colour = group), fill = NA, width = 0.62,
                 outlier.shape = NA, linewidth = 0.36) +
    geom_point(aes(colour = group), position = position_jitter(width = 0.16),
               size = 0.5, alpha = 0.5) +
    scale_colour_manual(values = GROUP_COLS, name = NULL, drop = FALSE) +
    scale_x_discrete(labels = labs_n) +
    coord_flip() +
    labs(x = NULL, y = if (show_y_title) "BLUP (SD units)" else NULL,
         title = TRAIT_LABELS_LONG[tr]) +
    theme_paper() +
    theme(axis.text.y = element_text(face = "italic", size = rel(0.85)),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = rel(1.0), face = "bold"))
}
pC <- species_panel("Callus-Stage")
pD <- species_panel("Post-transplant")

# ==============================================================================
# E -- correlations among trait BLUPs
# ==============================================================================
# Asterisk marks a correlation significant at alpha = 0.05 (two-sided).
cor_p <- trait_cor |>
  filter(method == "pearson") |>
  mutate(star = ifelse(!is.na(p) & p < 0.05, "*", "")) |>
  mutate(trait_x = factor(trait_x, TRAIT_LEVELS),
         trait_y = factor(trait_y, TRAIT_LEVELS)) |>
  filter(as.integer(trait_x) < as.integer(trait_y))

pE <- ggplot(cor_p, aes(trait_x, trait_y, fill = r)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = paste0(sub("^0\\.", ".", sprintf("%.2f", r)), star),
                colour = abs(r) > 0.55),
            size = BASE_SIZE / .pt * 1.05) +
  scale_fill_gradientn(colours = RAMP_DIV, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson r",
                       guide = guide_bar("Pearson r", width = 48)) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = INK[["primary"]]),
                      guide = "none") +
  scale_x_discrete(labels = \(x) TRAIT_LABELS_ALL[x], limits = TRAIT_LEVELS[-5]) +
  scale_y_discrete(labels = \(x) TRAIT_LABELS_ALL[x], limits = TRAIT_LEVELS[-1]) +
  labs(x = NULL, y = NULL, title = "Trait correlations") +
  theme_paper() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = rel(0.95), angle = 35, hjust = 1),
        axis.text.y = element_text(size = rel(0.95)),
        # only the lower triangle is drawn, so the bottom-right of the panel is
        # empty; the colourbar goes there instead of reserving a legend strip,
        # which would also leave a matching gap beside panel F
        legend.direction = "horizontal",
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.75)),
        legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size = rel(1.0), face = "bold", hjust = 0)) +
  legend_inside(0.99, 0.04, just = c(1, 0))

# ==============================================================================
# F -- repeatability
# ==============================================================================
# One estimate per trait, using every season that trait was measured in. Root
# weight and cutting weight were not recorded in 2026; that is stated in the
# figure legend rather than encoded as a second bar. The drop-2026 sensitivity
# fit remains in 02_clean_data/repeatability_sensitivity.csv.
rep_plot <- repeat_all |>
  transmute(trait, value = rep_generalised, lo = rep_gen_lo, hi = rep_gen_hi) |>
  mutate(trait = factor(trait, rev(TRAIT_LEVELS)))

# The bootstrap 95% intervals are in 02_clean_data/repeatability.csv and are
# quoted in the Results text; they are not drawn here.
pG <- ggplot(rep_plot, aes(value, trait)) +
  geom_col(fill = PAL[["blue"]], width = 0.60) +
  geom_text(aes(label = sprintf("%.2f", value), x = 0.02), hjust = 0,
            size = BASE_SIZE / .pt * 0.85, colour = "white", fontface = "bold") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(labels = \(x) TRAIT_LABELS_ALL[x]) +
  labs(x = "Repeatability (entry-mean basis)", y = NULL,
       title = "Repeatability") +
  theme_paper() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = rel(0.95)),
        axis.text.x = element_text(size = rel(0.95)),
        plot.title = element_text(size = rel(1.0), face = "bold", hjust = 0))

# ==============================================================================
# Assemble
# ==============================================================================
tagged <- function(p, tag) p + labs(tag = tag) +
  theme(plot.tag = element_text(size = rel(1.7), face = "bold"))

# patchwork aligns the left axis gutter across every row, so the very wide
# species labels in row 2 would push the y-axis titles of rows 1 and 3 far away
# from their own axes. Wrapping those rows makes them opaque to that alignment.
isolate_row <- function(p) patchwork::wrap_elements(full = p)

row1 <- (pA | tagged(pB, "B")) + plot_layout(widths = c(0.72, 1))
grp_legend <- ggplot(tibble(g = factor(GROUP_ORDER, GROUP_ORDER)),
                     aes(x = 1, y = 1, colour = g)) +
  geom_point(alpha = 0) +
  scale_colour_manual(values = GROUP_COLS, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(nrow = 1,
                               override.aes = list(alpha = 1, size = 2.2))) +
  theme_void(base_size = BASE_SIZE) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = rel(0.8)),
        legend.key.height = unit(8, "pt"))

row2 <- ((tagged(pC, "C") | tagged(pD, "D")) / grp_legend) +
  plot_layout(heights = c(1, 0.045))
# E carries tall rotated axis labels; isolating the two panels stops F's axis
# title being pushed down to align with them. The extra top margin keeps the
# E and F titles clear of the geographic-group legend sitting under C and D.
row3_margin <- theme(plot.margin = margin(t = 16, r = 4, b = 4, l = 4))
row3 <- (isolate_row(tagged(pE, "E") + row3_margin) |
         isolate_row(tagged(pG, "F") + row3_margin)) +
  plot_layout(widths = c(1, 1.08))

fig3 <- isolate_row(row1) / row2 / row3 +
  plot_layout(heights = c(0.80, 1.42, 0.96))

save_figure(fig3, "Figure2", width = 7.28, height = 8.70)



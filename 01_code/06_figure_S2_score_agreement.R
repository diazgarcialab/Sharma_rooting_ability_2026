# ==============================================================================
# 06_FigureS2.R
#
# TEST / SUPPLEMENTARY FIGURE -- season-to-season agreement matrices for the two
# ordinal rooting scores.
#
# Input   02_clean_data/pheno_long.csv   (from 03_phenotype_models.R)
# Output  03_figures/FigureS2.{pdf,png,tiff}
#         02_clean_data/score_agreement.csv
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/06_figure_S2_score_agreement.R"
  }))), "00_setup.R"))

suppressPackageStartupMessages(library(purrr))

N_BOOT     <- as.integer(Sys.getenv("ROOTING_NBOOT", unset = "1000"))
SCORE_MAX  <- 4L
SCORES     <- 0:SCORE_MAX

# Which trait x season values are already integers on the 0-4 scale.
NATIVE_INTEGER <- tibble::tribble(
  ~trait,            ~year,  ~native,
  "Callus-Stage",    "2024", TRUE,     # one score per bundle
  "Callus-Stage",    "2025", FALSE,    # mean of the individual cuttings
  "Callus-Stage",    "2026", TRUE,     # one score per bundle
  "Post-transplant", "2024", FALSE,    # mean of 6 plants
  "Post-transplant", "2025", FALSE,    # mean of 6 plants
  "Post-transplant", "2026", TRUE      # one score per bundle
)

TRAIT_LAB <- c("Callus-Stage" = "Callus-stage rooting",
               "Post-transplant" = "Post-transplant rooting")

# ---- 1. Score pairs ----------------------------------------------------------
pheno <- read_csv(file.path(DIR_CLEAN, "pheno_long.csv"), show_col_types = FALSE) |>
  filter(trait %in% names(TRAIT_LAB)) |>
  mutate(year = as.character(year))

wide <- pheno |>
  select(genotype, trait, year, value) |>
  pivot_wider(names_from = year, values_from = value)

pairs_tbl <- pheno |> distinct(trait, year) |>
  group_by(trait) |>
  summarise(y = list(sort(unique(year))), .groups = "drop") |>
  mutate(pp = map(y, \(v) as_tibble(t(combn(v, 2)), .name_repair = ~ c("y1", "y2")))) |>
  select(trait, pp) |> unnest(pp)

clip <- function(x) pmin(pmax(x, 0L), SCORE_MAX)
round_half_up <- function(x) clip(as.integer(floor(x + 0.5)))

native_of <- function(tr, yr) {
  v <- NATIVE_INTEGER$native[NATIVE_INTEGER$trait == tr & NATIVE_INTEGER$year == yr]
  stopifnot(length(v) == 1)
  v
}

pair_scores <- pmap(pairs_tbl, function(trait, y1, y2) {
  tr <- trait; ya <- y1; yb <- y2          # plain scalars, never column names
  d <- wide |> filter(.data$trait == tr) |>
    select(genotype, a = all_of(ya), b = all_of(yb)) |>
    filter(!is.na(a), !is.na(b))
  if (!nrow(d)) return(NULL)
  d$trait    <- tr
  d$y1       <- ya
  d$y2       <- yb
  d$a_native <- native_of(tr, ya)
  d$b_native <- native_of(tr, yb)
  d$a_int    <- clip(as.integer(round(d$a)))
  d$b_int    <- clip(as.integer(round(d$b)))
  d$a_up     <- round_half_up(d$a)
  d$b_up     <- round_half_up(d$b)
  as_tibble(d)
}) |> compact() |> bind_rows()

# ---- 2. Agreement statistics -------------------------------------------------

#' Quadratic-weighted Cohen's kappa on the 0..SCORE_MAX scale.
kappa_quad <- function(a, b, k = SCORE_MAX + 1L) {
  a <- factor(a, levels = SCORES); b <- factor(b, levels = SCORES)
  O <- table(a, b)
  n <- sum(O)
  if (n == 0) return(NA_real_)
  E <- outer(rowSums(O), colSums(O)) / n
  W <- outer(SCORES, SCORES, \(i, j) (i - j)^2)
  den <- sum(W * E)
  if (den == 0) return(NA_real_)
  1 - sum(W * O) / den
}

agreement <- pair_scores |>
  group_by(trait, y1, y2, a_native, b_native) |>
  group_modify(function(d, key) {
    kw   <- kappa_quad(d$a_int, d$b_int)
    kw_u <- kappa_quad(d$a_up,  d$b_up)
    boot <- if (N_BOOT > 0) replicate(N_BOOT, {
      i <- sample.int(nrow(d), replace = TRUE)
      kappa_quad(d$a_int[i], d$b_int[i])
    }) else NA_real_
    tibble(
      n            = nrow(d),
      exact        = mean(d$a_int == d$b_int),
      within_one   = mean(abs(d$a_int - d$b_int) <= 1),
      kappa_w      = kw,
      kappa_lo     = if (N_BOOT > 0) unname(quantile(boot, 0.025, na.rm = TRUE)) else NA_real_,
      kappa_hi     = if (N_BOOT > 0) unname(quantile(boot, 0.975, na.rm = TRUE)) else NA_real_,
      kappa_halfup = kw_u,
      spearman_raw = suppressWarnings(cor(d$a, d$b, method = "spearman")),
      pearson_raw  = cor(d$a, d$b)
    )
  }) |>
  ungroup() |>
  mutate(rounded = !(a_native & b_native),
         pair    = paste(y1, "vs", y2),
         # asterisk marks a pair where at least one season had to be rounded
         # from a mean; explained in the subtitle
         # Reviewer 1 asked where the Spearman correlation appears on the
         # figure. It is now printed on the strip beneath the kappa, on the
         # unrounded values, as the legend states.
         label   = sprintf("%s  \u00b7  %s%s\n\u03ba(w) %.2f  (%.2f\u2013%.2f)  \u00b7  \u03c1 %.2f\nn %d  \u00b7  exact %d%%  \u00b7  \u00b11 %d%%",
                           sub(" rooting$", "", TRAIT_LAB[trait]), pair,
                           ifelse(a_native & b_native, "", "*"),
                           kappa_w, kappa_lo, kappa_hi, spearman_raw, n,
                           round(100 * exact), round(100 * within_one)))

write_csv(agreement, file.path(DIR_CLEAN, "score_agreement.csv"))

# ---- 3. Matrices -------------------------------------------------------------
mats <- pair_scores |>
  count(trait, y1, y2, a_int, b_int, name = "count") |>
  group_by(trait, y1, y2) |>
  mutate(prop = count / sum(count)) |>
  ungroup() |>
  complete(nesting(trait, y1, y2), a_int = SCORES, b_int = SCORES,
           fill = list(count = 0L, prop = 0)) |>
  left_join(agreement |> select(trait, y1, y2, label, rounded),
            by = c("trait", "y1", "y2"))

lev <- agreement |> arrange(match(trait, names(TRAIT_LAB)), y1, y2) |> pull(label)
mats      <- mats      |> mutate(label = factor(label, levels = lev))
agreement <- agreement |> mutate(label = factor(label, levels = lev))

RAMP_COUNT <- c("#FFFFFF", "#E3EDF5", "#BBD5E8", "#83B4D6", "#4A8CBE", "#2A6296", "#153E66")

p <- ggplot(mats, aes(a_int, b_int)) +
  geom_tile(aes(fill = prop), colour = "white", linewidth = 0.7) +
  geom_text(aes(label = ifelse(count == 0, "", as.character(count)),
                colour = prop > 0.13),
            size = BASE_SIZE / .pt * 0.80) +
  geom_abline(slope = 1, intercept = 0, colour = alpha(INK[["primary"]], 0.45),
              linewidth = 0.35, linetype = "22") +
  scale_fill_gradientn(colours = RAMP_COUNT, limits = c(0, NA),
                       labels = scales::percent_format(accuracy = 1),
                       name = "Share of accessions",
                       guide = guide_bar("Share of accessions", width = 80)) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = INK[["primary"]]),
                      guide = "none") +
  scale_x_continuous(breaks = SCORES, expand = c(0, 0)) +
  scale_y_continuous(breaks = SCORES, expand = c(0, 0)) +
  facet_wrap(~ label, ncol = 3) +
  coord_fixed(clip = "off") +
  labs(x = "Score in the earlier season", y = "Score in the later season",
       subtitle = paste0("Test-retest agreement of the same accession between ",
                         "two seasons.\n",
                         "* at least one season was rounded to the nearest whole ",
                         "score from a mean; \u03ba(w) is computed on the rounded ",
                         "scores and \u03c1, the Spearman correlation, on the ",
                         "unrounded values.")) +
  theme_paper() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = INK[["faint"]],
                                    linewidth = 0.3),
        plot.subtitle = element_text(size = rel(0.8), colour = INK[["secondary"]],
                                     margin = margin(b = 8)),
        strip.text = element_text(size = rel(0.80), face = "bold", lineheight = 1.45,
                                  margin = margin(t = 2, b = 5)),
        panel.spacing.y = unit(14, "pt"),
        panel.spacing.x = unit(10, "pt"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.75)))

save_figure(p, "FigureS2", width = 7.28, height = 5.6)


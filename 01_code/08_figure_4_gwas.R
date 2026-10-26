# ==============================================================================
# 08_Figure4.R -- Figure 5: Manhattan, Q-Q and effect plots
#
#   Rscript 01_code/08_Figure4.R
#
# Reads whatever is in 02_clean_data/gwas_results/ and builds one column of
# panels per trait, so the figure grows as association tests finish. With one
# trait it is a single column; with five it is the full figure.
#
# Inputs (copied down from the server):
#   gwas_<slug>_blup_PCA10.assoc.txt   GEMMA output, one per trait
#   gwas_<slug>_blup_PCA10.log.txt     optional; PVE and n are read from it
#   peak_snps_annotated.tsv            lead SNPs, from extract_peak_genotypes.sh
#   peak_genotypes.raw                 PLINK --recode A calls at those SNPs
#
# The last two are only needed for the effect panels. Without them the script
# still produces the Manhattan and Q-Q rows and says what is missing.
#
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/08_figure_4_gwas.R"
  }))), "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(readr); library(stringr); library(tibble)
})

DIR_GWAS_RES <- file.path(DIR_CLEAN, "gwas_results")
if (!dir.exists(DIR_GWAS_RES))
  stop("no ", DIR_GWAS_RES, " -- copy the GEMMA output down first")

# ---- tunables ----------------------------------------------------------------
# Significance threshold. Bonferroni by default: alpha / (number of SNPs tested
# for that trait). "fdr" switches back to Benjamini-Hochberg.
#
#   THRESHOLD=fdr Rscript 01_code/08_Figure4.R
#
# Significance is recomputed HERE from the full p-value column, not read from
# peak_snps_annotated.tsv. That file records whatever threshold the extraction
# script used, and this way the figure cannot disagree with itself if the two
# were run under different settings.
THRESHOLD <- tolower(Sys.getenv("THRESHOLD", unset = "bonferroni"))
stopifnot(THRESHOLD %in% c("bonferroni", "fdr"))
ALPHA         <- 0.05
FDR_LEVEL     <- ALPHA
# Grid resolution is matched to the point size: bigger marks cover more area, so
# fewer cells are needed for the result to look identical to plotting everything.
GRID_NX       <- 1400      # plotting cells across the genome; lower these two
GRID_NY       <- 600       # to shrink the PDF further
SEED          <- 20260818
MAX_EFFECT_PANELS <- 6     # effect panels drawn per trait, when they are drawn

# Effect panels are OFF by default. With five traits the Manhattan and Q-Q rows
# already fill a page, and which loci deserve an effect plot is a judgement made
# after looking at the Manhattans, not before.
#
#   DRAW_EFFECTS=1 Rscript 01_code/08_Figure4.R
#       every lead SNP, into effect_plots
#
#   EFFECT_LOCI="chr02:19087797,chr08:14798595" Rscript 01_code/08_Figure4.R
#       only the named loci, matched on chromosome and position
DRAW_EFFECTS <- tolower(Sys.getenv("DRAW_EFFECTS", unset = "false")) %in%
                c("1", "true", "yes")
EFFECT_LOCI  <- trimws(strsplit(Sys.getenv("EFFECT_LOCI", unset = ""), ",")[[1]])
EFFECT_LOCI  <- EFFECT_LOCI[nzchar(EFFECT_LOCI)]
if (length(EFFECT_LOCI)) DRAW_EFFECTS <- TRUE

# MANHATTAN_ONLY=1 produces just the stacked Manhattan panels, one per trait, on
# a shared genome axis with tick labels on the bottom panel only. That is the
# layout for reading DOWN a column to see whether a peak in one trait sits at the
# same position in another; the Q-Q panels are not useful for that comparison and
# only push the panels apart.
MANHATTAN_ONLY <- tolower(Sys.getenv("MANHATTAN_ONLY", unset = "false")) %in%
                  c("1", "true", "yes")
# Width of the Q-Q column relative to the Manhattan beside it.
QQ_WIDTH <- as.numeric(Sys.getenv("QQ_WIDTH", unset = "0.20"))
# Width of the side Q-Q in Figure 4, as a fraction of its Manhattan.
QQ_SIDE_W <- as.numeric(Sys.getenv("QQ_SIDE_W", unset = "0.26"))
# Figure 4 is drawn at FIG_W inches and will be reduced to the 185 mm (7.28 in)
# two-column measure in print. Type is scaled by that ratio so it lands at the
# size theme_paper() intends AFTER the reduction, rather than 19% smaller.
FIG_W     <- as.numeric(Sys.getenv("FIG_W", unset = "9.00"))
TYPE_SCALE <- FIG_W / 7.28
BASE_SIZE <- BASE_SIZE * TYPE_SCALE
# With one or two traits everything fits on a page. With five it does not: five
# Manhattans and five Q-Q panels alone fill a full page, and up to thirty effect
# panels below them would be squeezed to a few millimetres each. So once the
# figure would exceed MAX_PAGE_IN the effect panels move to their own
# supplementary figure instead of being silently compressed.
MAX_PAGE_IN <- 9.6

# Peak labels. Significant SNPs arrive in clusters -- one QTL can contribute
# dozens -- so labels are placed on CLUMPED lead SNPs, not on every point:
# strongest first, then anything within LABEL_CLUMP_MB on the same chromosome is
# absorbed into it. MAX_LABELS_PER_CHR caps how many survive on one chromosome,
# which is what stops a dense region from collecting a stack of labels.
#
#   LABEL_LOCI="chr08:14798595,chr10:1210000"   label exactly these, nothing else
#   LABEL_CLUMP_MB=8 MAX_LABELS_PER_CHR=1        merge harder / label less
LABEL_CLUMP_MB     <- as.numeric(Sys.getenv("LABEL_CLUMP_MB", unset = "5"))
MAX_LABELS_PER_CHR <- as.integer(Sys.getenv("MAX_LABELS_PER_CHR", unset = "2"))
LABEL_LOCI <- trimws(strsplit(Sys.getenv("LABEL_LOCI", unset = ""), ",")[[1]])
LABEL_LOCI <- LABEL_LOCI[nzchar(LABEL_LOCI)]

TRAIT_LABELS <- c(
  callus_stage        = "Callus-stage rooting",
  post_transplant     = "Post-transplant rooting",
  root_weight         = "Root weight",
  cutting_weight      = "Cutting weight",
  transplant_response = "Transplant response index")

# Manhattan styling. The two chromosome colours are a POSITIONAL alternation --
# they carry no identity, only "this is a different chromosome from the last one"
# -- so they are deliberately recessive and close to neutral. The colour that has
# to carry meaning is the significant one, and that pair was validated rather
# than eyeballed: red against either band separates at deltaE 24.9 (deuteranopia)
# and 34.8 (tritanopia), well clear of the 8 threshold, and holds up in greyscale.
SIG_COL   <- "#C1272D"                      # significant SNPs and the threshold line
# Grey / orange alternation, validated against the red threshold line. At this
# lighter grey the worst adjacent separation is deltaE 12.5 under protanopia and
# 17.5 for normal vision -- still clear of the deltaE 8 floor. Going lighter
# again (#8A8A8A) drops normal vision to 15.1, right on the floor, so this is
# about as light as the grey can get and stay safe.
CHR_COLS  <- c("#7A7A7A", "#D97D2E")        # alternating chromosome bands
PT_SIZE   <- 0.85                           # null points
PT_SIZE_SIG <- PT_SIZE * 1.45               # significant points (ringed, not recoloured)

# ---- readers -----------------------------------------------------------------
read_assoc <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    out <- data.table::fread(path, showProgress = FALSE)
    tibble::as_tibble(out)
  } else {
    readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()),
                    progress = FALSE) |>
      mutate(across(c(ps, af, beta, se, p_wald), as.numeric))
  }
}

# Contig names look like VITVgdSC2_v1.0.hap1.chr08. Pull the number out; anything
# without a chrNN is an unplaced scaffold and is dropped from the Manhattan (it
# has no meaningful x position) but KEPT for the Q-Q and for lambda, where every
# test counts.
parse_chr <- function(x) {
  n <- suppressWarnings(as.integer(str_match(x, "(?i)chr0*([0-9]+)")[, 2]))
  n
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  stats::median(stats::qchisq(p, df = 1, lower.tail = FALSE)) / stats::qchisq(0.5, df = 1)
}

# ---- locate the association files --------------------------------------------
assoc_files <- list.files(DIR_GWAS_RES, pattern = "^gwas_.*\\.assoc\\.txt$",
                          full.names = TRUE)
if (!length(assoc_files))
  stop("no gwas_*.assoc.txt in ", DIR_GWAS_RES)

slug_of <- function(f) {
  s <- sub("^gwas_", "", basename(f)); s <- sub("\\.assoc\\.txt$", "", s)
  sub("_(blup|blue)_PCA[0-9]+$", "", s)
}
slugs <- vapply(assoc_files, slug_of, character(1), USE.NAMES = FALSE)
ord   <- order(match(slugs, names(TRAIT_LABELS)), na.last = TRUE)
assoc_files <- assoc_files[ord]; slugs <- slugs[ord]

message("association files: ", length(assoc_files))
for (i in seq_along(slugs)) message("  ", slugs[i], "  <- ", basename(assoc_files[i]))

# ---- read + annotate ---------------------------------------------------------
set.seed(SEED)
res <- list()
for (i in seq_along(assoc_files)) {
  a <- read_assoc(assoc_files[i]) |>
    filter(is.finite(p_wald), p_wald > 0, p_wald <= 1) |>
    mutate(chr_n = parse_chr(chr))
  a$adj <- if (THRESHOLD == "fdr") p.adjust(a$p_wald, method = "BH") else
                                   pmin(a$p_wald * nrow(a), 1)
  a$fdr <- a$adj                      # kept for downstream column compatibility

  lam <- lambda_gc(a$p_wald)
  n_sig <- sum(a$adj <= ALPHA)
  # Bonferroni gives a fixed line at alpha/m; BH gives the largest p that still
  # clears the level, so the line moves with each trait's p-distribution.
  p_thresh <- if (THRESHOLD == "bonferroni") ALPHA / nrow(a) else
              if (n_sig > 0) max(a$p_wald[a$adj <= ALPHA]) else NA_real_

  res[[slugs[i]]] <- list(a = a, lambda = lam, n_sig = n_sig,
                          p_thresh = p_thresh, n_snp = nrow(a),
                          n_unplaced = sum(is.na(a$chr_n)))
  message(sprintf("  %-20s %9d SNPs  lambda=%.3f  %d significant (%s, p <= %.2e)",
                  slugs[i], nrow(a), lam, n_sig, THRESHOLD, p_thresh))
}

# ---- GEMMA logs: analysed n and PVE ------------------------------------------
read_log <- function(slug) {
  f <- list.files(DIR_GWAS_RES, pattern = paste0("^gwas_", slug, ".*\\.log\\.txt$"),
                  full.names = TRUE)
  if (!length(f)) return(list(n = NA_integer_, pve = NA_real_, pve_se = NA_real_))
  tx <- readLines(f[1], warn = FALSE)
  grab <- function(pat) {
    ln <- grep(pat, tx, value = TRUE)
    if (!length(ln)) return(NA_real_)
    suppressWarnings(as.numeric(sub(".*=\\s*", "", ln[1])))
  }
  list(n = as.integer(grab("number of analyzed individuals")),
       pve = grab("pve estimate"), pve_se = grab("se\\(pve\\)"))
}
logs <- lapply(slugs, read_log); names(logs) <- slugs

# ---- Manhattan ---------------------------------------------------------------
# One shared x axis across traits: cumulative offsets computed from the union of
# chromosome extents, so panels stack in register.
all_pos <- bind_rows(lapply(res, \(r) r$a |> filter(!is.na(chr_n)) |>
                              select(chr_n, ps)))
chr_span <- all_pos |> group_by(chr_n) |> summarise(len = max(ps), .groups = "drop") |>
  arrange(chr_n) |>
  mutate(offset = cumsum(as.numeric(len)) - as.numeric(len),
         centre = offset + len / 2)
GENOME_LEN <- sum(as.numeric(chr_span$len))

# One SNP per occupied plotting cell, plus every FDR-significant SNP regardless.
# See the note at the top of the file for why this is not a p-value cutoff.
thin_for_plot <- function(d) {
  if (!nrow(d)) return(d)
  xr <- range(d$x); yr <- range(d$logp)
  dx <- diff(xr); dy <- diff(yr)
  if (dx <= 0 || dy <= 0) return(d)
  cell <- paste(floor((d$x - xr[1]) / dx * GRID_NX),
                floor((d$logp - yr[1]) / dy * GRID_NY), sep = "_")
  keep <- !duplicated(cell) | d$sig
  d[keep, , drop = FALSE]
}

manhattan_panel <- function(slug, show_x = TRUE, note_thresh = FALSE) {
  r <- res[[slug]]
  d <- r$a |> filter(!is.na(chr_n)) |>
    left_join(chr_span |> select(chr_n, offset), by = "chr_n") |>
    mutate(x = offset + ps, logp = -log10(p_wald),
           band = factor(chr_n %% 2),
           sig  = adj <= ALPHA) |>
    thin_for_plot()

  ymax <- max(d$logp, na.rm = TRUE)
  p <- ggplot() +
    geom_point(data = filter(d, !sig), aes(x, logp, colour = band),
               size = PT_SIZE, stroke = 0, alpha = 0.85) +
    scale_colour_manual(values = setNames(CHR_COLS, c("0", "1")), guide = "none")
  # Significant SNPs keep their chromosome's colour and are distinguished by the
  # MARK instead: a slightly larger filled circle with a dark ring. Hue already
  # carries chromosome identity here, so overloading it with significance made
  # the red blocks read as a separate data series.
  if (any(d$sig))
    p <- p + geom_point(data = filter(d, sig), aes(x, logp, fill = band),
                        shape = 21, size = PT_SIZE_SIG, stroke = 0.35,
                        colour = INK[["primary"]]) +
      scale_fill_manual(values = setNames(CHR_COLS, c("0", "1")), guide = "none")
  if (is.finite(r$p_thresh))
    p <- p + geom_hline(yintercept = -log10(r$p_thresh), colour = SIG_COL,
                        linewidth = 0.4, linetype = "22")

  # ---- peak labels ----------------------------------------------------------
  lab <- d[d$sig, , drop = FALSE]
  if (nrow(lab)) {
    lab <- lab[order(lab$p_wald), ]
    if (length(LABEL_LOCI)) {
      want <- do.call(rbind, lapply(strsplit(LABEL_LOCI, ":"), \(z)
        data.frame(c = sub("^chr", "", z[1]), p = as.numeric(z[2]))))
      keep <- vapply(seq_len(nrow(lab)), \(j)
        any(sprintf("%02d", lab$chr_n[j]) == sprintf("%02d", as.integer(want$c)) &
            abs(lab$ps[j] - want$p) < 1e5), logical(1))
      lab <- lab[keep, , drop = FALSE]
    } else {
      keep <- logical(nrow(lab))
      # per-chromosome counter; indexing a named vector with an absent name
      # returns NA rather than NULL, so `%||%` does not catch it
      ncr <- setNames(integer(length(unique(lab$chr_n))),
                      as.character(unique(lab$chr_n)))
      for (j in seq_len(nrow(lab))) {
        cj <- as.character(lab$chr_n[j])
        near <- any(keep & lab$chr_n == lab$chr_n[j] &
                    abs(lab$ps - lab$ps[j]) < LABEL_CLUMP_MB * 1e6)
        if (near) next
        if (ncr[[cj]] >= MAX_LABELS_PER_CHR) next
        keep[j] <- TRUE; ncr[[cj]] <- ncr[[cj]] + 1L
      }
      lab <- lab[keep, , drop = FALSE]
    }
  }
  if (nrow(lab)) {
    lab$txt <- sprintf("chr%02d: %.1f Mb", lab$chr_n, lab$ps / 1e6)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      # ggrepel resolves collisions properly and draws a leader line back to the
      # peak, which fixed placement cannot do once two peaks are close together.
      p <- p + ggrepel::geom_text_repel(
        data = lab, aes(x, logp, label = txt),
        size = 1.95 * TYPE_SCALE, colour = INK[["primary"]], seed = 11,
        segment.colour = INK[["muted"]], segment.size = 0.22,
        min.segment.length = 0, box.padding = 0.28, point.padding = 0.15,
        max.overlaps = Inf, direction = "both", force = 3,
        nudge_y = 0.35, xlim = c(0, GENOME_LEN))
    } else {
      # Fallback: stagger alternate labels vertically so at least neighbouring
      # ones do not sit on each other, and anchor edge labels inward.
      message("note: install.packages(\"ggrepel\") for collision-free peak labels")
      lab <- lab[order(lab$x), ]
      lab$dy <- c(0, 0.55)[seq_len(nrow(lab)) %% 2 + 1]
      lab$hj <- ifelse(lab$x < GENOME_LEN * 0.08, 0,
                ifelse(lab$x > GENOME_LEN * 0.92, 1, 0.5))
      p <- p + geom_text(data = lab, aes(x, logp + dy, label = txt, hjust = hj),
                         vjust = -0.75, size = 1.95 * TYPE_SCALE, colour = INK[["primary"]])
    }
  }

  # A red note on the threshold line, first panel only, so the reader learns what
  # the line is once instead of five times. Placed at the left, where the right
  # end of the line collided with a peak label.
  if (note_thresh && is.finite(r$p_thresh)) {
    # Mid-genome, over chromosome 10 where the line is usually clear of peaks.
    x_note <- chr_span$centre[chr_span$chr_n == 10]
    if (!length(x_note)) x_note <- GENOME_LEN / 2
    p <- p + annotate("text", x = x_note[1], y = -log10(r$p_thresh),
                      label = "Bonferroni, 0.05", hjust = 0.5, vjust = -0.55,
                      size = 1.95 * TYPE_SCALE, colour = SIG_COL)
  }

  p +
    scale_x_continuous(breaks = chr_span$centre, labels = chr_span$chr_n,
                       expand = expansion(mult = 0.005),
                       limits = c(0, GENOME_LEN)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20)),
                       limits = c(0, NA)) +
    labs(x = if (show_x) "Chromosome" else NULL,
         y = expression(-log[10] * italic(p)),
         title = TRAIT_LABELS[slug] %||% slug) +
    theme_paper() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = if (show_x) element_text(size = rel(0.7)) else element_blank(),
          axis.title.x = if (show_x) element_text() else element_blank())
}

# ---- Q-Q ---------------------------------------------------------------------
qq_panel <- function(slug) {
  r <- res[[slug]]
  p <- sort(r$a$p_wald)
  n <- length(p)
  exp_p <- (seq_len(n) - 0.5) / n
  d <- tibble(obs = -log10(p), exp = -log10(exp_p),
              lo = -log10(qbeta(0.975, seq_len(n), n - seq_len(n) + 1)),
              hi = -log10(qbeta(0.025, seq_len(n), n - seq_len(n) + 1)))
  # thin the bulk for drawing only; the tail is what the eye reads
  idx <- which(d$obs > 2)
  bulk <- setdiff(seq_len(n), idx)
  if (length(bulk) > 20000) bulk <- sort(sample(bulk, 20000))
  dd <- d[sort(c(idx, bulk)), ]

  ggplot(dd, aes(exp, obs)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = INK[["faint"]], alpha = 0.55) +
    geom_abline(slope = 1, intercept = 0, colour = INK[["muted"]],
                linewidth = 0.3, linetype = "22") +
    geom_point(size = 0.7, stroke = 0, colour = CHR_COLS[1]) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.6, size = 2.3,
             colour = INK[["primary"]],
             label = sprintf("lambda == %.3f", r$lambda), parse = TRUE) +
    labs(x = expression("Expected " * -log[10] * italic(p)),
         y = expression("Observed " * -log[10] * italic(p))) +
    theme_paper()
}

# A narrow Q-Q for the side column: same data as qq_panel, stripped down so it
# reads at a fifth of the width. Axis titles only on the bottom row.
qq_panel_narrow <- function(slug, show_x = FALSE) {
  r <- res[[slug]]
  pv <- sort(r$a$p_wald); n <- length(pv)
  exp_p <- (seq_len(n) - 0.5) / n
  d <- tibble(obs = -log10(pv), exp = -log10(exp_p),
              lo = -log10(qbeta(0.975, seq_len(n), n - seq_len(n) + 1)),
              hi = -log10(qbeta(0.025, seq_len(n), n - seq_len(n) + 1)))
  idx <- which(d$obs > 2); bulk <- setdiff(seq_len(n), idx)
  if (length(bulk) > 6000) bulk <- sort(sample(bulk, 6000))
  dd <- d[sort(c(idx, bulk)), ]
  ggplot(dd, aes(exp, obs)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = INK[["faint"]], alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, colour = INK[["muted"]],
                linewidth = 0.3, linetype = "22") +
    # Same point size as the Manhattan beside it, so the two panels read as one
    # figure. Lambda is printed on the panel: Reviewer 1 asked for it to be
    # visible rather than only in gwas_summary.csv.
    geom_point(size = PT_SIZE, stroke = 0, colour = CHR_COLS[1]) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.12, vjust = 1.5, size = 2.3 * TYPE_SCALE,
             colour = INK[["secondary"]],
             label = sprintf("lambda == %.2f", r$lambda), parse = TRUE) +
    labs(x = if (show_x) "Expected" else NULL, y = "Observed") +
    theme_paper() +
    theme(axis.text = element_text(size = rel(0.78)),
          axis.title = element_text(size = rel(0.82)),
          plot.margin = margin(2, 4, 2, 2))
}

# ---- effect panels -----------------------------------------------------------
# Variants are matched on chromosome and position, never on variant ID. This
# VCF has no IDs -- every .bim row is named "." -- so a name-based join silently
# matches everything or nothing. peak_variants.bim lists the extracted variants
# in the same order as the .raw genotype columns, which is the reliable mapping.
peaks_file <- file.path(DIR_GWAS_RES, "peak_snps_annotated.tsv")
raw_file   <- file.path(DIR_GWAS_RES, "peak_genotypes.raw")
bim_file   <- file.path(DIR_GWAS_RES, "peak_variants.bim")
have_eff   <- DRAW_EFFECTS && file.exists(peaks_file) && file.exists(raw_file)

# A .raw of a few hundred samples x a few variants is tens of KB. Anything near
# a gigabyte means the extraction matched far more than the lead SNPs, and
# reading it will take the session down instead of erroring cleanly.
if (have_eff) {
  raw_mb <- file.size(raw_file) / 1024^2
  if (raw_mb > 200) {
    message(sprintf(
      "\nSKIPPING the effect panels: %s is %.0f MB.\n%s\n%s\n%s",
      basename(raw_file), raw_mb,
      "  A file this size means --extract matched far more than the lead SNPs,",
      "  which happens when the VCF has no variant IDs and every .bim row is '.'.",
      "  Re-run the current extract_peak_genotypes.sh, which selects by position."))
    have_eff <- FALSE
  }
}

effect_panels <- list()
if (!DRAW_EFFECTS) {
  message("\nEffect panels are off. Look at the Manhattans first, then re-run with")
  message("  EFFECT_LOCI=\"chr02:19087797,chr08:14798595\" Rscript 01_code/08_Figure4.R")
  message("  (or DRAW_EFFECTS=1 for every lead SNP).")
} else if (!have_eff) {
  message("\nNOTE: no peak_snps_annotated.tsv / peak_genotypes.raw in ", DIR_GWAS_RES)
  message("      Run extract_peak_genotypes.sh on the server and copy both down.")
  message("      Building the figure without effect panels for now.")
} else {
  peaks <- read_tsv(peaks_file, col_types = cols(.default = col_character())) |>
    mutate(across(c(ps, af, beta, se, p_wald, fdr), as.numeric))

  raw <- read_table(raw_file, col_types = cols(.default = col_character()),
                    progress = FALSE)
  gcols <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))

  # Column -> (chr, pos) by ROW ORDER of the companion .bim, which is how PLINK
  # writes both files. Falls back to parsing chr:pos out of the column name when
  # the .bim was not copied down.
  if (file.exists(bim_file)) {
    bim <- read_tsv(bim_file, col_names = c("chr", "id", "cm", "ps", "a1", "a2"),
                    col_types = cols(.default = col_character())) |>
      mutate(ps = as.numeric(ps))
    if (nrow(bim) != length(gcols))
      stop("peak_variants.bim has ", nrow(bim), " variants but peak_genotypes.raw ",
           "has ", length(gcols), " genotype columns -- they are out of step; ",
           "re-run extract_peak_genotypes.sh and copy both down together.")
    colmap <- tibble(snp_allele = gcols, chr = bim$chr, ps = bim$ps)
  } else {
    message("note: peak_variants.bim not found; falling back to parsing the ",
            ".raw column names for chr:pos")
    stem <- sub("_[ACGT]+$", "", gcols)
    colmap <- tibble(snp_allele = gcols,
                     chr = sub(":[0-9]+$", "", stem),
                     ps  = suppressWarnings(as.numeric(sub("^.*:", "", stem))))
    if (any(is.na(colmap$ps)))
      stop("cannot recover positions from the .raw column names and there is no ",
           "peak_variants.bim -- re-run extract_peak_genotypes.sh")
  }

  geno <- raw |> select(IID, all_of(gcols)) |>
    mutate(across(all_of(gcols), as.numeric)) |>
    pivot_longer(all_of(gcols), names_to = "snp_allele", values_to = "dosage") |>
    left_join(colmap, by = "snp_allele") |>
    mutate(counted = sub("^.*_", "", snp_allele))

  # library ID -> accession -> BLUP
  fam_all <- read_tsv(file.path(DIR_CLEAN, "gwas", "fam_phenotypes_all.txt"),
                      col_names = FALSE, col_types = cols(.default = col_character()))
  lib2acc <- tibble(IID = fam_all$X1, genotype = fam_all$X3)
  blups <- read_csv(file.path(DIR_CLEAN, "blups_by_trait.csv"),
                    show_col_types = FALSE)

  slug2trait <- setNames(names(TRAIT_LABELS), names(TRAIT_LABELS))
  TRAIT_COLNAME <- c(
    callus_stage = "Callus-Stage", post_transplant = "Post-transplant",
    root_weight = "Root weight", cutting_weight = "Cutting weight",
    transplant_response = "Transplant response")

  for (slug in slugs) {
    pk <- peaks |> filter(trait == slug) |> arrange(p_wald)
    if (length(EFFECT_LOCI)) {
      want <- tibble(tag = EFFECT_LOCI) |>
        separate(tag, into = c("chr_tag", "ps_tag"), sep = ":", convert = TRUE)
      pk <- pk |>
        mutate(chr_n = parse_chr(chr)) |>
        filter(paste0("chr", sprintf("%02d", chr_n)) %in% want$chr_tag &
               ps %in% want$ps_tag)
    } else {
      pk <- pk |> slice_head(n = MAX_EFFECT_PANELS)
    }
    if (!nrow(pk)) next
    tcol <- paste0(TRAIT_COLNAME[[slug]], "_blup")
    if (!tcol %in% names(blups)) { message("no BLUP column ", tcol); next }

    ph <- blups |> select(genotype, value = all_of(tcol))

    for (k in seq_len(nrow(pk))) {
      d <- geno |> filter(chr == pk$chr[k], ps == pk$ps[k]) |>
        left_join(lib2acc, by = "IID") |>
        left_join(ph, by = "genotype") |>
        filter(!is.na(dosage), !is.na(value))
      if (!nrow(d)) {
        message("  no genotypes for ", pk$chr[k], ":", pk$ps[k], " -- skipping")
        next
      }

      cnt <- pk$allele1[k]; oth <- pk$allele0[k]
      # PLINK counts the A1 allele in --recode A; the column name records which
      lab <- function(g) c(paste0(oth, oth), paste0(oth, cnt), paste0(cnt, cnt))[g + 1]
      d <- d |> mutate(gt = factor(lab(dosage),
                                   levels = c(paste0(oth, oth), paste0(oth, cnt),
                                              paste0(cnt, cnt))))
      nlab <- d |> count(gt) |> mutate(txt = paste0("n=", n))
      # raw per-copy slope, for comparison with GEMMA's structure-corrected beta
      raw_beta <- tryCatch(unname(coef(lm(value ~ dosage, data = d))[2]),
                           error = function(e) NA_real_)
      chrn <- parse_chr(pk$chr[k])
      ttl <- sprintf("chr%s : %.2f Mb", ifelse(is.na(chrn), "Un", sprintf("%02d", chrn)),
                     pk$ps[k] / 1e6)
      # Plain text, not plotmath: an NA slope or an unusual number format makes
      # parse() fail at draw time, which is a silly way to lose a figure.
      # Two short lines rather than one long one: at 7.28 in across three
      # columns a single line of this overflows the panel and gets clipped.
      sub <- sprintf("beta %.3f (GEMMA), %.3f (raw)\np = %s",
                     pk$beta[k], raw_beta,
                     formatC(pk$p_wald[k], format = "e", digits = 1))

      effect_panels[[sprintf("%s|%s:%.0f", slug, pk$chr[k], pk$ps[k])]] <-
        ggplot(d, aes(gt, value)) +
        geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.28,
                     colour = INK[["primary"]], fill = "white") +
        geom_jitter(width = 0.16, height = 0, size = 0.85, stroke = 0,
                    alpha = 0.45, colour = CHR_COLS[1]) +
        stat_summary(fun = mean, geom = "point", size = 1.5, shape = 18,
                     colour = SIG_COL) +
        geom_text(data = nlab, aes(gt, y = -Inf, label = txt), vjust = -0.8,
                  size = 1.9, colour = INK[["secondary"]], inherit.aes = FALSE) +
        labs(x = NULL, y = paste0(TRAIT_LABELS[slug], " BLUP"),
             title = ttl, subtitle = sub) +
        theme_paper() +
        theme(plot.title = element_text(size = rel(0.95), margin = margin(b = 3)),
              plot.subtitle = element_text(size = rel(0.72), lineheight = 1.15,
                                           colour = INK[["secondary"]],
                                           margin = margin(b = 5)),
              axis.title.y = element_text(size = rel(0.8)))
    }
  }
  message("\neffect panels built: ", length(effect_panels))
}

# ---- assemble ----------------------------------------------------------------
# In the stacked layout only the bottom panel carries the chromosome axis.
man_row <- lapply(seq_along(slugs), \(i)
  manhattan_panel(slugs[i], show_x = !MANHATTAN_ONLY || i == length(slugs),
                  note_thresh = i == 1L))
qq_row  <- lapply(slugs, qq_panel)

if (MANHATTAN_ONLY) {
  # Faint vertical rules at every lead position, drawn on ALL panels. The point
  # of the stack is to see whether a peak in one trait sits where another trait
  # also has one, and reading that off unaided across 500 Mb of x axis is hard.
  # GUIDES=0 turns them off.
  show_guides <- !tolower(Sys.getenv("GUIDES", unset = "1")) %in% c("0", "false", "no")
  guide_x <- numeric(0)
  if (show_guides) {
    # Lead positions taken from THIS script's significance calls, clumped at
    # 1 Mb, so the guides always match the points that are marked significant.
    gp <- bind_rows(lapply(slugs, function(sl) {
      aa <- res[[sl]]$a
      aa <- aa[aa$adj <= ALPHA & !is.na(aa$chr_n), , drop = FALSE]
      if (!nrow(aa)) return(NULL)
      aa <- aa[order(aa$p_wald), ]
      keep <- logical(nrow(aa))
      for (j in seq_len(nrow(aa))) {
        if (!any(keep & aa$chr_n == aa$chr_n[j] &
                 abs(aa$ps - aa$ps[j]) < 1e6)) keep[j] <- TRUE
      }
      aa[keep, c("chr_n", "ps")]
    }))
    if (nrow(gp)) {
      gp <- gp |> left_join(chr_span |> select(chr_n, offset), by = "chr_n") |>
        mutate(x = offset + ps)
      guide_x <- sort(unique(round(gp$x)))
    }
  }
  add_guides <- function(p) {
    if (!length(guide_x)) return(p)
    p$layers <- c(list(geom_vline(xintercept = guide_x, colour = INK[["faint"]],
                                  linewidth = 0.35)), p$layers)
    p
  }
  man_row <- lapply(man_row, add_guides)
  mp <- Map(function(p, t) p + labs(tag = t), man_row, LETTERS[seq_along(man_row)])
  # tighter vertical spacing than the default so the eye can follow a locus down
  # the stack; the shared x scale is what makes the comparison valid.
  mp <- lapply(mp, \(p) p + theme(plot.margin = margin(2, 4, 2, 4)))

  # Each trait gets its own Q-Q in a narrow right-hand column. wrap_elements per
  # row keeps `/` from flattening the rows into one grid.
  qn <- lapply(seq_along(slugs), \(i)
    qq_panel_narrow(slugs[i], show_x = i == length(slugs)))
  rows <- lapply(seq_along(mp), \(i)
    wrap_elements(full = mp[[i]] + qn[[i]] + plot_layout(widths = c(1, QQ_WIDTH))))
  fig_m <- wrap_plots(rows, ncol = 1)
  H_m <- 1.45 * length(mp) + 0.55
  save_figure(fig_m, "Figure4_manhattans_only", width = 7.28, height = H_m)
  cat("\n----- stacked Manhattan panels -----\n")
  for (i in seq_along(slugs))
    cat(sprintf("  %s  %-22s lambda %.3f   %d SNPs (%s, p <= %.2e)\n",
                LETTERS[i], TRAIT_LABELS[slugs[i]] %||% slugs[i],
                res[[slugs[i]]]$lambda, res[[slugs[i]]]$n_sig, THRESHOLD,
                res[[slugs[i]]]$p_thresh))
  message("\ndone")
  quit(save = "no", status = 0)
}

panels <- c(man_row, qq_row, unname(effect_panels))
tags   <- LETTERS[seq_along(panels)]
panels <- Map(function(p, t) p + labs(tag = t), panels, tags)

n_tr  <- length(slugs)
n_eff <- length(effect_panels)

# Manhattans are wide, Q-Q and effect panels are square-ish. Manhattans stack one
# per row; the Q-Q plots share a row; the effect plots fill a grid below.
# Each Manhattan sits beside its own Q-Q rather than in a separate block, so
# that the trait title above the pair clearly governs both panels. The Q-Q keeps
# its own tag (F-J). wrap_elements per row stops `/` flattening the rows into
# one grid.
qq_side <- lapply(seq_along(slugs), \(i)
  qq_panel_narrow(slugs[i], show_x = i == n_tr) +
    labs(tag = LETTERS[n_tr + i]))
man_rows <- lapply(seq_len(n_tr), \(i)
  wrap_elements(full = panels[[i]] + qq_side[[i]] +
                  plot_layout(widths = c(1, QQ_SIDE_W))))
man_block <- wrap_plots(man_rows, ncol = 1)

# Does everything fit on one page?
eff_ncol_try <- if (n_eff <= 3) max(n_eff, 1L) else if (n_eff == 4) 2L else 3L
eff_rows_try <- if (n_eff) ceiling(n_eff / eff_ncol_try) else 0L
H_all <- 2.0 * n_tr + 2.3 + 2.15 * eff_rows_try
SPLIT <- n_eff > 0 && H_all > MAX_PAGE_IN

if (SPLIT) {
  message(sprintf(
    "\n%d traits and %d effect panels would need %.1f in of height.\n%s\n%s",
    n_tr, n_eff, H_all,
    "  -> Figure 5 keeps the Manhattan and Q-Q panels;",
    "     the effect panels go to effect_plots."))
}

if (n_eff && !SPLIT) {
  # Balance the grid rather than always using 3 columns: 4 panels read better as
  # 2 x 2 than as 3 + 1 stranded on its own row.
  eff_ncol  <- eff_ncol_try
  eff_rows  <- eff_rows_try
  eff_items <- panels[2 * n_tr + seq_len(n_eff)]
  # y-axis title on the left column only; repeating it in every cell is noise
  # The y-axis title is deliberately kept on EVERY effect panel, not just the
  # left column. Two reasons, one cosmetic and one substantive:
  #
  #   ggplot aligns plot.title to the panel, and these panels have single-digit
  #   y labels, so without the axis title the panel starts far enough left that
  #   the title collides with the tag ("Dchr08 : 14.80 Mb"). Widening
  #   plot.margin does not help: the tag is positioned in npc of the whole plot,
  #   so the margin moves tag and title together.
  #
  #   More importantly, once several traits are in the figure their effect
  #   panels share one grid, and a panel labelled only "chr08 : 14.80 Mb" gives
  #   no clue which trait it belongs to.
  eff_block <- wrap_plots(eff_items, ncol = eff_ncol)
  # wrap_elements around each block: `/` flattens nested patchworks, so with
  # five Manhattans the outer layout ends up with 6+ children while `heights`
  # supplies 3, and patchwork aborts with "Need 6 panels, but together nrow and
  # ncol only provide 2".
  fig <- wrap_elements(full = man_block) / wrap_elements(full = eff_block) +
    plot_layout(heights = c(1.15 * n_tr, 1.3 * eff_rows))
  H <- 2.0 * n_tr + 1.2 + 2.15 * eff_rows
} else {
  fig <- man_block
  H <- 2.0 * n_tr + 1.2
}
# Wider than the usual 185 mm two-column measure: the Q-Q now shares the row
# with the Manhattan, so the pair needs the extra width to keep both legible.
W <- FIG_W

# No silent capping. If the figure genuinely needs more than a page, say so
# rather than shrinking every panel until nothing is legible.
if (H > MAX_PAGE_IN) {
  warning(sprintf("figure is %.1f in tall, more than a printed page", H))
}
save_figure(fig, "Figure4", width = W, height = H)

# ---- supplementary: effect panels --------------------------------------------
if (SPLIT) {
  eff_ncol  <- 3L
  eff_rows  <- ceiling(n_eff / eff_ncol)
  eff_items <- panels[2 * n_tr + seq_len(n_eff)]
  eff_block <- wrap_plots(eff_items, ncol = eff_ncol)
  H_s <- min(2.15 * eff_rows + 0.3, 9.4)
  # More panels than fit on one page: split into numbered parts rather than
  # shrinking them.
  per_page <- floor((9.4 - 0.3) / 2.15) * eff_ncol
  if (n_eff > per_page) {
    parts <- split(seq_len(n_eff), ceiling(seq_len(n_eff) / per_page))
    for (i in seq_along(parts)) {
      blk <- wrap_plots(eff_items[parts[[i]]], ncol = eff_ncol)
      rws <- ceiling(length(parts[[i]]) / eff_ncol)
      save_figure(blk, sprintf("effect_plots_part%d", i),
                  width = W, height = 2.15 * rws + 0.3)
    }
    message("effect panels written as ", length(parts), " parts")
  } else {
    save_figure(eff_block, "effect_plots", width = W, height = H_s)
  }
}

# ---- machine-readable summary + caption --------------------------------------
summ <- tibble(
  trait      = slugs,
  label      = unname(TRAIT_LABELS[slugs]),
  n_snp      = vapply(res, \(r) r$n_snp, numeric(1)),
  n_unplaced = vapply(res, \(r) r$n_unplaced, numeric(1)),
  lambda     = vapply(res, \(r) r$lambda, numeric(1)),
  n_sig_fdr  = vapply(res, \(r) r$n_sig, numeric(1)),
  p_threshold= vapply(res, \(r) r$p_thresh, numeric(1)),
  n_analysed = vapply(slugs, \(s) as.numeric(logs[[s]]$n), numeric(1)),
  pve        = vapply(slugs, \(s) logs[[s]]$pve, numeric(1)),
  pve_se     = vapply(slugs, \(s) logs[[s]]$pve_se, numeric(1)))
write_csv(summ, file.path(DIR_CLEAN, "gwas_summary.csv"))

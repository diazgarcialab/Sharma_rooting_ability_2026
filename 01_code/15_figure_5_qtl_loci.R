# ==============================================================================
# 15_qtl_gene_context.R -- genes around each QTL: a supplementary table and a
# locus-zoom figure per QTL (local Manhattan + gene track).
#
#   Rscript 01_code/15_qtl_gene_context.R
#   WINDOW_KB=50 Rscript 01_code/15_qtl_gene_context.R
#
# Inputs
#   02_clean_data/gwas_results/gwas_*.assoc.txt      GEMMA output
#   00_raw_data/genomic_resources/VITVgdSC2_v1.0.gff3
#   00_raw_data/genomic_resources/VITVgdSC2_v1.0.annotation.txt
#
# Outputs
#   02_clean_data/qtl_genes.csv          one row per QTL x gene, for the supplement
#   02_clean_data/qtl_leads.csv          the QTL list the table was built from
#   03_figures/Figure5_locus_<trait>_chrNN_<Mb>.{pdf,png,tiff}
#
# ==============================================================================

CODE_DIR <- dirname(normalizePath(local({
  for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) f[1] else "01_code/15_figure_5_qtl_loci.R"
})))
source(file.path(CODE_DIR, "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(tibble)
  library(ggplot2); library(patchwork)
})

WINDOW_KB   <- as.numeric(Sys.getenv("WINDOW_KB", unset = "100"))
# The TABLE window and the FIGURE window are deliberately different. +/-100 kb is
# the conservative net for listing candidate genes, but at that width a 3 kb gene
# occupies 1.5% of the panel and its exon structure is invisible -- the models
# collapse into featureless blocks. The figure therefore zooms to the interval
# where LD is still appreciable, which is both the more defensible interval and
# the one at which a gene model can actually be read.
FIG_WINDOW_KB <- as.numeric(Sys.getenv("FIG_WINDOW_KB", unset = "30"))
LD_CORE_KB  <- as.numeric(Sys.getenv("LD_CORE_KB", unset = "25"))
# How many genes get a name on the panel. Every gene inside the LD core used to
# be labelled, which at the smaller type was fine; at the larger type the names
# collide with each other and with the models. The nearest few to the lead SNP
# are the ones a reader looks for, and the rest are in the table.
MAX_GENE_LABELS <- as.integer(Sys.getenv("MAX_GENE_LABELS", unset = "4"))
ALPHA       <- 0.05
THRESHOLD   <- tolower(Sys.getenv("THRESHOLD", unset = "bonferroni"))
CLUMP_MB    <- as.numeric(Sys.getenv("CLUMP_MB", unset = "5"))
MAX_PER_CHR <- as.integer(Sys.getenv("MAX_PER_CHR", unset = "2"))

DIR_RES  <- file.path(DIR_RAW, "genomic_resources")
DIR_GW   <- file.path(DIR_CLEAN, "gwas_results")
GFF      <- file.path(DIR_RES, "VITVgdSC2_v1.0.gff3")
ANN      <- file.path(DIR_RES, "VITVgdSC2_v1.0.annotation.txt")
CACHE_G  <- file.path(DIR_CLEAN, "genes_hap1.tsv")
CACHE_A  <- file.path(DIR_CLEAN, "gene_annotation.tsv")

PT <- 0.85; SIG <- "#C1272D"; GREY <- "#7A7A7A"; ORANGE <- "#D97D2E"

# Genes inside the LD core are coloured by trait. There is no legend: each panel
# names its trait in the title and every core gene carries its own ID, so colour
# is a secondary cue rather than the only way to tell panels apart. That matters
# for the palette check -- the worst adjacent pair separates at deltaE 8.0 under
# deuteranopia, which sits on the floor and is only acceptable BECAUSE of that
# redundant labelling. Genes outside the core stay grey in every panel.
# Exons are drawn in the trait colour and the rest of the gene body in a washed
# version of the same hue, so intron/exon structure is visible without a legend.
lighten <- function(col, f = 0.72) {
  m <- grDevices::col2rgb(col)
  grDevices::rgb(t(255 - (255 - m) * (1 - f)), maxColorValue = 255)
}

TRAIT_GENE_COLS <- c(callus_stage = "#1F4E79", post_transplant = "#C1272D",
                     root_weight = "#0F7B6C", cutting_weight = "#8B4A9C",
                     transplant_response = "#B07A00")
TRAIT_LABELS <- c(callus_stage = "Callus-stage rooting",
                  post_transplant = "Post-transplant rooting",
                  root_weight = "Root weight", cutting_weight = "Cutting weight",
                  transplant_response = "Transplant response index")

# ---- 1. cache the two big files down to what is needed ------------------------
# The GFF is ~300 MB and the annotation ~67 MB. Reading either into R directly is
# slow and pointless: only `gene` rows and one description per gene are used, so
# awk reduces them once and the result is reused.
if (!file.exists(CACHE_G)) {
  stopifnot(file.exists(GFF))
  message("extracting gene models from the GFF3 (once) ...")
  system(sprintf(
    "awk -F'\\t' 'BEGIN{OFS=\"\\t\"; print \"chr\",\"start\",\"end\",\"strand\",\"gene\"} $3==\"gene\" && $1 ~ /hap1\\.chr/ { id=$9; sub(/^.*ID=/,\"\",id); sub(/;.*$/,\"\",id); print $1,$4,$5,$7,id }' %s > %s",
    shQuote(GFF), shQuote(CACHE_G)))
}
if (!file.exists(CACHE_A)) {
  stopifnot(file.exists(ANN))
  message("extracting gene descriptions from the annotation (once) ...")
  # column 2 is the hap1 transcript, column 6 the description; one row per gene
  system(sprintf(
    "awk -F'\\t' 'BEGIN{OFS=\"\\t\"; print \"gene\",\"description\"} $2 ~ /hap1\\.chr/ { g=$2; sub(/\\.t[0-9]+$/,\"\",g); if (!(g in seen) && $6 != \"\") { seen[g]=1; print g,$6 } }' %s > %s",
    shQuote(ANN), shQuote(CACHE_A)))
}
genes <- read_tsv(CACHE_G, show_col_types = FALSE) |>
  mutate(chr_n = as.integer(str_match(chr, "chr0*([0-9]+)")[, 2]),
         gene_short = sub("^.*\\.(g[0-9]+)$", "\\1", gene))
annot <- read_tsv(CACHE_A, show_col_types = FALSE)
message(sprintf("gene models: %s on %d chromosomes | annotated: %s",
                format(nrow(genes), big.mark = ","), n_distinct(genes$chr_n),
                format(nrow(annot), big.mark = ",")))

# ---- 2. lead SNPs -------------------------------------------------------------
assoc_files <- list.files(DIR_GW, pattern = "^gwas_.*\\.assoc\\.txt$", full.names = TRUE)
if (!length(assoc_files)) stop("no gwas_*.assoc.txt in ", DIR_GW)
slug_of <- function(f) sub("_(blup|blue)_PCA[0-9]+$", "",
                           sub("\\.assoc\\.txt$", "", sub("^gwas_", "", basename(f))))

read_assoc <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE))
    tibble::as_tibble(data.table::fread(path, showProgress = FALSE))
  else read_tsv(path, col_types = cols(.default = col_character()), progress = FALSE) |>
    mutate(across(c(ps, af, beta, se, p_wald), as.numeric))
}

leads <- list(); assoc_cache <- list()
for (f in assoc_files) {
  sl <- slug_of(f)
  a <- read_assoc(f) |> filter(is.finite(p_wald), p_wald > 0, p_wald <= 1) |>
    mutate(chr_n = as.integer(str_match(chr, "(?i)chr0*([0-9]+)")[, 2]))
  a$adj <- if (THRESHOLD == "fdr") p.adjust(a$p_wald, "BH") else pmin(a$p_wald * nrow(a), 1)
  thr <- if (THRESHOLD == "bonferroni") ALPHA / nrow(a) else
         if (any(a$adj <= ALPHA)) max(a$p_wald[a$adj <= ALPHA]) else NA_real_
  assoc_cache[[sl]] <- list(a = a, thr = thr)
  sg <- a |> filter(adj <= ALPHA, !is.na(chr_n)) |> arrange(p_wald)
  if (!nrow(sg)) { message(sprintf("  %-20s no significant SNP", sl)); next }
  keep <- logical(nrow(sg))
  ncr <- setNames(integer(length(unique(sg$chr_n))), as.character(unique(sg$chr_n)))
  for (j in seq_len(nrow(sg))) {
    if (any(keep & sg$chr_n == sg$chr_n[j] & abs(sg$ps - sg$ps[j]) < CLUMP_MB * 1e6)) next
    cj <- as.character(sg$chr_n[j])
    if (ncr[[cj]] >= MAX_PER_CHR) next
    keep[j] <- TRUE; ncr[[cj]] <- ncr[[cj]] + 1L
  }
  leads[[length(leads) + 1]] <- sg[keep, ] |>
    transmute(trait = sl, chr, chr_n, ps, p_wald, beta, se, af, thr = thr)
  message(sprintf("  %-20s %d significant SNP(s), %d QTL after clumping",
                  sl, nrow(sg), sum(keep)))
}
leads <- bind_rows(leads)
if (!nrow(leads)) stop("no significant QTL to annotate")
write_csv(leads, file.path(DIR_CLEAN, "qtl_leads.csv"))

# ---- 3. genes in the window ---------------------------------------------------
W  <- WINDOW_KB * 1000
WF <- FIG_WINDOW_KB * 1000
qtl_genes <- leads |>
  rowwise() |>
  do({
    ld <- .
    g <- genes |>
      filter(chr_n == ld$chr_n, end >= ld$ps - W, start <= ld$ps + W) |>
      mutate(trait = ld$trait, lead_ps = ld$ps, lead_p = ld$p_wald,
             # 0 when the lead SNP falls inside the gene
             distance_bp = pmax(0, pmax(start - ld$ps, ld$ps - end)),
             position = case_when(ld$ps >= start & ld$ps <= end ~ "within",
                                  ld$ps < start ~ "downstream of SNP",
                                  TRUE ~ "upstream of SNP"))
    g
  }) |>
  ungroup() |>
  left_join(annot, by = "gene") |>
  arrange(trait, chr_n, lead_ps, distance_bp) |>
  transmute(trait, chromosome = chr_n, lead_snp_bp = lead_ps, lead_p,
            # gene_full FIRST: inside transmute a later expression sees the
            # newly created columns, so writing `gene = gene_short` before
            # `gene_full = gene` silently set gene_full to the short name -- which
            # broke the exon join and put a duplicated ID in the table.
            gene_full = gene, gene = gene_short, start, end, strand,
            distance_bp, within_ld_core = distance_bp <= LD_CORE_KB * 1000,
            position, description)
write_csv(qtl_genes, file.path(DIR_CLEAN, "qtl_genes.csv"))

cat("\n", strrep("=", 76), "\n", sep = "")
cat(sprintf("GENES WITHIN +/-%g kb OF %d QTL\n", WINDOW_KB, nrow(leads)))
cat(strrep("=", 76), "\n")
smry <- qtl_genes |>
  group_by(trait, chromosome, lead_snp_bp) |>
  summarise(genes = n(), annotated = sum(!is.na(description)),
            within_core = sum(within_ld_core), .groups = "drop")
print(as.data.frame(smry), row.names = FALSE)

# ---- 3b. exon structure, for the windows only ---------------------------------
# Gene models are drawn as exon boxes on an intron line rather than as a plain
# arrow. That needs exon rows, and pulling every exon in the genome into R would
# be wasteful, so awk is given the window list and returns only what is on
# screen. If the extraction yields nothing the drawing falls back to a single
# block per gene, so the figure never depends on it.
exons <- tibble(chr_n = integer(), start = numeric(), end = numeric(),
                gene = character())
if (file.exists(GFF) && nrow(leads)) {
  reg_file <- tempfile(); ex_file <- tempfile()
  leads |> transmute(chr = paste0("VITVgdSC2_v1.0.hap1.chr", sprintf("%02d", chr_n)),
                     lo = ps - W, hi = ps + W) |>
    write_tsv(reg_file, col_names = FALSE)
  system(sprintf(
    "awk -F'\t' 'NR==FNR{n++; c[n]=$1; lo[n]=$2; hi[n]=$3; next} $3==\"exon\" { for(i=1;i<=n;i++) if ($1==c[i] && $5>=lo[i] && $4<=hi[i]) { id=$9; sub(/^.*Parent=/,\"\",id); sub(/;.*$/,\"\",id); sub(/\\.t[0-9]+$/,\"\",id); print $1\"\t\"$4\"\t\"$5\"\t\"id; break } }' %s %s > %s",
    shQuote(reg_file), shQuote(GFF), shQuote(ex_file)))
  if (file.exists(ex_file) && file.size(ex_file) > 0) {
    exons <- read_tsv(ex_file, col_names = c("chr", "start", "end", "gene"),
                      show_col_types = FALSE) |>
      mutate(chr_n = as.integer(str_match(chr, "chr0*([0-9]+)")[, 2])) |>
      distinct(chr_n, start, end, gene)
  }
  message(sprintf("exon rows in the QTL windows: %s",
                  format(nrow(exons), big.mark = ",")))
}

# ---- 4. locus zoom per QTL ----------------------------------------------------
# A proper gene model: a thin intron line spanning the locus, exon boxes on top
# of it, and a solid head at the 3' end so strand is readable without a legend.
# Genes with no exon rows fall back to one block, which still carries a head.
gene_track <- function(g, xlim, lead_bp, ex = NULL, core_col = SIG) {
  if (!nrow(g)) return(ggplot() + theme_void())
  span <- diff(xlim)
  g <- g |> arrange(start)
  # stack into rows so overlapping models never draw on top of each other
  occupied <- numeric(0); g$row <- 1L
  for (i in seq_len(nrow(g))) {
    r <- 1L
    while (r <= length(occupied) && occupied[r] > g$start[i] - span * 0.015) r <- r + 1L
    g$row[i] <- r; occupied[r] <- g$end[i]
  }
  g <- g |> mutate(core = distance_bp <= LD_CORE_KB * 1000,
                   s_c = pmax(start, xlim[1]), e_c = pmin(end, xlim[2]))

  # Gene models follow the genome-browser idiom: a thin INTRON LINE spanning the
  # gene with EXON BOXES on top, small chevrons along the line for strand, a
  # solid head only where the true 3' end is inside the window, and a zig-zag
  # break where the gene runs past the panel edge.
  #
  # The previous drawing was a pale pentagon body with darker exons on top. Two
  # things went wrong with it, both reported in review. (1) The body was
  # lighten(col, 0.72), which at print size is close enough to the panel
  # background that a long multi-exon gene read as a row of disconnected ticks
  # with no indication they belonged to one gene (g248770, chr15: 17 exons, of
  # which 15 fall in the window). (2) A gene wider than the window was simply
  # cut at the edge with no marker, so there was no way to tell a truncated gene
  # from one that genuinely ends there. The line ties the exons together and the
  # zig-zag says "continues"; neither costs any vertical space.
  HH <- 0.26                                     # half-height of an exon box
  g$w   <- g$e_c - g$s_c
  g$hd  <- pmin(0.32 * g$w, span * 0.020)        # head length in bp
  g$clipL <- g$start < xlim[1]
  g$clipR <- g$end   > xlim[2]

  # strand chevrons along the intron line
  tick <- span * 0.006
  chev <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    s0 <- g$s_c[i]; e0 <- g$e_c[i]
    if (e0 - s0 < span * 0.03) return(NULL)
    at <- seq(s0, e0, length.out = max(2, round((e0 - s0) / (span * 0.05)) + 1))
    at <- at[-c(1, length(at))]
    if (!length(at)) return(NULL)
    d <- if (g$strand[i] != "-") 1 else -1
    do.call(rbind, lapply(at, function(x) data.frame(
      x = x - d * tick / 2, xend = x + d * tick / 2,
      y = c(g$row[i] - HH * 0.5, g$row[i] + HH * 0.5),
      yend = g$row[i], core = g$core[i])))
  }))

  # terminal head, only when the real 3' end is on screen
  tip <- g |> mutate(vis3 = ifelse(strand != "-", !clipR, !clipL)) |> filter(vis3)
  tip_poly <- if (nrow(tip)) do.call(rbind, lapply(seq_len(nrow(tip)), function(i) {
    fwd <- tip$strand[i] != "-"; h <- span * 0.009
    x0 <- if (fwd) tip$e_c[i] else tip$s_c[i]; d <- if (fwd) 1 else -1
    data.frame(id = i, core = tip$core[i], x = c(x0, x0, x0 + d * h),
               y = c(tip$row[i] - HH, tip$row[i] + HH, tip$row[i]))
  })) else NULL

  # zig-zag break where the model is truncated by the window
  tr <- rbind(g |> filter(clipR) |> transmute(x = e_c, row, core, d = 1),
              g |> filter(clipL) |> transmute(x = s_c, row, core, d = -1))
  trunc_path <- if (nrow(tr)) do.call(rbind, lapply(seq_len(nrow(tr)), function(k) {
    wz <- span * 0.0045; y0 <- tr$row[k]
    data.frame(id = k, core = tr$core[k],
               x = tr$x[k] - tr$d[k] * c(0, wz, 0, wz, 0) * 0.9,
               y = y0 + c(-HH * 1.15, -HH * 0.55, 0, HH * 0.55, HH * 1.15))
  })) else NULL

  ex_use <- NULL
  if (!is.null(ex) && nrow(ex)) {
    ex_use <- ex |> filter(gene %in% g$gene_full) |>
      left_join(g |> select(gene_full, row, core, s_c, e_c, strand),
                by = c("gene" = "gene_full")) |>
      # Clip to the gene only: nothing is carved out of the model any more, so
      # no exon can be clipped out of existence (this dropped exon 2 of g248760).
      mutate(start = pmax(start, s_c), end = pmin(end, e_c),
             core_ex = ifelse(core, "exon_core", "exon_other")) |>
      filter(end > start)
  }
  # a gene with no exon rows still needs a body, or it would be a bare line
  no_ex <- if (is.null(ex_use)) g else g |> filter(!(gene_full %in% ex_use$gene))

  pg <- ggplot() +
    annotate("rect", xmin = lead_bp - LD_CORE_KB * 1000,
             xmax = lead_bp + LD_CORE_KB * 1000, ymin = -Inf, ymax = Inf,
             fill = ORANGE, alpha = 0.07) +
    geom_segment(data = g, aes(x = s_c, xend = e_c, y = row, yend = row,
                               colour = core), linewidth = 0.4) +
    { if (!is.null(chev) && nrow(chev))
        geom_segment(data = chev, aes(x = x, xend = xend, y = y, yend = yend,
                                      colour = core),
                     linewidth = 0.28, lineend = "round") } +
    { if (nrow(no_ex))
        geom_rect(data = no_ex, aes(xmin = s_c, xmax = e_c, ymin = row - HH,
                                    ymax = row + HH, fill = core), colour = NA) } +
    { if (!is.null(ex_use) && nrow(ex_use))
        geom_rect(data = ex_use, aes(xmin = start, xmax = end,
                                     ymin = row - HH, ymax = row + HH,
                                     fill = core_ex), colour = NA) } +
    { if (!is.null(tip_poly))
        geom_polygon(data = tip_poly, aes(x, y, group = id, fill = core),
                     colour = NA) } +
    { if (!is.null(trunc_path))
        geom_path(data = trunc_path, aes(x, y, group = id, colour = core),
                  linewidth = 0.5, lineend = "round") } +
    scale_fill_manual(values = c(`TRUE` = core_col, `FALSE` = GREY,
                                 exon_core = core_col, exon_other = GREY),
                      guide = "none") +
    scale_colour_manual(values = c(`TRUE` = core_col, `FALSE` = GREY), guide = "none") +
    scale_x_continuous(limits = xlim, expand = c(0, 0),
                       labels = \(v) sprintf("%.2f", v / 1e6)) +
    scale_y_reverse(expand = expansion(add = 1.05)) +
    labs(x = NULL, y = NULL) +
    theme_paper() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          # Horizontal rules behind the gene models read as if the rows meant
          # something, so those go. The VERTICAL rules stay: with the same breaks
          # on both panels they run from the Manhattan straight down to the axis,
          # which is what ties a peak to the gene under it.
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(colour = INK[["faint"]], linewidth = 0.2),
          panel.border = element_rect(fill = NA, colour = INK[["faint"]],
                                      linewidth = 0.3))

  lab <- g |> filter(core) |> mutate(mid = (s_c + e_c) / 2) |>
    arrange(distance_bp) |> head(MAX_GENE_LABELS)
  if (nrow(lab)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      pg <- pg + ggrepel::geom_text_repel(
        data = lab, aes(mid, row, label = gene_short), size = 2.15, seed = 3,
        colour = INK[["primary"]], segment.colour = INK[["muted"]],
        segment.size = 0.18, min.segment.length = 0, box.padding = 0.16,
        max.overlaps = Inf, direction = "both", nudge_y = 0.70,
        point.padding = 0.22)
    } else {
      pg <- pg + geom_text(data = lab, aes(mid, row, label = gene_short),
                           size = 2.15, vjust = -1, colour = INK[["primary"]])
    }
  }
  pg
}

dir.create(DIR_FIG, showWarnings = FALSE)
# One figure per QTL is redundant now that they are all in one multipanel figure;
# INDIVIDUAL=1 brings them back for talks or a closer look at a single locus.
INDIVIDUAL <- tolower(Sys.getenv("INDIVIDUAL", unset = "false")) %in% c("1","true","yes")
made <- character(0)
for (i in if (INDIVIDUAL) seq_len(nrow(leads)) else integer(0)) {
  ld <- leads[i, ]
  ac <- assoc_cache[[ld$trait]]
  reg <- ac$a |> filter(chr_n == ld$chr_n, ps >= ld$ps - WF, ps <= ld$ps + WF) |>
    mutate(logp = -log10(p_wald), sig = adj <= ALPHA)
  if (nrow(reg) < 5) next
  xlim <- c(ld$ps - WF, ld$ps + WF)

  p_top <- ggplot(reg, aes(ps, logp)) +
    annotate("rect", xmin = ld$ps - LD_CORE_KB * 1000, xmax = ld$ps + LD_CORE_KB * 1000,
             ymin = -Inf, ymax = Inf, fill = ORANGE, alpha = 0.07) +
    geom_point(size = PT, stroke = 0, colour = GREY) +
    geom_point(data = filter(reg, sig), aes(ps, logp), shape = 21, fill = ORANGE,
               colour = INK[["primary"]], size = PT * 1.45, stroke = 0.35) +
    geom_point(data = ld, aes(ps, -log10(p_wald)), shape = 23, fill = SIG,
               colour = INK[["primary"]], size = 1.9, stroke = 0.4) +
    { if (is.finite(ac$thr)) geom_hline(yintercept = -log10(ac$thr), colour = SIG,
                                        linewidth = 0.4, linetype = "22") } +
    scale_x_continuous(limits = xlim, expand = c(0, 0),
                       labels = \(v) sprintf("%.2f", v / 1e6)) +
    labs(x = NULL, y = expression(-log[10] * italic(p)),
         title = sprintf("%s  |  chr%02d: %.3f Mb",
                         TRAIT_LABELS[ld$trait] %||% ld$trait, ld$chr_n, ld$ps / 1e6)) +
    theme_paper() +
    theme(plot.title = element_text(size = rel(0.95)),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())

  g <- qtl_genes |> filter(trait == ld$trait, chromosome == ld$chr_n,
                           lead_snp_bp == ld$ps,
                           end >= ld$ps - WF, start <= ld$ps + WF) |>
    rename(gene_short = gene)
  p_bot <- gene_track(g, xlim, ld$ps, ex = exons |> filter(chr_n == ld$chr_n),
                      core_col = unname(TRAIT_GENE_COLS[ld$trait]) %||% SIG) +
    labs(x = sprintf("Position on chromosome %d (Mb)", ld$chr_n))

  fig <- p_top / p_bot + plot_layout(heights = c(1, 0.45))
  nm <- sprintf("Figure5_locus_%s_chr%02d_%.2fMb", ld$trait, ld$chr_n, ld$ps / 1e6)
  save_figure(fig, nm, width = 3.6, height = 3.0)
  made <- c(made, nm)
}
cat(sprintf("\n%d locus-zoom figures written to %s\n", length(made), DIR_FIG))

# ---- 5. one multipanel figure with every QTL ---------------------------------
# Space is saved by dropping repeated furniture rather than by shrinking panels:
# the y-axis title appears once per row, the x-axis title once per column, and
# each locus carries a compact one-line title instead of a title plus a legend.
NCOL   <- as.integer(Sys.getenv("QTL_NCOL", unset = "3"))
# Height given to one row of panels. Panels are never squeezed to fit a page --
# if there are more rows than fit, the figure is split into parts instead.
ROW_IN <- as.numeric(Sys.getenv("QTL_ROW_IN", unset = "1.55"))
PAGE_IN <- 9.4

locus_panel <- function(ld, show_y, show_x) {
  ac <- assoc_cache[[ld$trait]]
  reg <- ac$a |> filter(chr_n == ld$chr_n, ps >= ld$ps - WF, ps <= ld$ps + WF) |>
    mutate(logp = -log10(p_wald), sig = adj <= ALPHA)
  if (nrow(reg) < 5) return(NULL)
  xlim <- c(ld$ps - WF, ld$ps + WF)
  ymax <- max(reg$logp, na.rm = TRUE)

  # The lead marker sits at the centre of every panel by construction, so its
  # coordinate is carried by a repelled label on the point itself rather than by
  # a subtitle and an axis title -- two lines of furniture per panel, times
  # eighteen panels, is most of a row of plots.
  lead_lab <- tibble(ps = ld$ps, logp = -log10(ld$p_wald),
                     txt = sprintf("chr%02d: %.3f Mb", ld$chr_n, ld$ps / 1e6))

  pt <- ggplot(reg, aes(ps, logp)) +
    annotate("rect", xmin = ld$ps - LD_CORE_KB * 1000,
             xmax = ld$ps + LD_CORE_KB * 1000, ymin = -Inf, ymax = Inf,
             fill = ORANGE, alpha = 0.07) +
    geom_point(size = PT * 1.7, stroke = 0, colour = GREY) +
    geom_point(data = filter(reg, sig), aes(ps, logp), shape = 21, fill = ORANGE,
               colour = INK[["primary"]], size = PT * 2.4, stroke = 0.56) +
    geom_point(data = ld, aes(ps, -log10(p_wald)), shape = 23,
               fill = unname(TRAIT_GENE_COLS[ld$trait]) %||% SIG,
               colour = INK[["primary"]], size = 2.7, stroke = 0.62) +
    { if (is.finite(ac$thr)) geom_hline(yintercept = -log10(ac$thr), colour = SIG,
                                        linewidth = 0.42, linetype = "22") } +
    # trait name inside the panel, top left
    annotate("text", x = xlim[1] + diff(xlim) * 0.02, y = ymax * 1.27,
             label = TRAIT_LABELS[ld$trait] %||% ld$trait, hjust = 0, vjust = 0.4,
             size = 2.45, fontface = "bold",
             colour = unname(TRAIT_GENE_COLS[ld$trait]) %||% INK[["primary"]]) +
    # The x expansion MUST match the gene panel below (expansion(mult = 0.04)).
    # When this was expand = c(0, 0) and the gene panel was 0.04, the same data
    # coordinates landed at different x positions in the two panels, so the
    # shaded +/-25 kb band and the vertical gridlines were visibly offset
    # between the Manhattan and the gene track (reported by NC).
    scale_x_continuous(limits = xlim, expand = expansion(mult = 0.04)) +
    # Just enough headroom for the two labels and no more. Earlier versions grew
    # this to 0.55 to guarantee the labels never collided, which bought the space
    # at the cost of a band of empty panel above every peak.
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.32)),
                       breaks = scales::pretty_breaks(3)) +
    labs(x = NULL, y = if (show_y) expression(-log[10] * italic(p)) else NULL) +
    theme_paper() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = rel(1.05)),
          axis.title.y = element_text(size = rel(1.0)),
          panel.border = element_rect(fill = NA, colour = INK[["faint"]],
                                      linewidth = 0.3),
          plot.margin = margin(2, 3, 0, 3))

  # The position label is placed deterministically rather than with ggrepel.
  # Three repel configurations were tried on the same panels: free placement
  # pulled the label onto the trait name, constraining it below the peak put it
  # on the threshold line, and pinning it to the right edge pushed it off the
  # panel with only its leader line visible. The geometry here is fixed and
  # known -- the lead marker is always dead centre and the trait name always
  # occupies the top left -- so the label is anchored to the top right and a
  # short leader is drawn to the marker.
  # Second line, right-aligned: sharing one line with the trait name works until
  # a wide position ("chr02: 19.088 Mb") meets a long trait name, and then they
  # touch. Both fit inside the existing headroom when stacked, so no extra
  # vertical space is needed.
  lab_x <- xlim[2] - diff(xlim) * 0.02
  lab_y <- ymax * 1.11
  pt <- pt +
    geom_segment(data = lead_lab,
                 aes(x = xlim[2] - diff(xlim) * 0.32, y = lab_y * 0.99,
                     xend = ps + diff(xlim) * 0.02, yend = logp + ymax * 0.02),
                 colour = INK[["muted"]], linewidth = 0.22) +
    geom_text(data = lead_lab, aes(x = lab_x, y = lab_y, label = txt),
              hjust = 1, vjust = 0.5, size = 2.45, colour = INK[["primary"]])

  g <- qtl_genes |> filter(trait == ld$trait, chromosome == ld$chr_n,
                           lead_snp_bp == ld$ps,
                           end >= ld$ps - WF, start <= ld$ps + WF) |>
    rename(gene_short = gene)
  pb <- gene_track(g, xlim, ld$ps, ex = exons |> filter(chr_n == ld$chr_n),
                   core_col = unname(TRAIT_GENE_COLS[ld$trait]) %||% SIG) +
    # No axis title: the marker position is on the point already, and the unit
    # rides on the last tick instead of a whole label line.
    # The unit rides on each tick and there is no axis title at all. Two breaks
    # rather than three, and a small x expansion, keep the outermost label from
    # running off the panel -- which is what clipped "1.24 Mb" to "1.24 M" when
    # the labels were wider than the space at the panel edge.
    # Breaks are placed relative to the lead SNP rather than by pretty(): a
    # pretty break can land on the panel edge, where the widened "19.12Mb" label
    # is clipped to "19.12M".
    scale_x_continuous(limits = xlim, expand = expansion(mult = 0.04),
                       breaks = ld$ps + c(-0.62, 0, 0.62) * WF,
                       labels = \(v) paste0(sprintf("%.2f", v / 1e6), "Mb")) +
    labs(x = NULL) +
    theme(axis.text.x = element_text(size = rel(1.05)),
          plot.margin = margin(0, 3, 2, 3))

  wrap_elements(full = pt / pb + plot_layout(heights = c(1, 0.58)))
}

n_q <- nrow(leads)
rows_per_page <- max(1L, floor(PAGE_IN / ROW_IN))
per_page <- rows_per_page * NCOL
pages <- split(seq_len(n_q), ceiling(seq_len(n_q) / per_page))
for (pi in seq_along(pages)) {
  ix <- pages[[pi]]
  nrow_p <- ceiling(length(ix) / NCOL)
  cells <- lapply(seq_along(ix), function(k)
    locus_panel(leads[ix[k], ],
                show_y = ((k - 1) %% NCOL) == 0,
                show_x = k > (nrow_p - 1) * NCOL))
  cells <- cells[!vapply(cells, is.null, logical(1))]
  if (!length(cells)) next
  fig <- wrap_plots(cells, ncol = NCOL)
  nm <- if (length(pages) > 1) sprintf("Figure5_part%d", pi) else
        "Figure5"
  save_figure(fig, nm, width = 7.28, height = ROW_IN * nrow_p + 0.15)
}
cat(sprintf("\nmultipanel figure: %d QTL, %d column(s), %d page(s)\n",
            n_q, NCOL, length(pages)))

# ---- 6. the supplementary table ----------------------------------------------
tab <- qtl_genes |>
  left_join(leads |> transmute(trait, chromosome = chr_n, lead_snp_bp = ps,
                               lead_beta = beta, lead_se = se, lead_af = af),
            by = c("trait", "chromosome", "lead_snp_bp")) |>
  transmute(Trait = TRAIT_LABELS[trait] %||% trait,
            Chromosome = chromosome,
            `Lead SNP (bp)` = lead_snp_bp,
            `Lead p` = signif(lead_p, 3),
            `Effect (beta)` = round(lead_beta, 4),
            `SE` = round(lead_se, 4),
            `Allele freq` = round(lead_af, 3),
            Gene = gene, `Gene ID` = gene_full,
            `Start (bp)` = start, `End (bp)` = end, Strand = strand,
            `Distance to lead SNP (kb)` = round(distance_bp / 1000, 2),
            `Within LD core` = ifelse(within_ld_core, "yes", "no"),
            `Relative position` = position,
            `Functional annotation` = ifelse(is.na(description),
                                             "no functional annotation", description))
write_csv(tab, file.path(DIR_CLEAN, "table_s_qtl_genes.csv"))
xl <- file.path(DIR_CLEAN, "table_s_qtl_genes.xlsx")
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(tab, xl); message("wrote ", xl)
} else if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(tab, xl, asTable = TRUE); message("wrote ", xl)
} else {
  message("note: install.packages(\"writexl\") to also get an .xlsx; the CSV is written either way")
}
cat(sprintf("supplementary table: %d rows, %d QTL -> table_s_qtl_genes.csv\n",
            nrow(tab), n_distinct(paste(tab$Trait, tab$Chromosome, tab$`Lead SNP (bp)`))))
message("done")

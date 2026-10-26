# ==============================================================================
# 00_setup.R
# Project-wide configuration: paths, packages, plotting theme, species look-up.
# Sourced by every analysis / figure script. Contains no analysis itself.
#
# Project : Harnessing Vitis germplasm diversity to dissect and predict
#           adventitious rooting traits in grapevine
# ==============================================================================

# ---- 1. Paths ----------------------------------------------------------------
# PROJ_DIR is the "final_paper" folder; everything else hangs off it.
# Resolution order:
#   1. PROJ_DIR already defined in the calling environment
#   2. environment variable ROOTING_PROJ_DIR
#   3. auto-detected from the location of this file (Rscript, source(), RStudio)
#   4. the working directory or its parent, whichever contains 00_raw_data/
.locate_this_file <- function() {
  for (i in seq_len(sys.nframe())) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) return(normalizePath(o, mustWork = FALSE))
  }
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], mustWork = FALSE))
  NA_character_
}

if (!exists("PROJ_DIR") || !dir.exists(file.path(PROJ_DIR, "00_raw_data"))) {
  PROJ_DIR <- Sys.getenv("ROOTING_PROJ_DIR", unset = NA_character_)
  if (is.na(PROJ_DIR) || !dir.exists(file.path(PROJ_DIR, "00_raw_data"))) {
    .this <- .locate_this_file()                       # .../01_code/00_setup.R
    PROJ_DIR <- if (!is.na(.this)) dirname(dirname(.this)) else NA_character_
  }
  if (is.na(PROJ_DIR) || !dir.exists(file.path(PROJ_DIR, "00_raw_data"))) {
    PROJ_DIR <- if (dir.exists("00_raw_data")) normalizePath(".") else normalizePath("..")
  }
}
if (!dir.exists(file.path(PROJ_DIR, "00_raw_data"))) {
  stop("Could not locate the project folder. Set it explicitly, e.g.\n",
       '  PROJ_DIR <- "~/Documents/ucdavis/rooting_ability/final_paper"')
}
PROJ_DIR <- normalizePath(PROJ_DIR)

DIR_RAW    <- file.path(PROJ_DIR, "00_raw_data")
DIR_CODE   <- file.path(PROJ_DIR, "01_code")
DIR_CLEAN  <- file.path(PROJ_DIR, "02_clean_data")
DIR_FIG    <- file.path(PROJ_DIR, "03_figures")
DIR_CACHE  <- file.path(PROJ_DIR, "00_raw_data", "spatial_cache")  # WorldClim / GADM

for (d in c(DIR_CLEAN, DIR_FIG, DIR_CACHE)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- 2. Packages -------------------------------------------------------------
.pkgs <- c("dplyr", "tidyr", "readr", "stringr", "tibble", "ggplot2", "scales",
           "patchwork", "sf", "terra", "ragg")

.missing <- .pkgs[!vapply(.pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing)) {
  stop("Missing packages: ", paste(.missing, collapse = ", "),
       "\n  install.packages(c(", paste0('"', .missing, '"', collapse = ", "), "))")
}
invisible(lapply(.pkgs, library, character.only = TRUE))

sf::sf_use_s2(FALSE)      # planar ops are fine at this scale and much faster

# R ships a 60 s download timeout, which aborts the ~130 MB WorldClim archive
# partway through and surfaces as a confusing "cannot open URL" / corrupt-zip
# error. Raise it before anything tries to download.
options(timeout = max(1800, getOption("timeout")))
options(dplyr.summarise.inform = FALSE)
set.seed(20260818)        # deterministic jitter etc.

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 3. Study extent ---------------------------------------------------------
# Southwestern USA + northern Mexico, matching the germplasm collection range.
STUDY_EXTENT <- c(xmin = -125, xmax = -94, ymin = 25.5, ymax = 43)

# ---- 4. Categorical palette --------------------------------------------------
# Okabe-Ito derived, validated for CVD separation (deutan/protan/tritan) and for
# lightness/chroma against a light map background. Species identity is carried by
# colour x shape jointly, so each colour is reused at most 4 times but never
# within the same shape family.
PAL <- c(
  blue     = "#0072B2",
  vermilion= "#D55E00",
  green    = "#009E73",
  purple   = "#CC79A7",
  skyblue  = "#56B4E9",
  amber    = "#E69F00"
)

# Neutral ink ramp (text, axes, borders) -- never a series colour.
INK <- c(primary = "#1A1A1A", secondary = "#4D4D4D",
         muted   = "#8A8A8A", faint     = "#D9D9D9")

# Sequential ramps for the climate rasters. One hue each, light -> dark, low
# chroma at the light end so the categorical points stay on top of the stack.
RAMP_PRECIP <- c("#FBF7EC", "#DCEBE2", "#AFD7DA", "#79BBD0", "#4A93BF", "#2E639C", "#1F3B70")
RAMP_TEMP   <- c("#FEF6E3", "#FDE2B8", "#FAC38C", "#F09A64", "#DC6E4A", "#B8453C", "#7F2733")

# Diverging ramp for correlations: two hues either side of a neutral grey
# midpoint (never a hue at the middle, never a rainbow).
RAMP_DIV <- c("#1F3B70", "#4A93BF", "#AFD7DA", "#F0EFEE",
              "#FAC38C", "#DC6E4A", "#7F2733")

# Year is ordinal, so it gets a single-hue SEQUENTIAL ramp (light = earliest),
# not categorical hues -- and a hue distinct from the species palette so a year
# is never mistaken for a species. Checked for lightness monotonicity, which is
# the right criterion for a sequential scale.
YEAR_COLS <- c("2024" = "#9CC0DE", "2025" = "#3C7DBE", "2026" = "#153E66")

OCEAN_FILL  <- "#EDF1F4"
LAND_FILL   <- "#F2F0EC"   # land outside the raster footprint

# ---- 5. Species look-up ------------------------------------------------------
# `species` matches the values in 00_raw_data/BLUPs_by_year_and_trait.csv.
# `label` is the typeset name used in figures (italics applied at plot time).
# `group` = broad geographic origin of the species; it drives the point shape.
# `colour` is assigned within a group so that the four hues in any one shape
# family are maximally separated.
SPECIES_LUT <- tibble::tribble(
  ~species,             ~label,             ~group,                        ~colour,
  # -- Southwestern deserts & California ---------------------------------------
  "V. arizonica",       "V. arizonica",     "Southwest & California",      PAL[["blue"]],
  "V. girdiana",        "V. girdiana",      "Southwest & California",      PAL[["skyblue"]],
  "V. californica",     "V. californica",   "Southwest & California",      PAL[["green"]],
  "V. treleasei",       "V. treleasei",     "Southwest & California",      PAL[["purple"]],
  # -- Southern Great Plains ---------------------------------------------------
  "V. acerifolia",      "V. acerifolia",    "Southern Plains",             PAL[["vermilion"]],
  "V. doaniana",        "V. \u00d7 doaniana",    "Southern Plains",             PAL[["amber"]],
  "V. candicans",       "V. candicans",     "Southern Plains",             PAL[["purple"]],
  "V. champinii",       "V. champinii",     "Southern Plains",             PAL[["green"]],
  # -- Edwards Plateau & Ozarks ------------------------------------------------
  "V. berlandieri",     "V. berlandieri",   "Edwards Plateau & Ozarks",    PAL[["blue"]],
  "V. monticola",       "V. monticola",     "Edwards Plateau & Ozarks",    PAL[["skyblue"]],
  "V. rupestris",       "V. rupestris",     "Edwards Plateau & Ozarks",    PAL[["vermilion"]],
  "V. cinerea",         "V. cinerea",       "Edwards Plateau & Ozarks",    PAL[["green"]],
  # -- Northern & Eastern ------------------------------------------------------
  "V. riparia",         "V. riparia",       "Northern & Eastern",          PAL[["blue"]],
  "V. vulpina",         "V. vulpina",       "Northern & Eastern",          PAL[["skyblue"]],
  "V. aestivalis",      "V. aestivalis",    "Northern & Eastern",          PAL[["amber"]],
  "V. labrusca",        "V. labrusca",      "Northern & Eastern",          PAL[["purple"]],
  # -- Southeastern ------------------------------------------------------------
  "V. simpsonii",       "V. simpsonii",     "Southeastern",                PAL[["vermilion"]],
  "V. shuttleworthii",  "V. shuttleworthii","Southeastern",                PAL[["green"]]
)

# Shape per geographic group. Solid pch codes are used so the `fill` aesthetic
# stays free for the climate raster; points are drawn over a white halo layer,
# which gives the same "filled symbol with a white outline" look without needing
# a second fill scale (i.e. without ggnewscale).
GROUP_SHAPES <- c(
  "Southwest & California"   = 16,  # circle
  "Southern Plains"          = 17,  # triangle up
  "Edwards Plateau & Ozarks" = 15,  # square
  "Northern & Eastern"       = 18,  # diamond
  "Southeastern"             =  8   # star
)
# pch 18/17 render optically smaller/larger than a circle of the same cex.
SHAPE_SIZE_ADJ <- c("16" = 1.00, "17" = 1.05, "15" = 0.90, "18" = 1.30, "8" = 1.05)
GROUP_ORDER <- names(GROUP_SHAPES)

# Colour per geographic group, for panels where species sits on an axis and
# colouring by species would be redundant. Subset of the validated PAL.
GROUP_COLS <- c(
  "Southwest & California"   = PAL[["blue"]],
  "Southern Plains"          = PAL[["vermilion"]],
  "Edwards Plateau & Ozarks" = PAL[["green"]],
  "Northern & Eastern"       = PAL[["purple"]],
  "Southeastern"             = PAL[["amber"]]
)

SPECIES_LUT <- SPECIES_LUT |>
  mutate(group = factor(group, levels = GROUP_ORDER),
         shape = unname(GROUP_SHAPES[as.character(group)])) |>
  arrange(group, label)

SPECIES_LEVELS <- SPECIES_LUT$species
SPECIES_COLS   <- setNames(SPECIES_LUT$colour, SPECIES_LUT$species)
SPECIES_SHAPES <- setNames(SPECIES_LUT$shape,  SPECIES_LUT$species)
SPECIES_SIZEADJ <- setNames(unname(SHAPE_SIZE_ADJ[as.character(SPECIES_LUT$shape)]),
                            SPECIES_LUT$species)
SPECIES_LABELS <- setNames(SPECIES_LUT$label,  SPECIES_LUT$species)

# ---- 6. Plot theme -----------------------------------------------------------
BASE_SIZE <- 8   # pt; journal figures are read at ~100% of print size

theme_paper <- function(base_size = BASE_SIZE, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text             = element_text(colour = INK[["primary"]]),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(size = rel(1.05), face = "bold", hjust = 0,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = rel(0.95), colour = INK[["secondary"]],
                                      margin = margin(b = 6)),
      plot.tag         = element_text(size = rel(1.7), face = "bold"),
      plot.tag.position = c(0.012, 0.985),
      # axis titles and tick labels are both full-strength ink; the grey ramp is
      # reserved for grid lines, captions and other non-data furniture
      axis.title       = element_text(size = rel(0.95), colour = INK[["primary"]]),
      axis.text        = element_text(size = rel(0.85), colour = INK[["primary"]]),
      legend.text      = element_text(size = rel(0.85), colour = INK[["primary"]]),
      panel.grid.major = element_line(colour = INK[["faint"]], linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(size = rel(0.9), face = "bold"),
      legend.key.height = unit(9, "pt"),
      strip.text       = element_text(size = rel(0.95), face = "bold",
                                      colour = INK[["primary"]]),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

theme_map <- function(base_size = BASE_SIZE) {
  theme_paper(base_size) +
    theme(
      panel.background = element_rect(fill = OCEAN_FILL, colour = NA),
      panel.border     = element_rect(fill = NA, colour = INK[["muted"]],
                                      linewidth = 0.3),
      panel.grid.major = element_line(colour = alpha("white", 0.45), linewidth = 0.22),
      axis.title       = element_blank(),
      axis.text        = element_text(size = rel(0.8), colour = INK[["primary"]]),
      axis.ticks       = element_line(colour = INK[["secondary"]], linewidth = 0.25),
      axis.ticks.length = unit(1.5, "pt")
    )
}

theme_set(theme_paper())

# ---- 7. Small helpers --------------------------------------------------------

#' Save a figure at journal resolution in several formats.
save_figure <- function(plot, name, width, height, dir = DIR_FIG,
                        formats = c("pdf", "png", "tiff"), dpi = 600) {
  paths <- character(0)
  for (fmt in formats) {
    f <- file.path(dir, paste0(name, ".", fmt))
    if (fmt == "pdf") {
      ggsave(f, plot, width = width, height = height, units = "in",
             device = grDevices::cairo_pdf)
    } else if (fmt == "tiff") {
      ggsave(f, plot, width = width, height = height, units = "in", dpi = dpi,
             device = ragg::agg_tiff, compression = "lzw")
    } else {
      ggsave(f, plot, width = width, height = height, units = "in", dpi = dpi,
             device = ragg::agg_png, background = "white")
    }
    paths <- c(paths, f)
  }
  message("saved: ", paste(basename(paths), collapse = ", "))
  invisible(paths)
}

# ggplot2 changed its guide API in 3.5 and moved themes to S7 in 4.0, so a few
# calls below are version-guarded. patchwork < 1.3.2 cannot combine with an S7
# theme (its `&` operator fails to dispatch), which is why none of these scripts
# use `&`.
GG_VER <- utils::packageVersion("ggplot2")
PW_VER <- utils::packageVersion("patchwork")
GG35   <- GG_VER >= "3.5.0"
if (GG_VER >= "4.0.0" && PW_VER < "1.3.2") {
  message("note: ggplot2 ", GG_VER, " with patchwork ", PW_VER,
          " -- consider install.packages(\"patchwork\") for full compatibility")
}

#' Place a legend inside the panel, compatible with ggplot2 < 3.5 and >= 3.5.
legend_inside <- function(x, y, just = c(0, 0)) {
  if (GG35) {
    theme(legend.position = "inside",
          legend.position.inside = c(x, y),
          legend.justification = just)
  } else {
    theme(legend.position = c(x, y), legend.justification = just)
  }
}

#' Compact horizontal colourbar guide, ditto.
guide_bar <- function(title, width = 52, height = 4.5) {
  if (GG35) {
    guide_colourbar(theme = theme(
      legend.title.position = "top",
      legend.title    = element_text(hjust = 0),
      legend.key.width  = unit(width,  "pt"),
      legend.key.height = unit(height, "pt"),
      legend.ticks = element_line(colour = "white", linewidth = 0.4),
      legend.frame = element_blank()))
  } else {
    guide_colourbar(title.position = "top", title.hjust = 0,
                    barwidth = unit(width, "pt"), barheight = unit(height, "pt"),
                    ticks.colour = "white", ticks.linewidth = 0.4,
                    frame.colour = NA)
  }
}

#' Degree labels for map axes.
DEG <- "\u00b0"   # escaped so the scripts are locale-independent
lab_lon <- function(x) paste0(abs(x), DEG, ifelse(x < 0, "W", "E"))
lab_lat <- function(x) paste0(abs(x), DEG, ifelse(x < 0, "S", "N"))

message("setup complete  |  PROJ_DIR = ", PROJ_DIR)
message("  R ", getRversion(), " | ggplot2 ", GG_VER, " | patchwork ", PW_VER,
        " | sf ", utils::packageVersion("sf"),
        " | terra ", utils::packageVersion("terra"))

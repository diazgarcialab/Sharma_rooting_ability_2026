# ==============================================================================
# 02_Figure3.R
#
# FIGURE 1 -- Geographic distribution and climatic context of the Vitis panel.
#
#   Panel A : accessions over WorldClim 2.1 annual precipitation (bio12)
#   Panel B : accessions over WorldClim 2.1 annual mean temperature (bio1)
#   Legend  : hand-built panel, species grouped by geographic origin, with the
#             number of mapped accessions per species
#
# Species identity is encoded jointly by colour and shape: shape = broad
# geographic origin of the species, colour = species within that origin. The six
# hues are an Okabe-Ito subset validated for colour-vision-deficiency separation.
#
# Climate rasters are resolved in three tiers (see get_climate_panels() below):
#   1. a cached GeoTIFF in 00_raw_data/spatial_cache/
#   2. a fresh WorldClim 2.1 download via geodata (needs internet, once)
#   3. an interpolated stand-in built from the accession bioclim values, used
#      only so the figure still renders offline -- it is watermarked as a draft
#      and must NOT be used in the manuscript
#
# Input   02_clean_data/accessions_georeferenced.csv   (from 01_prepare_accessions.R)
# Output  03_figures/Figure3.{pdf,png,tiff}
# ==============================================================================

source(file.path(dirname(normalizePath(
  local({
    for (i in seq_len(sys.nframe())) { o <- sys.frame(i)$ofile; if (!is.null(o)) return(o) }
    a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) f[1] else "01_code/02_figure_3_accession_map.R"
  }))), "00_setup.R"))

WC_RES     <- 2.5      # WorldClim resolution in arc-minutes (2.5 or 10)
POINT_SIZE <- 1.55
JITTER_DEG <- 0.055    # tiny deterministic jitter; many accessions share a site

# ---- 1. Accessions -----------------------------------------------------------
acc_file <- file.path(DIR_CLEAN, "accessions_georeferenced.csv")
if (!file.exists(acc_file)) {
  source(file.path(DIR_CODE, "01_prepare_accessions.R"))
}
acc <- read_csv(acc_file, show_col_types = FALSE) |>
  mutate(species = factor(species, levels = SPECIES_LEVELS))

# Species actually represented on the map (three species have no georeferenced
# accession and are reported in the caption instead of shown as empty keys).
sp_n <- acc |> count(species, .drop = FALSE, name = "n_mapped")
sp_on_map <- SPECIES_LUT |>
  left_join(sp_n, by = "species") |>
  filter(n_mapped > 0) |>
  mutate(group = droplevels(group))

acc <- acc |>
  filter(species %in% sp_on_map$species) |>
  mutate(species = droplevels(species),
         size_adj = unname(SPECIES_SIZEADJ[as.character(species)]),
         lon_j = longitude + runif(dplyr::n(), -JITTER_DEG, JITTER_DEG),
         lat_j = latitude  + runif(dplyr::n(), -JITTER_DEG, JITTER_DEG))

message(sprintf("mapping %d accessions, %d species",
                nrow(acc), dplyr::n_distinct(acc$species)))

# ---- 2. Basemap --------------------------------------------------------------
bbox_sf <- sf::st_as_sfc(sf::st_bbox(c(STUDY_EXTENT), crs = 4326))

get_basemap <- function() {
  if (requireNamespace("rnaturalearthdata", quietly = TRUE)) {
    ctry  <- sf::st_as_sf(rnaturalearthdata::countries50)
    admin <- sf::st_as_sf(rnaturalearthdata::states50)
  } else if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    ctry  <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
    admin <- rnaturalearth::ne_states(country = c("United States of America", "Mexico"),
                                      returnclass = "sf")
  } else {
    stop("install.packages('rnaturalearthdata')  # basemap polygons")
  }
  admin_col <- intersect(c("admin", "iso_a2", "sov_a3"), names(admin))[1]
  admin <- admin[admin[[admin_col]] %in%
                   c("United States of America", "Mexico", "US", "MX", "USA", "MEX"), ]
  list(
    countries = sf::st_make_valid(sf::st_crop(sf::st_make_valid(ctry),  sf::st_bbox(bbox_sf))),
    admin     = sf::st_make_valid(sf::st_crop(sf::st_make_valid(admin), sf::st_bbox(bbox_sf)))
  )
}
base <- get_basemap()

# ---- 3. Climate rasters ------------------------------------------------------
CACHE_TIF <- file.path(DIR_CACHE, sprintf("worldclim_bio1_bio12_%smin.tif", WC_RES))

#' Return a 2-layer SpatRaster (precip = bio12 mm, temp = bio1 degC) cropped to
#' STUDY_EXTENT, plus a flag saying whether it is the real WorldClim product.
get_climate_panels <- function() {

  ## tier 1 -- cached crop
  if (file.exists(CACHE_TIF)) {
    r <- terra::rast(CACHE_TIF)
    names(r) <- c("precip", "temp")
    message("climate: cached WorldClim crop")
    return(list(r = r, real = TRUE))
  }

  ## tier 2 -- download WorldClim 2.1 via geodata
  if (requireNamespace("geodata", quietly = TRUE)) {
    message("climate: fetching WorldClim 2.1 (", WC_RES,
            " arc-min). First run only; ~130 MB, a few minutes.")
    ok <- try({
      wc <- geodata::worldclim_global(var = "bio", res = WC_RES, path = DIR_CACHE)
      idx <- function(b) grep(sprintf("bio_%d$", b), names(wc))
      r <- c(wc[[idx(12)]], wc[[idx(1)]])
      r <- terra::crop(r, terra::ext(STUDY_EXTENT[c("xmin","xmax","ymin","ymax")]))
      names(r) <- c("precip", "temp")
      terra::writeRaster(r, CACHE_TIF, overwrite = TRUE,
                         gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
      r
    }, silent = TRUE)
    if (!inherits(ok, "try-error")) {
      message("climate: downloaded WorldClim 2.1 (", WC_RES, " arc-min) -> cached")
      return(list(r = ok, real = TRUE))
    }
    warning("WorldClim download failed: ",
            conditionMessage(attr(ok, "condition")),
            "\n  A partial download leaves a broken archive behind. Clear it with:",
            "\n    unlink(file.path(DIR_CACHE, \"climate\"), recursive = TRUE)",
            "\n  then re-run. See 01_code/00b_check_environment.R.",
            call. = FALSE)
  } else {
    warning("package 'geodata' not installed -- install.packages('geodata')")
  }

  ## tier 3 -- offline stand-in, interpolated from the accession bioclim values
  warning("Falling back to an INTERPOLATED climate surface. ",
          "Draft only -- do not submit this version.", call. = FALSE)
  cl <- read_csv(file.path(DIR_RAW, "climatic_data_allaccessions.csv"),
                 show_col_types = FALSE) |>
    filter(!is.na(longitude), !is.na(latitude)) |>
    select(x = longitude, y = latitude,
           precip = Annual_Precipitation, temp = Annual_Mean_Temperature)

  grid <- terra::rast(terra::ext(STUDY_EXTENT[c("xmin","xmax","ymin","ymax")]),
                      resolution = 1 / 12, crs = "EPSG:4326")
  xy   <- terra::xyFromCell(grid, seq_len(terra::ncell(grid)))
  idw  <- function(v, p = 2.5, k = 12) {
    out <- numeric(nrow(xy))
    for (i in seq_len(nrow(xy))) {
      d <- sqrt((cl$x - xy[i, 1])^2 + (cl$y - xy[i, 2])^2) + 1e-6
      o <- order(d)[seq_len(min(k, length(d)))]
      w <- 1 / d[o]^p
      out[i] <- sum(w * v[o]) / sum(w)
    }
    out
  }
  r <- c(terra::setValues(grid, idw(cl$precip)),
         terra::setValues(grid, idw(cl$temp)))
  names(r) <- c("precip", "temp")
  r <- terra::mask(r, terra::vect(sf::st_union(base$countries)))
  list(r = r, real = FALSE)
}

clim_panels <- get_climate_panels()
clim_r      <- clim_panels$r
IS_REAL_WC  <- clim_panels$real

ras_df <- terra::as.data.frame(clim_r, xy = TRUE, na.rm = TRUE) |>
  as_tibble()

# ---- 4. Map furniture --------------------------------------------------------
# Scale bar + north arrow on a translucent plate in the lower-right corner, so
# they stay legible whatever the raster does underneath.
map_furniture <- function(km = 400, x_right = -94.6, y_bottom = 25.9) {
  deg  <- km / (111.32 * cos((y_bottom + 0.6) * pi / 180))   # bar length in degrees
  bx1  <- x_right - 0.35
  bx0  <- bx1 - deg
  ax   <- bx0 - 1.25                                          # arrow centre
  list(
    annotate("rect", xmin = ax - 1.05, xmax = x_right, ymin = y_bottom - 0.35,
             ymax = y_bottom + 2.35, fill = alpha("white", 0.86), colour = NA),
    # north arrow
    annotate("polygon",
             x = c(ax, ax - 0.42, ax, ax + 0.42),
             y = y_bottom + c(2.00, 0.62, 0.82, 0.62),
             fill = INK[["primary"]], colour = INK[["primary"]], linewidth = 0.2),
    annotate("text", x = ax, y = y_bottom + 0.22, label = "N",
             size = BASE_SIZE / .pt * 0.8, fontface = "bold",
             colour = INK[["primary"]]),
    # scale bar
    annotate("rect", xmin = bx0, xmax = bx1, ymin = y_bottom + 0.35,
             ymax = y_bottom + 0.72,
             fill = "white", colour = INK[["primary"]], linewidth = 0.25),
    annotate("rect", xmin = bx0, xmax = bx0 + deg / 2, ymin = y_bottom + 0.35,
             ymax = y_bottom + 0.72,
             fill = INK[["primary"]], colour = INK[["primary"]], linewidth = 0.25),
    annotate("text", x = (bx0 + bx1) / 2, y = y_bottom + 1.35,
             label = paste0(km, " km"),
             size = BASE_SIZE / .pt * 0.75, colour = INK[["primary"]])
  )
}

# ---- 5. Panel builder --------------------------------------------------------
build_panel <- function(var, ramp, legend_title, breaks, labels = waiver(),
                        show_y = TRUE, furniture = FALSE) {

  d <- ras_df |> select(x, y, value = all_of(var))

  ggplot() +
    # land outside the raster footprint, so coastlines stay legible
    geom_sf(data = base$countries, fill = LAND_FILL, colour = NA) +
    geom_raster(data = d, aes(x, y, fill = value), interpolate = TRUE) +
    scale_fill_gradientn(
      colours = ramp, name = legend_title,
      breaks = breaks, labels = labels,
      limits = range(breaks), oob = scales::squish,
      guide = guide_bar(legend_title, width = 66)
    ) +
    # administrative boundaries drawn over the raster in white / grey
    geom_sf(data = base$admin, fill = NA,
            colour = alpha("white", 0.75), linewidth = 0.16) +
    geom_sf(data = base$countries, fill = NA,
            colour = INK[["secondary"]], linewidth = 0.30) +
    # accessions
    # white halo, then the coloured symbol: reads as a filled point with a
    # white outline while leaving `fill` to the raster scale
    geom_point(data = acc, aes(lon_j, lat_j, shape = species,
                               size = I(POINT_SIZE * size_adj + 0.55)),
               colour = "white") +
    geom_point(data = acc, aes(lon_j, lat_j, colour = species, shape = species,
                               size = I(POINT_SIZE * size_adj)),
               alpha = 0.95) +
    scale_colour_manual(values = SPECIES_COLS, guide = "none") +
    scale_shape_manual(values = SPECIES_SHAPES, guide = "none") +
    { if (furniture) map_furniture() } +
    coord_sf(xlim = STUDY_EXTENT[c("xmin", "xmax")],
             ylim = STUDY_EXTENT[c("ymin", "ymax")],
             expand = FALSE, crs = 4326) +
    scale_x_continuous(breaks = seq(-120, -95, by = 5), labels = lab_lon) +
    scale_y_continuous(breaks = seq(25, 40, by = 5), labels = lab_lat) +
    theme_map() +
    legend_inside(0.015, 0.02) +
    theme(
      legend.direction     = "horizontal",
      legend.background    = element_rect(fill = alpha("white", 0.88), colour = NA),
      legend.margin        = margin(3, 5, 3, 5),
      legend.title         = element_text(size = rel(0.8), face = "plain",
                                          colour = INK[["primary"]]),
      legend.text          = element_text(size = rel(0.72), colour = INK[["secondary"]]),
      axis.text.y          = if (show_y) element_text(size = rel(0.8),
                                                      colour = INK[["primary"]])
                             else element_blank(),
      axis.ticks.y         = if (show_y) element_line(colour = INK[["secondary"]],
                                                      linewidth = 0.25)
                             else element_blank()
    )
}

pA <- build_panel(
  "precip", RAMP_PRECIP, "Annual precipitation (mm)",
  breaks = c(0, 300, 600, 900, 1200), furniture = TRUE
) +
  # The letter is a title, not a tag: a tag is drawn inside the plot area and a
  # map fills its panel edge to edge, so the letter would sit on the raster.
  labs(title = "A") +
  theme(plot.title = element_text(size = rel(1.7), face = "bold", hjust = 0,
                                  margin = margin(b = 2)),
        plot.title.position = "plot")

pB <- build_panel(
  "temp", RAMP_TEMP, paste0("Annual mean temperature (", DEG, "C)"),
  breaks = c(-5, 0, 5, 10, 15, 20, 25), show_y = FALSE
) +
  labs(title = "B") +
  theme(plot.title = element_text(size = rel(1.7), face = "bold", hjust = 0,
                                  margin = margin(b = 2)),
        plot.title.position = "plot")

# ---- 6. Species legend panel -------------------------------------------------
# Built by hand so species can be grouped under their geographic-origin heading
# and annotated with the number of mapped accessions.
legend_panel <- function(lut, ncol_pad = 0.06) {
  groups <- levels(droplevels(lut$group))
  d <- lut |>
    mutate(size_adj = unname(SPECIES_SIZEADJ[as.character(species)])) |>
    group_by(group) |>
    mutate(row = dplyr::row_number()) |>
    ungroup() |>
    mutate(col = match(as.character(group), groups),
           x   = col, y = -row)
  heads <- data.frame(group = groups, x = seq_along(groups), y = 0)
  n_row <- max(d$row)

  ggplot() +
    geom_text(data = heads, aes(x = x - 0.02, y = y, label = toupper(group)),
              hjust = 0, vjust = 0.5, size = BASE_SIZE / .pt * 0.70,
              fontface = "bold", colour = INK[["secondary"]]) +
    geom_segment(data = heads,
                 aes(x = x - 0.02, xend = x + 0.93, y = y - 0.42, yend = y - 0.42),
                 colour = INK[["faint"]], linewidth = 0.3) +
    geom_text(data = heads, aes(x = x + 0.93, y = y, label = "n"),
              hjust = 1, vjust = 0.5, size = BASE_SIZE / .pt * 0.70,
              fontface = "bold", colour = INK[["secondary"]]) +
    geom_point(data = d, aes(x = x + 0.03, y = y, shape = species,
                             size = I((POINT_SIZE + 0.35) * size_adj + 0.85)),
               colour = "white") +
    geom_point(data = d, aes(x = x + 0.03, y = y, colour = species, shape = species,
                             size = I((POINT_SIZE + 0.35) * size_adj))) +
    geom_text(data = d, aes(x = x + 0.10, y = y, label = label),
              hjust = 0, vjust = 0.5, size = BASE_SIZE / .pt * 0.85,
              fontface = "italic", colour = INK[["primary"]]) +
    geom_text(data = d, aes(x = x + 0.93, y = y, label = n_mapped),
              hjust = 1, vjust = 0.5, size = BASE_SIZE / .pt * 0.78,
              colour = INK[["muted"]]) +
    scale_colour_manual(values = SPECIES_COLS, guide = "none") +
    scale_shape_manual(values = SPECIES_SHAPES, guide = "none") +
    coord_cartesian(xlim = c(1 - ncol_pad, length(groups) + 1 - ncol_pad),
                    ylim = c(-n_row - 0.55, 0.55), expand = FALSE) +
    theme_void(base_size = BASE_SIZE) +
    theme(plot.margin = margin(2, 2, 2, 2))
}

pL <- legend_panel(sp_on_map)

# ---- 7. Assemble -------------------------------------------------------------
# NB: no `&` here. patchwork's `&` operator dispatches through S3 and fails
# against ggplot2 >= 4.0 themes (which are S7): "Can't find method for generic
# `&(e1, e2)`". Panel-level theming is applied to pA / pB directly instead.
# ---- 7b. Climate-of-origin association panels --------------------------------
# Written by 14_climate_rooting.R. If that has not been run the figure falls
# back to the two maps, so this script never depends on it.
assoc_file <- file.path(DIR_CLEAN, "climate_rooting_associations.csv")
if (!file.exists(assoc_file)) {
  cr_path <- file.path(DIR_CODE, "14_climate_rooting.R")
  if (file.exists(cr_path)) source(cr_path, local = new.env())
}
HAVE_ASSOC <- file.exists(assoc_file)

if (HAVE_ASSOC) {
  suppressPackageStartupMessages(library(tidyr))
  assoc <- readr::read_csv(assoc_file, show_col_types = FALSE)

  # Tidy the bioclim names for an axis: "Precipitation_Warmest_Quarter" is a
  # column header, not a label.
  pretty_clim <- function(x) {
    x |> str_replace_all("_", " ") |> str_replace_all("\\.", " ") |>
      str_replace("Mean Temperature", "Mean temp.") |>
      str_replace("Min Temperature", "Min temp.") |>
      str_replace("Max Temperature", "Max temp.") |>
      str_replace("Annual Mean Temperature", "Annual mean temp.") |>
      str_replace("Precipitation", "Precip.") |>
      str_replace("Mean Diurnal Range", "Diurnal range") |>
      str_replace("^altitude$", "Altitude") |>
      str_squish() |>
      sentence_case()
  }
  # Sentence case, so panel C's tick labels and panel D's axis title agree.
  # The bioclim column headers are Title_Case ("Precipitation_Driest_Month"),
  # which came through as "Precip. Driest Month" and sat next to a y-axis title
  # written "Root weight BLUP" -- two capitalisation styles in one panel.
  # Only the first word keeps its capital; all-caps tokens (BLUP, SD) are left
  # alone in case one ever reaches this function.
  sentence_case <- function(x) {
    vapply(x, function(s) {
      w <- strsplit(s, " ", fixed = TRUE)[[1]]
      if (length(w) < 2) return(s)
      i <- seq_along(w)[-1]
      acro <- grepl("^[A-Z]{2,}$", w)
      w[i][!acro[i]] <- tolower(w[i][!acro[i]])
      paste(w, collapse = " ")
    }, character(1), USE.NAMES = FALSE)
  }
  TRAIT_ORDER <- c("Callus-Stage", "Post-transplant", "Root weight",
                   "Cutting weight", "Transplant response")
  TRAIT_NICE  <- c("Callus-stage\nrooting", "Post-transplant\nrooting", "Root\nweight",
                   "Cutting\nweight", "Transplant response\nindex")

  aa <- assoc |>
    mutate(trait = factor(trait, levels = TRAIT_ORDER,
                          labels = TRAIT_NICE),
           clim_lab = pretty_clim(climate),
           star = case_when(q_within < 0.01 ~ "**",
                            q_within < 0.05 ~ "*",
                            q_within < 0.10 ~ "\u00b7",
                            TRUE ~ "")) |>
    filter(!is.na(trait))
  # order climate variables by how strongly they associate overall, so the
  # panel reads left to right instead of in file order
  clim_ord <- aa |> group_by(clim_lab) |>
    summarise(m = max(abs(beta_within)), .groups = "drop") |>
    arrange(desc(m)) |> pull(clim_lab)
  aa <- aa |> mutate(clim_lab = factor(clim_lab, levels = clim_ord))

  # Reviewer 1 asked for the n column to be populated. The number of
  # georeferenced accessions contributing to a row is constant within a trait,
  # so it is carried in the axis label rather than as a separate column.
  n_by_trait <- aa |>
    group_by(trait) |>
    summarise(n = max(n, na.rm = TRUE), .groups = "drop")
  TRAIT_LAB <- setNames(
    sprintf("%s (n = %d)", as.character(n_by_trait$trait), n_by_trait$n),
    as.character(n_by_trait$trait))

  lim_c <- max(abs(aa$beta_within), na.rm = TRUE)
  pC <- ggplot(aa, aes(clim_lab, trait, fill = beta_within)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = star), size = 2.5, vjust = 0.78,
              colour = INK[["primary"]]) +
    scale_fill_gradientn(colours = RAMP_DIV, limits = c(-lim_c, lim_c),
                         name = "Within-species slope\n(change in BLUP per SD)",
                         guide = guide_colourbar(title.position = "top")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), limits = rev(TRAIT_NICE),
                     labels = TRAIT_LAB) +
    labs(x = NULL, y = NULL, title = "C") +
    theme_paper() +
    theme(plot.title = element_text(face = "bold", size = rel(1.7), hjust = 0,
                                    margin = margin(b = 2)),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
          axis.text.y = element_text(size = rel(0.75), lineheight = 1.05),
          panel.grid = element_blank(),
          legend.position = "bottom",
          legend.key.width = unit(30, "pt"), legend.key.height = unit(4, "pt"),
          legend.title = element_text(size = rel(0.75)),
          legend.margin = margin(t = -2),
          # top margin so the panel tag sits ABOVE the first row label instead
          # of on top of it
          # The tag is pinned near the top-left of the whole plot, which is
          # exactly where the first row's axis label sits. Reserving a top strip
          # and putting the tag in it is the only placement that reliably keeps
          # the two apart.
          # 45-degree tick labels overhang to the left of the panel; without
          # this the first one ("Precip. Driest Month") is clipped.
          plot.margin = margin(t = 4, r = 4, b = 4, l = 26))

  # Panel D: the single strongest association, drawn as data rather than asserted
  best <- assoc |> slice_min(q_within, n = 1)
  bx <- best$climate[1]; bt <- best$trait[1]
  pan <- readr::read_csv(file.path(DIR_CLEAN, "accessions_panel.csv"),
                         show_col_types = FALSE)
  bl <- readr::read_csv(file.path(DIR_CLEAN, "blups_by_trait.csv"),
                        show_col_types = FALSE)
  dd <- pan |> filter(georeferenced) |>
    select(genotype, species, all_of(bx)) |>
    inner_join(bl |> select(genotype, value = all_of(paste0(bt, "_blup"))),
               by = "genotype") |>
    filter(is.finite(.data[[bx]]), is.finite(value)) |>
    mutate(species = factor(species, levels = SPECIES_LEVELS),
           species = droplevels(species))

  pD <- ggplot(dd, aes(.data[[bx]], value)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = INK[["primary"]],
                fill = INK[["faint"]], linewidth = 0.5) +
    geom_point(aes(colour = species, shape = species),
               size = 1.15, stroke = 0.3, alpha = 0.9) +
    scale_colour_manual(values = SPECIES_COLS[levels(dd$species)], guide = "none") +
    scale_shape_manual(values = SPECIES_SHAPES[levels(dd$species)], guide = "none") +
    labs(x = paste0(pretty_clim(bx), " (mm)"),
         y = paste0(bt, " BLUP")) +
    theme_paper()
}

# ---- 7. Assemble -------------------------------------------------------------
# NB: no `&` here. patchwork's `&` operator dispatches through S3 and fails
# against ggplot2 >= 4.0 themes (which are S7): "Can't find method for generic
# `&(e1, e2)`". Panel-level theming is applied to pA / pB directly instead.
if (HAVE_ASSOC) {
  # wrap_elements isolates the bottom row: without it patchwork aligns D's x-axis
  # title to the bottom of C's rotated tick labels and leaves a large gap.
  pD <- pD + labs(title = "D") +
    theme(plot.title = element_text(face = "bold", size = rel(1.7), hjust = 0,
                                    margin = margin(b = 2)),
          plot.title.position = "plot")
  row_cd <- wrap_elements(full = pC) | wrap_elements(full = pD)
  fig1 <- (pA | pB) / pL / row_cd +
    plot_layout(heights = c(1, 0.34, 0.95))
  FIG_H <- 6.5
} else {
  fig1 <- (pA | pB) / pL +
    plot_layout(heights = c(1, 0.36))
  FIG_H <- 3.35
}

if (!IS_REAL_WC) {
  fig1 <- fig1 +
    plot_annotation(caption = paste(
      "DRAFT: interpolated climate surface.",
      "install.packages('geodata') and re-run for the WorldClim 2.1 background."),
      theme = theme(plot.caption = element_text(colour = "#B8453C", face = "bold",
                                                size = rel(0.9), hjust = 0)))
}

# Two-column journal figure: 180 mm wide.
save_figure(fig1, "Figure3", width = 7.09, height = FIG_H)
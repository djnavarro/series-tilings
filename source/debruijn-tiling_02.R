# set up ------------------------------------------------------------------
#
# R translation of the de Bruijn grid method for generating aperiodic rhombus
# tilings, ported from https://github.com/samacqua/tilings (Python/matplotlib)
# into a de Bruijn-inspired R/ggplot2 version.
#
# The construction follows de Bruijn (1981): N families of parallel lines are
# drawn at angles j*pi/N, each family offset by a small "phase" sigma_j. Here
# the phases are drawn at random for every piece, which is the standard way to
# get a *random* instance of a de Bruijn tiling (for N = 5 this gives a random
# Penrose tiling; other N give generalized aperiodic tilings). Rhombi sit at
# the pairwise intersections of the grid lines and are then glued edge-to-edge
# via breadth-first search to assemble the final gap-free tiling.

name <- "debruijn-tiling" # only change this for an entirely new system
version <- 2              # increment to 2 for "debruijn-tiling_02.R", etc
format <- "png"           # png is usually a good choice!

library(ggplot2)

# define common helper functions
source(here::here("source", "common.R"), echo = FALSE)

# make sure we haven't accidentally messed up the versions
assert_version_consistency(name)


# grid construction --------------------------------------------------------

# Random phase offsets in (0, 1), nudged apart so that (a) no two offsets
# coincide and (b) no two offsets sum to an integer -- both are required so
# that no three grid lines ever meet at a single point.
sample_offsets <- function(n, eps = 1e-6) {
  offsets <- runif(n)
  repeat {
    bad <- outer(offsets, offsets, function(a, b) abs(a - b) < eps) |
      outer(offsets, offsets, function(a, b) abs(a + b - 1) < eps)
    diag(bad) <- FALSE
    if (!any(bad)) break
    idx <- which(rowSums(bad) > 0)
    offsets[idx] <- runif(length(idx))
  }
  offsets
}

# N families of parallel lines filling the bounding box. Returns one row per
# clipped line segment; each family shares a common angle, and lines within a
# family are indexed by an integer k.
generate_grid_lines <- function(xlim, ylim, offsets, n_line_families, ang_offset = 0) {
  rows <- vector("list", n_line_families)

  for (j in seq_len(n_line_families)) {
    angle <- (j - 1) * pi / n_line_families + ang_offset
    cos_a <- cos(angle)
    sin_a <- sin(angle)
    offset <- offsets[j]

    corners_x <- c(xlim[1], xlim[2], xlim[1], xlim[2])
    corners_y <- c(ylim[1], ylim[1], ylim[2], ylim[2])
    k_vals <- corners_x * cos_a + corners_y * sin_a + offset
    k_min <- floor(min(k_vals))
    k_max <- ceiling(max(k_vals))
    ks <- (k_min - 1):(k_max + 1)

    fam_rows <- lapply(ks, function(k) {
      pts <- matrix(nrow = 0, ncol = 2)

      if (abs(sin_a) > EPS) {
        for (y in ylim) {
          x <- (k - offset - y * sin_a) / cos_a
          if (xlim[1] <= x && x <= xlim[2]) pts <- rbind(pts, c(x, y))
        }
      }
      if (abs(cos_a) > EPS) {
        for (x in xlim) {
          y <- (k - offset - x * cos_a) / sin_a
          if (ylim[1] <= y && y <= ylim[2]) pts <- rbind(pts, c(x, y))
        }
      }
      if (nrow(pts) < 2) return(NULL)
      pts <- unique(pts)[1:2, , drop = FALSE]
      data.frame(
        line_family = j, k = k, angle = angle,
        x1 = pts[1, 1], y1 = pts[1, 2], x2 = pts[2, 1], y2 = pts[2, 2]
      )
    })
    rows[[j]] <- do.call(rbind, fam_rows)
  }

  do.call(rbind, rows)
}

# Pairwise intersections between *different* line families, vectorised over
# every (k1, k2) combination within each family pair.
find_intersections <- function(lines, xlim, ylim, offsets) {
  families <- sort(unique(lines$line_family))
  results <- list()

  for (pair in combn(families, 2, simplify = FALSE)) {
    f1 <- pair[1]; f2 <- pair[2]
    l1 <- lines[lines$line_family == f1, ]
    l2 <- lines[lines$line_family == f2, ]
    a1 <- l1$angle[1]; a2 <- l2$angle[1]
    cos1 <- cos(a1); sin1 <- sin(a1)
    cos2 <- cos(a2); sin2 <- sin(a2)
    det <- cos1 * sin2 - cos2 * sin1
    if (abs(det) < EPS) next

    k1_adj <- l1$k - offsets[f1]
    k2_adj <- l2$k - offsets[f2]
    grid <- expand.grid(i1 = seq_along(k1_adj), i2 = seq_along(k2_adj))
    k1 <- k1_adj[grid$i1]; k2 <- k2_adj[grid$i2]

    x <- (k1 * sin2 - k2 * sin1) / det
    y <- (k2 * cos1 - k1 * cos2) / det
    keep <- xlim[1] <= x & x <= xlim[2] & ylim[1] <= y & y <= ylim[2]
    if (!any(keep)) next

    results[[length(results) + 1]] <- data.frame(
      x = x[keep], y = y[keep],
      f1 = f1, f2 = f2, a1 = a1, a2 = a2,
      k1 = l1$k[grid$i1[keep]], k2 = l2$k[grid$i2[keep]]
    )
  }

  do.call(rbind, results)
}


# rhombus geometry ----------------------------------------------------------

# Vertices of the rhombus sitting at `center`, with edges perpendicular to the
# two crossing grid-line angles, plus its smallest internal angle (radians) --
# used later to decide the tile's colour.
create_rhombus <- function(center, angle1, angle2, scale = 0.5) {
  angle1 <- angle1 + pi / 2
  angle2 <- angle2 + pi / 2
  d1 <- c(-sin(angle1), cos(angle1))
  d2 <- c(-sin(angle2), cos(angle2))

  verts <- rbind(
    center + scale * d1 + scale * d2,
    center - scale * d1 + scale * d2,
    center - scale * d1 - scale * d2,
    center + scale * d1 - scale * d2
  )

  diff <- abs((angle2 - angle1) %% (2 * pi))
  if (diff > pi) diff <- 2 * pi - diff
  smallest_angle <- if (diff > pi / 2) pi - diff else diff

  list(vertices = verts, smallest_angle = smallest_angle)
}

# Vertices with the min/max projection onto `direction`.
extreme_vertices <- function(vertices, direction, farthest = TRUE) {
  proj <- as.numeric(vertices %*% direction)
  extreme <- if (farthest) max(proj) else min(proj)
  vertices[abs(proj - extreme) < 1e-8, , drop = FALSE]
}


# assembling rhombi into a data structure ------------------------------------

# Two rhombi are neighbours if they sit next to each other along one of their
# two shared grid lines (i.e. share the same (family, k) pair).
compute_neighbors <- function(rhombi) {
  n <- length(rhombi)
  groups <- list()
  for (idx in seq_len(n)) {
    r <- rhombi[[idx]]
    for (side in 1:2) {
      key <- paste(r$fam[side], r$k[side], sep = "_")
      groups[[key]] <- c(groups[[key]], idx)
    }
  }

  neighbours <- vector("list", n)
  for (key in names(groups)) {
    grp <- groups[[key]]
    fam_k <- strsplit(key, "_")[[1]]
    fam <- as.integer(fam_k[1])
    angle <- rhombi[[grp[1]]]$angle[rhombi[[grp[1]]]$fam == fam][1]
    d <- c(-sin(angle), cos(angle))
    proj <- sapply(grp, function(idx) sum(rhombi[[idx]]$center * d))
    grp <- grp[order(proj)]
    for (j in seq_along(grp)) {
      if (j > 1) neighbours[[grp[j]]] <- union(neighbours[[grp[j]]], grp[j - 1])
      if (j < length(grp)) neighbours[[grp[j]]] <- union(neighbours[[grp[j]]], grp[j + 1])
    }
  }
  neighbours
}

# Glue every rhombus edge-to-edge via breadth-first search, starting from an
# arbitrary seed tile, producing the final assembled (gap-free) tiling.
build_tiling <- function(rhombi, neighbours, scale = 1) {
  n <- length(rhombi)
  centred <- lapply(rhombi, function(r) {
    v <- r$vertices
    sweep(v, 2, colMeans(v)) * scale
  })

  place_neighbor <- function(base_idx, neigh_idx, placed) {
    r_base <- rhombi[[base_idx]]; r_neigh <- rhombi[[neigh_idx]]
    shared <- NA
    for (side in 1:2) {
      match_side <- which(r_neigh$fam == r_base$fam[side] & r_neigh$k == r_base$k[side])
      if (length(match_side) > 0) { shared <- r_base$angle[side]; break }
    }
    d <- c(-sin(shared), cos(shared))

    sgn <- sign(sum((r_neigh$center - r_base$center) * d))
    if (sgn == 0) sgn <- 1

    front <- extreme_vertices(placed[[base_idx]], d, farthest = (sgn > 0))
    rear <- extreme_vertices(centred[[neigh_idx]], d, farthest = (sgn < 0))
    shift <- colMeans(front) - colMeans(rear)
    placed[[neigh_idx]] <- sweep(centred[[neigh_idx]], 2, shift, "+")
    placed
  }

  placed <- vector("list", n)
  seed <- 1
  placed[[seed]] <- centred[[seed]]
  queue <- seed
  seen <- rep(FALSE, n); seen[seed] <- TRUE
  head <- 1

  while (head <= length(queue)) {
    idx <- queue[head]; head <- head + 1
    for (nbr in neighbours[[idx]]) {
      if (seen[nbr]) next
      placed <- place_neighbor(idx, nbr, placed)
      seen[nbr] <- TRUE
      queue <- c(queue, nbr)
    }
  }
  placed
}


# top-level pipeline ---------------------------------------------------------

EPS <- 1e-10

# Default tile palette: qualitative Set3, extended by interpolation if more
# than 12 distinct rhombus shapes are present. Any function(n) -> character
# vector of n colours can be substituted via the `palette` argument below.
default_debruijn_palette <- function(n) {
  if (n <= 12) {
    RColorBrewer::brewer.pal(max(3, n), "Set3")[seq_len(n)]
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
  }
}

# Build the full data structure for one random de Bruijn tiling: draws random
# phase offsets, computes the grid-line intersections, turns each one into a
# rhombus, and glues all rhombi into a single assembled tiling. Returns a
# tidy data frame (one row per tile vertex) ready for ggplot2.
generate_debruijn_tiling <- function(n_line_families = 5, x_scale = 5, y_scale = 5,
                                      rhombus_scale = 1, palette = default_debruijn_palette) {
  xlim <- c(-x_scale, x_scale)
  ylim <- c(-y_scale, y_scale)
  offsets <- sample_offsets(n_line_families)

  lines <- generate_grid_lines(xlim, ylim, offsets, n_line_families)
  intersections <- find_intersections(lines, xlim, ylim, offsets)

  rhombi <- lapply(seq_len(nrow(intersections)), function(i) {
    row <- intersections[i, ]
    rh <- create_rhombus(c(row$x, row$y), row$a1, row$a2, scale = rhombus_scale / 2)
    list(
      vertices = rh$vertices,
      smallest_angle = rh$smallest_angle,
      center = c(row$x, row$y),
      fam = c(row$f1, row$f2),
      k = c(row$k1, row$k2),
      angle = c(row$a1, row$a2)
    )
  })

  # colour each tile by its (rounded) smallest internal angle, so each
  # distinct rhombus shape gets its own colour -- mirrors the Python original
  bucketed <- round(sapply(rhombi, function(r) r$smallest_angle), 4)
  unique_angles <- sort(unique(bucketed))
  tile_palette <- palette(length(unique_angles))
  tile_colour <- tile_palette[match(bucketed, unique_angles)]

  neighbours <- compute_neighbors(rhombi)
  placed <- build_tiling(rhombi, neighbours, scale = rhombus_scale)

  # a handful of rhombi right at the box edge can end up with no shared-line
  # neighbours (especially with a very asymmetric x_scale/y_scale) and are
  # never reached by the BFS; drop those rather than erroring
  reached <- !vapply(placed, is.null, logical(1))
  if (!all(reached)) {
    message(sum(!reached), " isolated edge tile(s) dropped from the assembled tiling")
  }

  do.call(rbind, lapply(which(reached), function(idx) {
    v <- placed[[idx]]
    data.frame(
      tile_id = idx,
      vertex = 1:4,
      x = v[, 1], y = v[, 2],
      colour = tile_colour[idx]
    )
  }))
}

# Render an assembled tiling with each rhombus shaded according to its shape
# (smallest internal angle), styled after the Python project's renders:
# transparent background, no axes, equal aspect ratio.
plot_debruijn_tiling <- function(tiling) {
  ggplot(tiling, aes(x = x, y = y, group = tile_id, fill = colour)) +
    geom_polygon(colour = "black", linewidth = 0.15) +
    scale_fill_identity() +
    coord_equal() +
    theme_void()
}

# Probabilistic palette functions via scico
sample_palette <- function(seed) {
  known_palettes <- scico::scico_palette_names()
  palette_name <- sample(known_palettes, 1)
  palette_fn <- function(n) {
    scico::scico(n = n, palette = palette_name)
  }
  return(palette_fn)
}


# make the art ---------------------------------------------------------------

art_generator <- function(seed) {

  set.seed(seed)

  # specify the output path and message the user
  output <- output_path(name, version, seed, format)
  message("generating art at ", output)

  # generate high-level parameters for this piece
  n_line_families <- sample(5:9, size = 1)
  y_scale <- 15
  x_scale <- 15
  palette <- sample_palette(seed)

  # generate the tiling data structure and render it
  tiling <- generate_debruijn_tiling(
    n_line_families = n_line_families, 
    x_scale = x_scale, 
    y_scale = y_scale,
    palette = palette
  )
  p <- plot_debruijn_tiling(tiling)

  ggsave(output, plot = p, width = 10, height = 10, units = "in", dpi = 300, bg = "#ffffff00")
}

for (seed in 101:200) art_generator(seed)

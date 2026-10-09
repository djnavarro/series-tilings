# Setup
name    <- "debruijn-tiling"
version <- 10
format  <- "png"
source(here::here("source", "common.R"), echo = FALSE)
assert_version_consistency(name)

# Packages
library(ggplot2)
library(tibble)
library(purrr)
library(tidyr)
library(dplyr)

# Every intersection of two grid lines becomes a rhombus, in the usual
# "dual" construction for these grid methods: the rhombus edges run
# perpendicular to the two line directions that cross at this point, and
# its two pairs of opposite edges have length `2 * scale`.
build_rhombus <- function(x_mid,
                          y_mid,
                          angle_1,
                          angle_2,
                          offset_1 = NA_real_,
                          offset_2 = NA_real_,
                          family_1 = NA_real_,
                          family_2 = NA_real_,
                          line_1 = NA_real_,
                          line_2 = NA_real_,
                          scale = 1) {

  # Rhombus edge directions are rotated 90 degrees from the line angles; 
  # this ensures that rhombi that sit on the same line will have matching
  # faces, and can be tiled next to one another without gaps when the
  # tiling is assembled
  rh_angle_1 <- angle_1 + pi / 2
  rh_angle_2 <- angle_2 + pi / 2

  # These are unit vectors pointing along the two pairs of edges
  d1 <- c(-sin(rh_angle_1), cos(rh_angle_1))
  d2 <- c(-sin(rh_angle_2), cos(rh_angle_2))

  # Centroid coordinates as a named vector; names become the column names
  # for the vertices matrix, and later become column names for the data
  # frame
  centre <- c(x = x_mid, y = y_mid)

  # The four corners, reached from the centre by moving +/- scale along
  # each of the two edge directions in turn; it is convenient to represent
  # vertices as a matrix rather than a data frame, even though plotting is
  # easier in the latter format. 
  vertices <- rbind(
    centre + scale * d1 + scale * d2,
    centre - scale * d1 + scale * d2,
    centre - scale * d1 - scale * d2,
    centre + scale * d1 - scale * d2
  )

  # The angle *between* the two grid-line directions determines the shape
  # of the rhombus (a 72-degree gap between line families gives "thin"
  # Penrose rhombi, a 36-degree gap gives "fat" ones, etc). Reducing this
  # angle to its smallest equivalent (the acute angle, <= pi/2) gives a
  # single number that identifies the tile's shape regardless of its
  # orientation: it's not strictly needed to construct the tiling, but gets 
  # used later to decide the fill colour of a tile
  diff <- abs((rh_angle_2 - rh_angle_1) %% (2 * pi))
  if (diff > pi) diff <- 2 * pi - diff
  tile_acute <- if (diff > pi / 2) pi - diff else diff

  # Return the built rhombus as a list containing a matrix of vertices,
  # plus a one-row tibble of metadata: the acute angle (for colouring)
  # and, crucially, *which* two grid lines (family + line index) this
  # tile sits on. Keeping track of this is critical for the tiling process,
  # since compute_neighbours() and build_tiling() later match tiles up by 
  # looking for shared (family, line) pairs
  list(
    vertices = vertices,
    properties = tibble(
      tile_acute = tile_acute,
      x_mid = x_mid,
      y_mid = y_mid,
      angle_1 = angle_1,
      angle_2 = angle_2,
      offset_1 = offset_1,
      offset_2 = offset_2,
      family_1 = family_1,
      family_2 = family_2,
      line_1 = line_1,
      line_2 = line_2
    )
  )
}

# There's no probabilistic component to setting the angles that define
# the line families; only to the horizontal offsets that displace them.
# The default value to `angle_offset` is deliberately slightly non-zero,
# which means that angles are not *strictly* defined relative to north,
# but are ever so slightly rotated. This is purely a numerical convenience,
# to avoid issues that pop up with vertical lines (infinite slope)
build_families <- function(n_families, 
                           offsets = NULL, 
                           angle_offset = 0.01,
                           seed = NULL) {

  # if the line family offsets aren't set, generate them randomly
  if (is.null(offsets)) offsets <- sample_offsets(n_families, seed = seed)

  # return a tibble with summary information for each line family
  tibble(
    family  = seq_len(n_families),
    offset = offsets,
    angle  = (family - 1) * pi / n_families + angle_offset
  ) 
}

# Each of the line families is offset horizontally from the origin by a phase
# These phases are what make the resulting tiling look "random"
sample_offsets <- function(n_families, eps = 1e-6, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  offsets <- runif(n_families)
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

# Handy little infix operator to check if a scalar value falls within the 
# range specified by a length-2 vector
`%within%` <- function(value, range) range[1] <= value & value <= range[2]

# Build the families of parallel lines, and populate each family with a 
# sufficient number of lines to fill the bounding box (defined by the 
# x_scale and y_scale parameters). Returns a data frame with one row for
# each line that passes through the bounding box, with coordinates and 
# other useful information.
build_grid_lines <- function(x_scale, 
                             y_scale, 
                             offsets, 
                             eps = 1e-6) {

  # Set up the line families
  n_families <- length(offsets)
  families <- build_families(n_families, offsets)
    
  # Define the bounding box; keep this internal because it's not really a 
  # plot limit; it behaves like a plot limit for the purposes of this 
  # function, but when the final tiling gets assembled the rhombi get moved 
  # around, so the final plot limits are different
  xlim <- c(-x_scale, x_scale)
  ylim <- c(-y_scale, y_scale)

  # For each family, work out which integer `line` values could possibly
  # intersect the box at all, so later steps don't waste time clipping
  # lines that are nowhere near it. `corner_k` evaluates the line equation
  # `x*cos(angle) + y*sin(angle) + offset` at each of the box's 4 corners.
  # The resulting range of `line` must cover every integer between the
  # smallest and largest corner value (with a 1-unit pad on each side, to
  # be safe about rounding) for the family's lines to span the box
  line_cases <- families |> 
    cross_join(expand_grid(xlim, ylim)) |> 
    mutate(corner_k = xlim * cos(angle) + ylim * sin(angle) + offset) |> 
    reframe(
      .by = family,
      across(c(offset, angle), first),
      line = seq(
        from = floor(min(corner_k)) - 1,
        to = ceiling(max(corner_k)) + 1
      )
    )

  # Function to clip one specific line (one `family`/`line` combination) 
  # to the bounding box. A line can leave the box through a horizontal edge 
  # (at y = ylim) or a vertical edge (at x = xlim); this solves for both and
  # keeps only the crossing points that are actually within the box 
  # (`x_ok`/`y_ok`), then de-duplicates in case a line passes exactly
  # through a corner. Near-horizontal/near-vertical lines need the
  # opposite equation (division by a near-zero cos/sin), which is why both
  # cases are computed and only the valid one(s) kept via `eps`.
  build_grid_line <- function(family, offset, angle, line) {

    points <- bind_rows(
      tibble(
        y = ylim, 
        x = (line - offset - y * sin(angle)) / cos(angle),
        x_ok = x %within% xlim,
        y_ok = abs(sin(angle)) > eps
      ),
      tibble(
        x = xlim, 
        y = (line - offset - x * cos(angle)) / sin(angle),
        y_ok = y %within% ylim,
        x_ok = abs(cos(angle)) > eps
      )
    )

    points <- points |> 
      filter(x_ok & y_ok) |> 
      distinct()

    # A line that misses the box entirely contributes no segment. This
    # isn't expected to happen normally, but retained just in case
    if (nrow(points) == 0) {
      return(tibble(
        family = numeric(0),
        line = numeric(0),
        angle = numeric(0),
        offset = numeric(0),
        x1 = numeric(0),
        y1 = numeric(0),
        x2 = numeric(0),
        y2 = numeric(0)
      ))
    }

    tibble(
      family = family, 
      line = line, 
      angle = angle,
      offset = offset,
      x1 = points$x[1], 
      y1 = points$y[1], 
      x2 = points$x[2], 
      y2 = points$y[2]
    )
  }

  # Clip every candidate line in every family, and return a tibble
  line_cases |> pmap_dfr(build_grid_line)
}

# Find every point where a line from one family crosses a line from
# another family, within the bounding box. Lines in the same family 
# are parallel, so the search is performed across pairs of lines that
# belong to different families.
find_intersections <- function(lines, 
                               x_scale = NULL, 
                               y_scale = NULL, 
                               eps = 1e-10) {

  # Detect the box scale from the lines if the user doesn't specify
  if (is.null(x_scale)) x_scale <- max(abs(c(lines$x1, lines$x2)))
  if (is.null(y_scale)) y_scale <- max(abs(c(lines$y1, lines$y2)))

  # As before: keep this internal because it's not really a plot limit
  xlim <- c(-x_scale, x_scale)
  ylim <- c(-y_scale, y_scale)

  # Find all unique pairs of line families: the base R combn() function
  # is more efficient in general but expand_grid() then filter() is 
  # perfectly fine here
  family <- sort(unique(lines$family))
  family_pairs <- expand_grid(family_1 = family, family_2 = family) |> 
    filter(family_2 > family_1)

  # Helper function used to detect all the intersections between two 
  # line families that fall inside the plotting area
  find_family_intersections <- function(family_1, family_2) {

    # Empty tibble to return if there is no intersection
    null_intersection <- tibble(
      x = numeric(0),
      y = numeric(0),
      angle_1 = numeric(0),
      angle_2 = numeric(0),
      offset_1 = numeric(0),
      offset_2 = numeric(0),
      family_1 = numeric(0),
      family_2 = numeric(0),
      line_1 = numeric(0),
      line_2 = numeric(0)
    )

    # Sets of lines associated with each of the two families
    line_set_1 <- lines |> filter(family == family_1)
    line_set_2 <- lines |> filter(family == family_2)

    # The angles specifying the two line families
    angle_1 <- line_set_1$angle[1]
    angle_2 <- line_set_2$angle[1]

    # The offsets associated with the two line families
    offset_1 <- line_set_1$offset[1]
    offset_2 <- line_set_2$offset[1]

    # Take sin and cosine of both angles
    cos_a1 <- cos(angle_1) 
    sin_a1 <- sin(angle_1)
    cos_a2 <- cos(angle_2) 
    sin_a2 <- sin(angle_2)

    # Return an empty tibble if the determinant is too low. When 
    # det = 0 the lines are exactly parallel and no intersection
    # exists. A threshold is applied to avoid numerical issues
    # that arise with nearly-parallel lines
    det <- cos_a1 * sin_a2 - cos_a2 * sin_a1
    if (abs(det) < eps) return(null_intersection)

    # Determine an "adjusted" coordinate for each line, after taking
    # family-specific offset values into account; then compute the x
    # and y values at which the intersection occurs; retaining only
    # those intersections that remain within the plot limits
    grid <- expand_grid(
      line_1 = line_set_1$line,
      line_2 = line_set_2$line
    ) |> 
      mutate(
        line_coord_1 = line_1 - offset_1,
        line_coord_2 = line_2 - offset_2,
        x = (line_coord_1 * sin_a2 - line_coord_2 * sin_a1) / det,
        y = (line_coord_2 * cos_a1 - line_coord_1 * cos_a2) / det,
        keep = (x %within% xlim) & (y %within% ylim)
      ) |> 
      filter(keep)

    # If none of the pairs are retained return an empty tibble
    if (nrow(grid) == 0) return(null_intersection)

    # Return a tibble containing the retained intersections
    grid |> 
      transmute(
        x = x, 
        y = y,
        angle_1 = angle_1,
        angle_2 = angle_2,
        offset_1 = offset_1,
        offset_2 = offset_2,
        family_1 = family_1,
        family_2 = family_2,
        line_1 = line_1,
        line_2 = line_2
      )
  }

  family_pairs |> pmap_dfr(find_family_intersections)    
}

# Find consecutive pairs for every grid line, and returns them as a tidy 
# edge list (idx, neighbour), with each pair listed in both directions so 
# the breadth-first search that comes later in build_tiling() can look up 
# neighbours from either side.
compute_neighbours <- function(rhombi_tbl) {
  props <- rhombi_tbl$properties

  # Pivot from one-row-per-tile to one-row-per-(tile, side): every tile
  # contributes two rows, one for each of the two grid lines it sits on.
  # In this context bind_rows() is simpler than using pivot_longer()
  sides <- bind_rows(
    props |> select(
      idx, family = family_1, line = line_1, 
      angle = angle_1, x_mid, y_mid
    ),
    props |> select(
      idx, family = family_2, line = line_2, 
      angle = angle_2, x_mid, y_mid
    )
  )

  edges <- sides |>
    # Tiles strung along a grid line can be fully ordered by their position 
    # along the line: `proj` takes the centre of the tile and projects it 
    # onto the direction perpendicular to the line, and increases monotonically
    # as you move along the line
    mutate(proj = -sin(angle) * x_mid + cos(angle) * y_mid) |>

    # Sorting tiles within each (family, line) group by the projected position, 
    # ensures that tiles adjacent in the sorted order are exactly the tiles
    # adjacent in physical space along the line
    arrange(family, line, proj) |>

    # Pair each tile with the very next one in sorted order, within its
    # group. These are edge-sharing neighbour pairs.
    mutate(neighbour = lead(idx), .by = c(family, line)) |>
    filter(!is.na(neighbour)) |>
    select(idx, neighbour)

  # Each pair was only recorded in one direction (idx -> neighbour); add
  # the mirror-image rows so the relation is symmetric
  edges_reversed <- edges |> rename(idx = neighbour, neighbour = idx)

  bind_rows(edges, edges_reversed) |> distinct()
}

# Enumerate every maximal clique in the `neighbours` graph: a clique is a
# set of tiles that are *all* mutually edge-adjacent to one another (not
# just chained together). Uses the Bron-Kerbosch algorithm with pivoting
# (https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm). Because
# each tile sits on only two grid lines and can have at most two neighbours
# along each one, the graph is sparse (max degree 4) and most maximal
# cliques will simply be the individual edges (pairs) recorded in
# `neighbours`; cliques of 3+ tiles occur only when several tiles happen to
# be pairwise adjacent, which is possible but not guaranteed by the tiling
# geometry. Returns a data frame with one row per clique/tile pair.
find_cliques <- function(neighbours) {

  # Every tile that appears anywhere in the edge list is a vertex in the
  # graph; tiles with no recorded edges never show up here at all, and are
  # not part of any clique of size > 1
  all_idx <- sort(unique(c(neighbours$idx, neighbours$neighbour)))

  # Adjacency list keyed by tile index, as in build_tiling()
  adjacency <- split(neighbours$neighbour, neighbours$idx)
  neighbours_of <- function(i) {
    nbrs <- adjacency[[as.character(i)]]
    if (is.null(nbrs)) integer(0) else nbrs
  }

  cliques <- list()

  # Classic recursive Bron-Kerbosch with pivoting. `R` is the clique built
  # so far, `P` the candidates that could still extend it, and `X` the
  # vertices already excluded because any clique containing them would have
  # already been reported by an earlier branch
  bron_kerbosch <- function(R, P, X) {

    # A maximal clique is found once there are no candidates left to add,
    # and no excluded vertex that could have extended it either
    if (length(P) == 0 && length(X) == 0) {
      cliques[[length(cliques) + 1]] <<- R
      return(invisible())
    }

    # Pick the pivot from P union X with the most neighbours inside P; only
    # candidates that are *not* neighbours of the pivot need to be branched
    # on, since cliques containing the pivot itself are guaranteed to be
    # found via one of those other branches. This pruning is what keeps the
    # algorithm efficient relative to naive enumeration
    PX <- union(P, X)
    pivot <- PX[which.max(
      vapply(PX, \(v) length(intersect(P, neighbours_of(v))), integer(1))
    )]
    candidates <- setdiff(P, neighbours_of(pivot))

    for (v in candidates) {
      nbrs_v <- neighbours_of(v)
      bron_kerbosch(
        R = c(R, v),
        P = intersect(P, nbrs_v),
        X = intersect(X, nbrs_v)
      )
      P <- setdiff(P, v)
      X <- union(X, v)
    }
  }

  # Populate the cliques list via Bron-Kerbosch algorithm
  bron_kerbosch(R = integer(0), P = all_idx, X = integer(0))

  # Return a long-format tibble with one row for each tile/clique pair
  tidy_cliques <- imap_dfr(
    .x = cliques, 
    .f = \(cl, idx) tibble(
      clique_id = idx, 
      tile_id = cl,
      n_tiles = length(cl)
    )
  )
  tidy_cliques
}

# Centre a set of vertices on the origin
centre_vertices <- function(vertices, scale = 1) {
  sweep(vertices, 2, colMeans(vertices)) * scale
}

# Project every row of `vertices` (a matrix of x,y points) onto a direction
# vector; this gives, for each vertex, its signed distance along `direction` 
# (vertices further "ahead" in that direction get larger values). If we set 
# `farthest = TRUE` the function picks out the vertex (or vertices, because
# for our rhombi there will always be two of these) with the largest projection
# along this direction. If we set `farthest = FALSE` it does the opposite,
# and returns the vertices with the smallest projection
extreme_vertices <- function(vertices, direction, farthest = TRUE) {
  proj <- as.numeric(vertices %*% direction)
  extreme <- if (farthest) max(proj) else min(proj)
  vertices[abs(proj - extreme) < 1e-8, , drop = FALSE]
}

# Slot the neighbour tile ("nbr") into place immediately next to the current
# tile ("cur"), using only the information available at placement time. That
# includes the properties for both tiles (`props_cur` and `props_nbr`), as 
# well as the already-placed location of the current tile (`placed_cur`) and
# the "centred" location for the neighbour tile (`centred_nbr`; i.e., the 
# neighbour tile centred on the origin).
shift_neighbour <- function(props_cur, props_nbr, placed_cur, centred_nbr) {
  
  # Each tile sits on exactly two grid lines (its two "sides"); a pair of
  # neighbouring tiles will have one of these in common. Find the shared
  # (family, line) pair and pull out the angle associated with that family. 
  # This is the direction of the grid line that the shared edge lies on.
  sides_cur <- tibble(
    family = c(props_cur$family_1, props_cur$family_2),
    line   = c(props_cur$line_1,   props_cur$line_2),
    angle  = c(props_cur$angle_1,  props_cur$angle_2)
  )
  sides_nbr <- tibble(
    family = c(props_nbr$family_1, props_nbr$family_2),
    line   = c(props_nbr$line_1,   props_nbr$line_2)
  )
  shared_angle <- sides_cur |>
    semi_join(sides_nbr, by = c("family", "line")) |>
    slice(1) |>
    pull(angle)

  # Direction vector `d` that points along the shared grid line
  d <- c(-sin(shared_angle), cos(shared_angle))

  # Decide whether the current tile was originally projected farther
  # along direction `d` (i.e., `sgn > 0`) than the neighbouring tile, 
  # or whether it was projected less far along that direction 
  centre_cur <- c(props_cur$x_mid, props_cur$y_mid)
  centre_nbr <- c(props_nbr$x_mid, props_nbr$y_mid)
  sgn <- sign(sum((centre_nbr - centre_cur) * d))
  if (sgn == 0) sgn <- 1

  # Use extreme_vertices() to detect the "trailing" edge of the already-placed
  # current tile, and the "leading" edge of the to-be-placed neighbouring tile.
  # For our rhombi, which always have edges that are perpendicular to the 
  # relevant direction vector `d`, this will return two points, not one.
  trailing_cur <- extreme_vertices(placed_cur,  d, farthest = sgn > 0)
  leading_nbr  <- extreme_vertices(centred_nbr, d, farthest = sgn < 0)

  # Calculate the shift needed to slide the neighbour tile so that the 
  # midpoint of its "leading" edge lands exactly on the midpoint of the 
  # current tile's "trailing" edge, thereby glueing the two edges together.
  shift_nbr <- colMeans(trailing_cur) - colMeans(leading_nbr)
  shift_nbr
}

# Assemble the complete tiling using breadth-first search
build_tiling <- function(rhombi_tbl, neighbours) {

  # Pull the key information out of the rhombus list
  vertices <- rhombi_tbl$vertices
  props <- rhombi_tbl$properties
  n <- nrow(props)

  # Take each tile's vertex matrix and re-centre it on the origin by 
  # subtracting the column means (i.e. the centroid)
  centred <- map(vertices, centre_vertices)

  # `neighbours` is an edge list (idx, neighbour); split() turns it into a
  # lookup table keyed by tile index, so neighbours_of(i) returns the
  # vector of tile indices adjacent to tile i. split()'s keys become
  # character strings, hence the as.character() when looking one up
  adjacency <- split(neighbours$neighbour, neighbours$idx)
  neighbours_of <- function(i) adjacency[[as.character(i)]]

  # Define a seed tile from which to start the tiling
  seed <- 1

  # Define a `placed` list that will hold the final tiled locations of each
  # of the rhombi. Initially it contains only the seed tile, in centred 
  # coordinates
  placed <- vector("list", n)
  placed[[seed]] <- centred[[seed]]

  # The `seen` vector exists to keep track of which tiles have already been
  # placed: an already-placed tile does not need to be moved again, even if
  # it is a neighbour of the tile currently at the head of the queue
  seen <- rep(FALSE, n)
  seen[seed] <- TRUE

  # Set up for the search: `queue` holds the tiles waiting to have their
  # neighbours placed; `head` specifies the position within the queue that
  # we are currently looking at 
  queue <- seed
  head <- 1

  # Breadth-first search; terminates when there is no new item at the "head"
  # of the search queue. In this context the "head" is the current position
  # in the queue, but the queue can grow by multiple elements in each cycle
  # because one "current" tile might be used to connect more than one new
  # "neighbour" tiles.
  while (head <= length(queue)) {

    # Read the current tile by inspecting the "head" of the queue, then
    # advance the head location forward by 1 in anticipation of the next 
    # cycle
    idx_cur <- queue[head]
    head <- head + 1 

    # Place every neighbour of the current tile
    for (idx_nbr in neighbours_of(idx_cur)) {

      # If we've already seen the neighbour, skip it and move on
      if (seen[idx_nbr]) next

      # Rhombus vertices as they are currently known: the current tile
      # has already been placed, but the neighbour is only centred
      placed_cur  <- placed[[idx_cur]]
      centred_nbr <- centred[[idx_nbr]]

      # Rhombus properties for both tiles
      props_cur <- slice(props, idx_cur)
      props_nbr <- slice(props, idx_nbr)

      # Compute the shift required to move the centred neighbour tile
      # into its placed location in the tiling
      shift_nbr <- shift_neighbour(props_cur, props_nbr, placed_cur, centred_nbr)
      
      # Move the neighbour tile into its placed location
      placed[[idx_nbr]] <- sweep(centred[[idx_nbr]], 2, shift_nbr, "+")

      # Note that the neighbour
      seen[idx_nbr] <- TRUE
      queue <- c(queue, idx_nbr)
    }
  }

  # Return the list of placed vertex matrices, all expressed in the same 
  # edge-matched co-ordinate system. If a tile is never reached (which can
  # happen if it has no neighbours due to how the grid lines crossed the 
  # bounding box earlier) it is left as a NULL entry in the list
  placed
}

# Builds a random instance of Inigo Quilez's cosine-based colour palette
# (https://iquilezles.org/articles/palettes/): each RGB channel varies
# smoothly as `a + b * cos(2*pi*(c*t + d))` for t in [0, 1], where a, b, c,
# d are length-3 vectors. Returns a function of n (the number of colours 
# needed) rather than a fixed palette, so that calling it again with a 
# different n re-samples the same underlying curve.
cosine_palette <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  base <- colors(distinct = TRUE)
  function(n) {
    a <- c(0.5, 0.5, 0.5)
    b <- (sample(base, 1) |> col2rgb() |> as.vector()) / 255
    c <- (sample(base, 1) |> col2rgb() |> as.vector()) / 255
    d <- (sample(base, 1) |> col2rgb() |> as.vector()) / 255
    # evaluate the palette curve at n evenly-spaced points along [0, 1]
    pal <- vapply(
      seq(0, 1, length.out = n), 
      function(t) a + b * cos(2 * pi * (c * t + d)), 
      double(3)
    )
    # the cosine formula can occasionally push a channel outside [0, 1];
    # clip before converting to a valid hex colour
    pal[pal > 1] <- 1
    rgb(t(abs(pal)))
  }
}

# Build the full data structure for one random de Bruijn tiling: draws random
# phase offsets, computes the grid-line intersections, turns each one into a
# rhombus, and glues all rhombi into a single assembled tiling. Returns a
# tidy data frame (one row per tile vertex) ready for ggplot2.
build_debruijn_tiling <- function(n_families, palette, scale) {

  # `scale` is the overall bounding box size; dividing by n_families keeps
  # the *density* of grid lines (and hence the number of tiles) roughly
  # comparable across different n_families, since more line families
  # means more lines crossing any given region. `tile_scale` separately
  # controls the physical size of each individual rhombus and is left
  # fixed at 1. build_tiling() has its own `scale` argument for this,
  # unused here, since build_rhombus() already bakes the size in
  x_scale <- scale / n_families
  y_scale <- scale / n_families
  tile_scale <- 1
 
  # The grid construction pipeline: random phases -> clipped grid lines ->
  # pairwise intersections of those lines (one per rhombus)
  offsets <- sample_offsets(n_families)
  lines <- build_grid_lines(x_scale, y_scale, offsets)
  inter <- find_intersections(lines)

  # Turn every intersection into its own independent rhombus (vertices +
  # metadata); build_rhombus()'s arguments are named to match
  # `intersections`'s columns, so pmap() can call it row-by-row directly
  rhombi_raw <- inter |> 
    mutate(scale = tile_scale) |> 
    pmap(build_rhombus)

  # Re-shape the list of individually-built rhombi into the two parallel
  # structures that compute_neighbours() and build_tiling() expect: a plain
  # list of vertex matrices, and a single properties tibble with one row
  # per tile (tagged with a matching `idx`, 1-based in list order)
  rhombi_tbl <- list(
    vertices = rhombi_raw |> map(\(r) r$vertices),
    properties = rhombi_raw |> 
      imap_dfr(\(r, idx) r$properties, .id = "idx") |> 
      mutate(idx = as.numeric(idx))
  )

  # Figure out which tiles are edge-to-edge neighbours, then glue them all
  # together into one consistent, gap-free assembly
  neighbours <- compute_neighbours(rhombi_tbl)
  placed <- build_tiling(rhombi_tbl, neighbours)

  # Find the cliques of connected tiles in the placed tiling
  cliques <- find_cliques(neighbours)

  # Unique values for the internal acute angles
  bucketed <- round(rhombi_tbl$properties$tile_acute, 4)
  unique_angles <- sort(unique(bucketed))

  # Palettes
  full_palette <- sample(palette(length(unique_angles) + 4))
  tile_palette <- full_palette[seq_along(unique_angles)]
  edge_palette <- full_palette[-seq_along(unique_angles)]

  # Assign tiles
  tile_colour <- tile_palette[match(bucketed, unique_angles)]

  # A handful of rhombi right at the box edge can end up with no shared-line
  # neighbours and are never reached by the BFS; drop those rather than erroring
  reached <- !map_lgl(placed, is.null)
  if (!all(reached)) {
    message(
      sum(!reached), 
      " isolated edge tile(s) dropped from the assembled tiling"
    )
  }

  # Flatten the list of per-tile vertex matrices into one tidy data frame,
  # one row per vertex, ready for plotting
  tiling <- which(reached) |> 
    map_dfr(\(idx) {
      r <- placed[[idx]]
      tibble(
        tile_id = idx,
        vertex = 1:4,
        x = r[, 1],
        y = r[, 2],
        colour = tile_colour[idx]
      )
    })
  
  # Determine the point at which each clique meets
  clique_locs <- cliques |> 
    filter(n_tiles > 2) |> 
    split(~clique_id) |> 
    map_dfr(\(cl) {
      tiling |> 
        filter(tile_id %in% cl$tile_id) |> 
        mutate(x = round(x, 3), y = round(y, 3)) |> 
        count(x, y) |> 
        filter(n == max(n)) |> 
        select(-n) |> 
        mutate(clique_id = cl$clique_id[1])
    })
  
  clique_links <- cliques |> 
    filter(n_tiles > 2) |> 
    left_join(clique_locs, by = "clique_id") |> 
    rename(clique_x = x, clique_y = y) |> 
    left_join(
      tiling |> 
        summarise(.by = tile_id, tile_x = mean(x), tile_y = mean(y)),
      by = "tile_id"
    )
  
  return(list(tiling = tiling, clique = clique_links, edge_palette = edge_palette))
}

# Render an assembled tiling with each rhombus shaded according to its shape
# (based on smallest internal angle)
plot_debruijn_tiling <- function(tiling_info) {

  tiling <- tiling_info$tiling
  clique <- tiling_info$clique
  edge_palette <- tiling_info$edge_palette

  # Horizontal range spanned by the tiling
  x_min <- min(tiling$x)
  x_max <- max(tiling$x)

  # Vertical range spanned by the tiling
  y_min <- min(tiling$y)
  y_max <- max(tiling$y) 

  # Spans along both dimension
  x_extent <- x_max - x_min
  y_extent <- y_max - y_min

  # Decide on a common width and cropping for the plot
  extent <- min(c(x_extent, y_extent))
  crop <- extent * .2

  # Axis limits
  xlim <- c(x_min + crop, x_min + extent - crop)
  ylim <- c(y_min + crop, y_min + extent - crop) 

  # Build and return the plot 
  ggplot() +
    ggforce::geom_shape(
      mapping = aes(x = x, y = y, group = tile_id, fill = colour2),
      data = tiling |> 
        mutate(
          .by = tile_id,
          colour2 = if_else(runif(1) < .5, edge_palette[1], first(colour))
        ),
      colour = edge_palette[2], 
      linewidth = 0.15, 
      expand = -.002
    ) +
    geom_segment(
      mapping = aes(
        x = tile_x, 
        y = tile_y, 
        xend = clique_x,
        yend = clique_y 
      ),
      data = clique,
      color = edge_palette[3]
    ) + 
    geom_point(
      mapping = aes(
        x = clique_x,
        y = clique_y 
      ),
      data = clique |> distinct(clique_id, clique_x, clique_y),
      color = edge_palette[3],
      size = 2
    ) + 
    geom_point(
      mapping = aes(
        x = tile_x, 
        y = tile_y, 
      ),
      data = clique,
      color = edge_palette[3],
      size = 2
    ) +     
    scale_fill_identity() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    coord_equal(xlim = xlim, ylim = ylim) +
    theme_void() + 
    theme(panel.background = element_rect(fill = edge_palette[1]))
}

# Top-level assembly
art_generator <- function(seed, axis_scale = 35) {
  set.seed(seed)

  # Specify the output path and message the user
  output <- output_path(name, version, seed, format)
  message("generating art at ", output)

  # High level parameters
  n_families <- sample(5:10, size = 1)
  palette <- cosine_palette(seed)
 
  # Generate the tiling data structure
  tiling <- build_debruijn_tiling(
    n_families = n_families, 
    palette = palette, 
    scale = axis_scale
  ) 

  # The corresponding ggplot2 object
  p <- plot_debruijn_tiling(tiling)

  # Save result
  ggsave(
    filename = output, 
    plot = p, 
    width = 10, 
    height = 10, 
    units = "in", 
    dpi = 300, 
    bg = "#ffffff00"
  )
}

# Production
seeds <- (version * 100) + 0:99
for (seed in seeds) art_generator(seed)

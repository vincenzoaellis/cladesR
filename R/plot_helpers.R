#' Rectangular phylogram coordinates
#'
#' Internal helper. Returns node, edge, and tip coordinates for a rectangular
#' (right-facing) phylogram, used to lay out the clade-tracker plot.
#'
#' Coordinates come from \code{ape::plotPhyloCoor()}, which returns the same
#' x/y node positions \code{ape::plot.phylo()} would use but without drawing
#' anything. (The earlier approach called \code{ape::plot.phylo(plot = FALSE)}
#' for its side effect and read ape's internal \code{.PlotPhyloEnv}; this is
#' equivalent in output but touches no graphics device and no non-exported
#' environment.)
#'
#' @param tree A \code{phylo} object with branch lengths.
#' @return A list with \code{nodes} (node, x, y), \code{edges} (parent, child
#'   and their x/y), \code{tips} (species, tip, x, y), and \code{tree_height}.
#' @keywords internal
#' @noRd
make_rect_tree_coords <- function(tree) {

  ####################
  ## CHECK PACKAGES ##
  ####################

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }


  ######################
  ## CHECK INPUT DATA ##
  ######################

  if (!inherits(tree, "phylo")) {
    stop("tree must be an object of class 'phylo'.")
  }

  if (is.null(tree$edge.length)) {
    stop("tree must have branch lengths.")
  }


  ############################
  ## GET TREE COORDINATES   ##
  ############################

  ## plotPhyloCoor() returns an (Ntip + Nnode) x 2 matrix of x (= depth from the
  ## root, i.e. node.depth.edgelength) and y positions; row i gives node i.
  coords <- ape::plotPhyloCoor(tree)
  node_x <- coords[, 1]
  node_y <- coords[, 2]

  tree_height <- max(node_x)

  n_tip     <- ape::Ntip(tree)
  all_nodes <- seq_len(n_tip + tree$Nnode)

  node_df <- data.frame(
    node = all_nodes,
    x    = node_x[all_nodes],
    y    = node_y[all_nodes]
  )

  edge_df <- data.frame(
    parent   = tree$edge[, 1],
    child    = tree$edge[, 2],
    x_parent = node_x[tree$edge[, 1]],
    x_child  = node_x[tree$edge[, 2]],
    y_parent = node_y[tree$edge[, 1]],
    y_child  = node_y[tree$edge[, 2]]
  )

  tip_df <- data.frame(
    species = tree$tip.label,
    tip     = seq_along(tree$tip.label),
    x       = node_x[seq_along(tree$tip.label)],
    y       = node_y[seq_along(tree$tip.label)]
  )


  ####################
  ## RETURN RESULTS ##
  ####################

  list(
    nodes       = node_df,
    edges       = edge_df,
    tips        = tip_df,
    tree_height = tree_height
  )
}


#' Truncated tree segments for each time slice
#'
#' Internal helper. For each requested slice age, builds the horizontal and
#' vertical line segments of the tree truncated at that slice, rescaled into a
#' per-slice panel. Used by \code{build_collapsed_clade_data()}.
#'
#' @param tree A \code{phylo} object with branch lengths.
#' @param slice_ages Numeric vector of ages (before present), each strictly
#'   between 0 and the tree height.
#' @param tree_width,gap_width,triangle_width Panel layout widths.
#' @return A list with \code{horizontal}, \code{vertical}, \code{panels}, and
#'   \code{tree_height}.
#' @keywords internal
#' @noRd
build_truncated_tree_segments <- function(tree, slice_ages,
                                          tree_width = 1,
                                          gap_width = 0.65,
                                          triangle_width = 0.25) {

  ####################
  ## CHECK PACKAGES ##
  ####################

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }


  ############################
  ## TREE COORDINATES       ##
  ############################

  coords <- make_rect_tree_coords(tree)
  edge_df <- coords$edges
  tree_height <- coords$tree_height


  ######################
  ## CHECK INPUT DATA ##
  ######################

  ## Slice depth from the root is tree_height - age; an age at or beyond the
  ## tree height gives slice_depth <= 0 and divides by zero in the rescaling.
  if (any(slice_ages <= 0 | slice_ages >= tree_height)) {
    stop("All slice_ages must be strictly between 0 and the tree height (",
         round(tree_height, 4), ").")
  }


  ############################
  ## ORDER SLICES OLD->YOUNG ##
  ############################

  slice_ages <- sort(slice_ages, decreasing = TRUE)

  panel_step <- tree_width + triangle_width + gap_width

  panel_df <- data.frame(
    age = slice_ages,
    panel = seq_along(slice_ages),
    tree_xmin = (seq_along(slice_ages) - 1) * panel_step
  )

  panel_df$tree_xmax <- panel_df$tree_xmin + tree_width
  panel_df$tri_xbase <- panel_df$tree_xmax
  panel_df$tri_xmax <- panel_df$tri_xbase + triangle_width


  ############################
  ## HORIZONTAL SEGMENTS    ##
  ############################

  horiz <- lapply(seq_along(slice_ages), function(i) {

    age <- slice_ages[i]
    slice_depth <- tree_height - age
    panel <- panel_df[panel_df$age == age, ]

    e <- edge_df[edge_df$x_parent < slice_depth, ]
    e$x_child_clip <- pmin(e$x_child, slice_depth)

    out <- data.frame(
      age = age,
      x = panel$tree_xmin + tree_width * (e$x_parent / slice_depth),
      xend = panel$tree_xmin + tree_width * (e$x_child_clip / slice_depth),
      y = e$y_child,
      yend = e$y_child
    )

    out
  })

  horiz <- do.call(rbind, horiz)


  ############################
  ## VERTICAL SEGMENTS      ##
  ############################

  vert <- lapply(seq_along(slice_ages), function(i) {

    age <- slice_ages[i]
    slice_depth <- tree_height - age
    panel <- panel_df[panel_df$age == age, ]

    internal_nodes <- unique(edge_df$parent)

    out <- lapply(internal_nodes, function(nd) {

      child_edges <- edge_df[edge_df$parent == nd, ]

      if (nrow(child_edges) < 2) {
        return(NULL)
      }

      node_x <- child_edges$x_parent[1]

      if (node_x > slice_depth) {
        return(NULL)
      }

      data.frame(
        age = age,
        x = panel$tree_xmin + tree_width * (node_x / slice_depth),
        xend = panel$tree_xmin + tree_width * (node_x / slice_depth),
        y = min(child_edges$y_child),
        yend = max(child_edges$y_child)
      )
    })

    do.call(rbind, out)
  })

  vert <- do.call(rbind, vert)


  ####################
  ## RETURN RESULTS ##
  ####################

  list(
    horizontal = horiz,
    vertical = vert,
    panels = panel_df,
    tree_height = tree_height
  )
}


#' Collapsed clade triangles and connector lines
#'
#' Internal helper. Builds the collapsed-clade triangles (one per clade per
#' slice, area scaled by clade size), their metadata, and the connector lines
#' linking each younger clade to its parent at the previous slice. Consumed by
#' \code{plot_cladetracker()}.
#'
#' @param tree A \code{phylo} object with branch lengths.
#' @param slice_list A list of \code{extract_clade_sizes()} outputs, one per
#'   slice age (must all come from \code{tree}).
#' @param slice_ages Numeric vector of slice ages (one per \code{slice_list}
#'   element).
#' @param tree_width,gap_width,triangle_width,min_triangle_height Layout
#'   parameters.
#' @return A list with \code{tree_horizontal}, \code{tree_vertical},
#'   \code{panels}, \code{triangles}, \code{triangle_meta}, and
#'   \code{connectors}.
#' @keywords internal
#' @noRd
build_collapsed_clade_data <- function(tree, slice_list, slice_ages,
                                       tree_width = 1,
                                       gap_width = 0.65,
                                       triangle_width = 0.25,
                                       min_triangle_height = 0.35) {

  ####################
  ## CHECK PACKAGES ##
  ####################

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }


  ######################
  ## CHECK INPUT DATA ##
  ######################

  if (length(slice_list) != length(slice_ages)) {
    stop("slice_list and slice_ages must be same length.")
  }

  if (!all(vapply(slice_list, function(s) "species" %in% names(s), logical(1)))) {
    stop("Each element of slice_list must be output from extract_clade_sizes().")
  }

  ## All clade species must be tips of the supplied tree, else the coordinate
  ## lookups below silently return NA (slice_list came from a different tree).
  all_slice_species <- unique(unlist(lapply(slice_list, function(s) unlist(s$species))))
  if (!all(all_slice_species %in% tree$tip.label)) {
    stop("slice_list contains species not present in tree; ",
         "slice_list must come from the same tree.")
  }


  ############################
  ## ORDER SLICES OLD->YOUNG ##
  ############################

  ord <- order(slice_ages, decreasing = TRUE)
  slice_list <- slice_list[ord]
  slice_ages <- slice_ages[ord]


  ############################
  ## TREE COORDINATES       ##
  ############################

  coords <- make_rect_tree_coords(tree)
  tip_y  <- stats::setNames(coords$tips$y,  coords$tips$species)
  node_y <- stats::setNames(coords$nodes$y, coords$nodes$node)

  tree_segments <- build_truncated_tree_segments(
    tree = tree,
    slice_ages = slice_ages,
    tree_width = tree_width,
    gap_width = gap_width,
    triangle_width = triangle_width
  )

  panels <- tree_segments$panels


  ############################
  ## OLDEST LINEAGE COLORS  ##
  ############################

  oldest_species <- slice_list[[1]]$species
  oldest_age_val <- slice_ages[1]

  ## For each species set, find which oldest-slice clade contains all of them.
  ## Returns a lineage ID string used for optional coloring in the plot function.
  get_oldest_lineage <- function(sp_set) {
    hits <- sapply(oldest_species, function(x) all(sp_set %in% x))
    hit_names <- names(oldest_species)[hits]
    if (length(hit_names) == 0) return(NA_character_)
    paste0(oldest_age_val, "Ma_", hit_names[1])
  }


  ############################
  ## TRIANGLES              ##
  ############################

  max_size <- max(unlist(lapply(slice_list, function(x) {
    sapply(x$species, length)
  })))

  tri_meta <- list()
  tri_poly <- list()
  counter <- 1

  for (i in seq_along(slice_list)) {

    age <- slice_ages[i]
    sp_list <- slice_list[[i]]$species
    panel <- panels[panels$age == age, ]

    for (clade_name in names(sp_list)) {

      sp_set <- sp_list[[clade_name]]
      ys <- tip_y[sp_set]
      ys <- ys[!is.na(ys)]

      n_sp <- length(sp_set)

      ## Use the MRCA node's ape y-coordinate as ymid so the tree branch
      ## endpoint and triangle base center are on exactly the same y.
      ## For singletons getMRCA() is undefined, so fall back to the tip y.
      if (n_sp > 1) {
        mrca_node <- ape::getMRCA(tree, sp_set)
        ymid <- node_y[[as.character(mrca_node)]]
      } else {
        ymid <- ys[[1]]
      }

      ymin <- min(ys) - 0.20
      ymax <- max(ys) + 0.20

      if ((ymax - ymin) < min_triangle_height) {
        ymin <- ymid - min_triangle_height / 2
        ymax <- ymid + min_triangle_height / 2
      }

      tri_width_i <- triangle_width * sqrt(n_sp / max_size)
      clade_id <- paste0(age, "Ma_", clade_name)
      lineage  <- get_oldest_lineage(sp_set)

      tri_meta[[counter]] <- data.frame(
        age = age,
        clade = clade_name,
        clade_id = clade_id,
        lineage = lineage,
        n_species = n_sp,
        x_base = panel$tri_xbase,
        x_tip = panel$tri_xbase + tri_width_i,
        ymid = ymid,
        ymin = ymin,
        ymax = ymax
      )

      tri_poly[[counter]] <- data.frame(
        age = age,
        clade = clade_name,
        clade_id = clade_id,
        lineage = lineage,
        n_species = n_sp,
        x = c(panel$tri_xbase,
              panel$tri_xbase,
              panel$tri_xbase + tri_width_i),
        y = c(ymin, ymax, ymid)
      )

      counter <- counter + 1
    }
  }

  tri_meta <- do.call(rbind, tri_meta)
  tri_poly <- do.call(rbind, tri_poly)


  ############################
  ## CONNECTORS             ##
  ############################

  connectors <- list()
  counter <- 1

  for (i in seq_len(length(slice_list) - 1)) {

    old_age <- slice_ages[i]
    young_age <- slice_ages[i + 1]

    old_species <- slice_list[[i]]$species
    young_species <- slice_list[[i + 1]]$species

    for (young_clade in names(young_species)) {

      y_sp <- young_species[[young_clade]]

      parent_hits <- sapply(old_species, function(x) {
        all(y_sp %in% x)
      })

      if (any(parent_hits)) {
        old_clade <- names(old_species)[which(parent_hits)[1]]
      } else {
        overlap <- sapply(old_species, function(x) {
          length(intersect(y_sp, x))
        })

        old_clade <- names(which.max(overlap))
      }

      old_tri <- tri_meta[
        tri_meta$age == old_age &
          tri_meta$clade == old_clade,
      ]

      young_tri <- tri_meta[
        tri_meta$age == young_age &
          tri_meta$clade == young_clade,
      ]

      if (nrow(old_tri) == 1 && nrow(young_tri) == 1) {

        connectors[[counter]] <- data.frame(
          connector_id = paste(old_age, old_clade,
                               young_age, young_clade,
                               sep = "_"),
          old_age = old_age,
          young_age = young_age,
          old_clade = old_clade,
          young_clade = young_clade,
          lineage = young_tri$lineage,
          n_species = young_tri$n_species,
          x = old_tri$x_tip,
          xend = young_tri$x_base,
          y = old_tri$ymid,
          yend = young_tri$ymid
        )

        counter <- counter + 1
      }
    }
  }

  connectors <- do.call(rbind, connectors)


  ####################
  ## RETURN RESULTS ##
  ####################

  list(
    tree_horizontal = tree_segments$horizontal,
    tree_vertical = tree_segments$vertical,
    panels = panels,
    triangles = tri_poly,
    triangle_meta = tri_meta,
    connectors = connectors
  )
}

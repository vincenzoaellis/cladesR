#' Clade-tracker plot: collapsed phylogeny panels linked across time slices
#'
#' Draws a multi-panel "clade tracker": for each time slice (oldest on the left)
#' the phylogeny is collapsed into one triangle per clade, with triangle area
#' scaled by clade size, and connector lines link each clade to its parent at
#' the previous (older) slice. The full truncated tree is drawn faintly behind
#' the oldest panel for context. This visualizes how clades split, persist, or
#' diversify across slices --- the ancestor-descendant dynamics central to the
#' diversification-pulse framework (Ricklefs 2014).
#'
#' Two optional color modes (mutually exclusive; \code{heat_splits} wins if
#' both are set):
#' \itemize{
#'   \item \code{color_clades = TRUE} colors triangles/connectors by the
#'     oldest-slice lineage each clade descends from (Okabe-Ito palette).
#'   \item \code{heat_splits = TRUE} colors by the number of descendant clades
#'     a clade splits into at the next slice (blue = passes through, red =
#'     splits the most); \code{scale_by_slice} normalizes this within each panel.
#' }
#'
#' @param tree An ultrametric \code{phylo} object with branch lengths.
#' @param slice_list A list of \code{\link{extract_clade_sizes}} outputs, one
#'   per slice age (all from \code{tree}).
#' @param slice_ages Numeric vector of slice ages (one per \code{slice_list}
#'   element), each strictly between 0 and the tree height.
#' @param tree_width,gap_width,triangle_width Panel layout widths.
#' @param connector_alpha,tree_alpha Opacity of the connector lines and the
#'   background tree.
#' @param color_clades Logical; color by oldest-slice lineage. Default
#'   \code{FALSE}.
#' @param heat_splits Logical; color by number of descendant clades at the next
#'   slice. Default \code{FALSE}. Takes priority over \code{color_clades}.
#' @param scale_by_slice Logical; when \code{heat_splits} is on, normalize the
#'   split counts within each slice panel. Default \code{FALSE}.
#' @param show_clade_labels Logical; label each triangle with its clade ID and
#'   size. Default \code{FALSE}.
#' @param show_split_counts Logical; print the number of descendants at the apex
#'   of each splitting triangle. Default \code{FALSE}.
#' @param flip_triangles Logical; point triangle apices left instead of right.
#'   Default \code{FALSE}.
#' @param tree_color,tree_linewidth Color and line width of the background tree.
#' @param show_slice_ages Logical; label each panel with its slice age. Default
#'   \code{TRUE}.
#' @param show_legend Logical; show the color legend (when a color mode is on).
#'   Default \code{TRUE}.
#' @param legend_title Logical; show the legend title. Default \code{TRUE}.
#' @return A \code{ggplot} object.
#' @references
#' Ricklefs RE (2014) Reconciling diversification: random pulse models of
#' speciation and extinction. Am. Nat. 184:268-276. \doi{10.1086/676642}
#' @seealso \code{\link{extract_clade_sizes}}, \code{\link{ancestor_descendant_stats}}
#' @examples
#' \dontrun{
#' set.seed(1)
#' tree   <- ape::rcoal(40)
#' ages   <- c(0.9, 0.6, 0.3)
#' slices <- lapply(ages, function(a) extract_clade_sizes(tree, a))
#' plot_cladetracker(tree, slices, ages)
#' plot_cladetracker(tree, slices, ages, heat_splits = TRUE)
#' }
#' @export
plot_cladetracker <- function(tree, slice_list, slice_ages,
                              tree_width = 1,
                              gap_width = 0.65,
                              triangle_width = 0.25,
                              connector_alpha = 0.55,
                              tree_alpha = 0.45,
                              color_clades = FALSE,
                              heat_splits = FALSE,
                              scale_by_slice = FALSE,
                              show_clade_labels = FALSE,
                              show_split_counts = FALSE,
                              flip_triangles = FALSE,
                              tree_color = "gray25",
                              tree_linewidth = 0.25,
                              show_slice_ages = TRUE,
                              show_legend = TRUE,
                              legend_title = TRUE) {

  ####################
  ## CHECK PACKAGES ##
  ####################

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }



  ############################
  ## BUILD DATA             ##
  ############################

  plot_data <- build_collapsed_clade_data(
    tree = tree,
    slice_list = slice_list,
    slice_ages = slice_ages,
    tree_width = tree_width,
    gap_width = gap_width,
    triangle_width = triangle_width
  )

  ## Show the full tree only in the oldest (leftmost) panel
  oldest_age <- max(slice_ages)
  h <- plot_data$tree_horizontal[plot_data$tree_horizontal$age == oldest_age, ]
  v <- plot_data$tree_vertical[plot_data$tree_vertical$age == oldest_age, ]
  panels <- plot_data$panels

  triangles     <- plot_data$triangles
  triangle_meta <- plot_data$triangle_meta
  connectors    <- plot_data$connectors

  ## Flip triangles so the apex faces left instead of right.
  ## Reflects each triangle's x-coordinates around its own midpoint.
  if (flip_triangles) {
    triangles$x <- stats::ave(triangles$x, triangles$clade_id,
                       FUN = function(xs) min(xs) + max(xs) - xs)
  }


  ###############################
  ## CONDITIONAL GEOM LAYERS  ##
  ###############################

  ## heat_splits and color_clades are mutually exclusive; heat_splits takes priority.
  if (heat_splits && color_clades) {
    warning("heat_splits and color_clades cannot both be TRUE; using heat_splits.")
    color_clades <- FALSE
  }

  ## Build connector and triangle layers depending on color mode.
  if (heat_splits) {

    ## Count outgoing connectors per source clade per age
    ## (= number of descendants the clade splits into at the next slice).
    n_ch <- stats::aggregate(young_clade ~ old_age + old_clade, data = connectors,
                      FUN = length)
    names(n_ch) <- c("age", "clade", "n_children")

    ## When scale_by_slice = TRUE, normalize n_children within each age panel
    ## so the most-splitting clade in each slice maps to 1. This reveals
    ## within-slice variation even when one early panel dominates globally.
    ## When FALSE, use raw counts so the global maximum drives the color range.
    if (scale_by_slice) {
      ## Normalize within each age panel using (n - 1) / (max - 1) so that
      ## singletons (n = 1, passing through) always map to exactly 0 (blue)
      ## and the most-splitting clade maps to exactly 1 (red).
      ## Edge case: if all clades pass through (max = 1), set fill_value to 0.
      age_max <- tapply(n_ch$n_children, n_ch$age, max)
      n_ch$fill_value <- ifelse(
        age_max[as.character(n_ch$age)] == 1,
        0,
        (n_ch$n_children - 1) / (age_max[as.character(n_ch$age)] - 1)
      )
    } else {
      n_ch$fill_value <- n_ch$n_children
    }

    ## Join fill_value to triangles and connectors.
    ## Clades at the youngest slice have no outgoing connectors (NA -> gray).
    triangles  <- merge(triangles,  n_ch[, c("age", "clade", "fill_value")],
                        by = c("age", "clade"), all.x = TRUE)
    connectors <- merge(connectors, n_ch[, c("age", "clade", "fill_value")],
                        by.x = c("old_age", "old_clade"),
                        by.y = c("age", "clade"), all.x = TRUE)

    conn_layer <- ggplot2::geom_segment(
      data = connectors,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                   color = .data$fill_value),
      alpha = connector_alpha, lineend = "round", linewidth = 0.6
    )

    tri_layer <- ggplot2::geom_polygon(
      data = triangles,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$clade_id,
                   fill = .data$fill_value),
      color = "black", linewidth = 0.20, alpha = 0.90
    )

  } else if (color_clades) {

    ## Okabe-Ito 8-color palette; recycled with a warning if there are more lineages.

    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")

    lineage_levels <- unique(triangles$lineage[!is.na(triangles$lineage)])
    n_lin <- length(lineage_levels)

    if (n_lin > length(okabe_ito)) {
      warning("More lineages (", n_lin, ") than colors (", length(okabe_ito),
              "). Colors will be recycled.")
    }

    fill_colors <- stats::setNames(rep_len(okabe_ito, n_lin), lineage_levels)

    conn_layer <- ggplot2::geom_segment(
      data = connectors,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                   color = .data$lineage),
      alpha = connector_alpha, lineend = "round", linewidth = 0.6
    )

    tri_layer <- ggplot2::geom_polygon(
      data = triangles,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$clade_id,
                   fill = .data$lineage),
      color = "black", linewidth = 0.20, alpha = 0.90
    )

  } else {

    ## Default: gray triangles and connectors, no color mapping.

    conn_layer <- ggplot2::geom_segment(
      data = connectors,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "gray40", alpha = connector_alpha,
      lineend = "round", linewidth = 0.6
    )

    tri_layer <- ggplot2::geom_polygon(
      data = triangles,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$clade_id),
      fill = "gray70", color = "black", linewidth = 0.20, alpha = 0.90
    )

  }


  ############################
  ## BASE PLOT              ##
  ############################

  age_label_layer <- if (show_slice_ages) {
    ggplot2::geom_text(
      data = panels,
      ggplot2::aes(x = .data$tree_xmax, y = Inf,
                   label = paste0(.data$age, " Ma")),
      vjust = 1.5,
      size = 3
    )
  } else {
    NULL
  }

  p <- ggplot2::ggplot() +

    ## horizontal branches
    ggplot2::geom_segment(
      data = h,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = tree_color,
      alpha = tree_alpha,
      linewidth = tree_linewidth
    ) +

    ## vertical branches
    ggplot2::geom_segment(
      data = v,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = tree_color,
      alpha = tree_alpha,
      linewidth = tree_linewidth
    ) +

    conn_layer +

    tri_layer +

    age_label_layer +

    ggplot2::labs(x = NULL, y = NULL) +

    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )


  ##############################
  ## OPTIONAL CLADE LABELS   ##
  ##############################

  ## Labels are centered above each triangle to avoid clipping on the rightmost panel.
  if (show_clade_labels) {

    p <- p +
      ggplot2::geom_text(
        data = triangle_meta,
        ggplot2::aes(x = (.data$x_base + .data$x_tip) / 2, y = .data$ymax,
                     label = paste0(.data$clade, " (", .data$n_species, ")")),
        vjust = -0.3,
        hjust = 0.5,
        size = 2.4
      )
  }


  ##############################
  ## OPTIONAL SPLIT COUNTS    ##
  ##############################

  ## For each triangle that splits into 2+ descendants, print the count of
  ## outgoing connector lines at the triangle apex. Useful when lines are
  ## too numerous to count visually.
  if (show_split_counts) {

    n_children_df <- stats::aggregate(young_clade ~ old_age + old_clade,
                               data = connectors, FUN = length)
    names(n_children_df) <- c("age", "clade", "n_children")

    ## Only label clades that actually split (>= 2 children)
    n_children_df <- n_children_df[n_children_df$n_children >= 2, ]

    if (nrow(n_children_df) > 0) {

      split_labels <- merge(triangle_meta, n_children_df, by = c("age", "clade"))

      ## Place the count just beyond the apex; apex is x_tip normally,
      ## x_base when triangles are flipped to point leftward.
      if (flip_triangles) {
        split_labels$apex_x <- split_labels$x_base
        label_hjust <- 1.2
      } else {
        split_labels$apex_x <- split_labels$x_tip
        label_hjust <- -0.2
      }

      p <- p +
        ggplot2::geom_text(
          data = split_labels,
          ggplot2::aes(x = .data$apex_x, y = .data$ymid, label = .data$n_children),
          hjust = label_hjust,
          size = 2.5
        )
    }
  }



  ##############################
  ## OPTIONAL COLOR SCALES    ##
  ##############################

  if (heat_splits) {

    ## Blue (1 descendant, stagnant) -> yellow -> red (many descendants, splitting).
    ## NA (youngest-panel clades with no outgoing connectors) rendered in gray.
    heat_colors <- c("#4575B4", "#FEE090", "#D73027")

    ## Derive legend title strings (NULL suppresses the title when legend_title = FALSE)
    heat_name <- if (legend_title) {
      if (scale_by_slice) "Splitting intensity\n(scaled per slice)" else "Descendants\nat next slice"
    } else {
      NULL
    }

    if (scale_by_slice) {

      ## Per-slice normalized scale: 0 = fewest splits in this panel, 1 = most.
      p <- p +
        ggplot2::scale_fill_gradientn(
          colors   = heat_colors,
          na.value = "gray80",
          limits   = c(0, 1),
          breaks   = c(0, 0.5, 1),
          labels   = c("min", "mid", "max"),
          name     = heat_name
        ) +
        ggplot2::scale_color_gradientn(
          colors   = heat_colors,
          na.value = "gray60",
          limits   = c(0, 1),
          guide    = "none"
        )

    } else {

      ## Global raw-count scale; integer breaks only.
      p <- p +
        ggplot2::scale_fill_gradientn(
          colors   = heat_colors,
          na.value = "gray80",
          name     = heat_name,
          breaks   = function(lims) seq(ceiling(lims[1]), floor(lims[2]), by = 1)
        ) +
        ggplot2::scale_color_gradientn(
          colors   = heat_colors,
          na.value = "gray60",
          guide    = "none"
        )

    }

  } else if (color_clades) {
    p <- p +
      ggplot2::scale_fill_manual(values = fill_colors, guide = "none") +
      ggplot2::scale_color_manual(values = fill_colors, guide = "none")
  }


  ############################
  ## OPTIONAL LEGEND TOGGLE ##
  ############################

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }


  ####################
  ## RETURN PLOT    ##
  ####################

  p
}

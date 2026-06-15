#' Build a species-by-time-slice clade membership table
#'
#' Assembles a taxonomy-style data frame from a set of tree slices produced by
#' \code{\link{extract_clade_sizes}}: one row per species and one column per
#' time slice (ordered oldest to youngest). Each cell holds the clade ID that
#' the species belongs to at that slice, analogous to a column of taxonomic
#' rank.
#'
#' The result is directly usable for variance partitioning across phylogenetic
#' levels --- a clade-based, equal-aged alternative to a nested ANOVA on Linnaean
#' ranks (Ricklefs 2005). Because the clades are a random sample of possible
#' clades and each level should be tested against the next level down (not the
#' residual), the statistically appropriate model is a mixed model with nested
#' random effects, which yields the variance component at each level directly,
#' e.g. \code{lme4::lmer(trait ~ 1 + (1 | clade_15Ma / clade_10Ma / clade_5Ma))}
#' followed by \code{VarCorr()}. The same decomposition is available in base R
#' via error strata: \code{aov(trait ~ Error(clade_15Ma / clade_10Ma /
#' clade_5Ma))}.
#'
#' @param slice_list A list whose elements are outputs of
#'   \code{\link{extract_clade_sizes}} (each must contain a \code{species}
#'   element).
#' @param slice_ages Numeric vector of slice ages (one per element of
#'   \code{slice_list}, same order). Used to name and order the columns.
#' @return A data frame with a \code{species} column followed by one
#'   \code{clade_<age>Ma} column per slice (oldest first). Rows are sorted
#'   hierarchically (oldest clade, then progressively younger, then species) so
#'   the table reads like a taxonomy.
#' @references
#' Ricklefs RE (2005) Small clades at the periphery of passerine morphological
#' space. Am. Nat. 165:651-659. \doi{10.1086/429676}
#' @seealso \code{\link{extract_clade_sizes}},
#'   \code{\link{ancestor_descendant_stats}}
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(40)
#' ages <- c(1.0, 0.6, 0.3)
#' slices <- lapply(ages, function(a) extract_clade_sizes(tree, a))
#' tab <- build_clade_table(slices, ages)
#' head(tab)
#' @export
build_clade_table <- function(slice_list, slice_ages) {

  ######################
  ## CHECK INPUT DATA ##
  ######################

  if (length(slice_list) != length(slice_ages)) {
    stop("slice_list and slice_ages must be the same length.")
  }

  if (!all(sapply(slice_list, function(s) "species" %in% names(s)))) {
    stop("Each element of slice_list must be output from extract_clade_sizes().")
  }

  ############################
  ## ORDER SLICES OLD->YOUNG ##
  ############################

  ord        <- order(slice_ages, decreasing = TRUE)
  slice_list <- slice_list[ord]
  slice_ages <- slice_ages[ord]

  ####################################
  ## COLLECT ALL SPECIES            ##
  ####################################

  ## All species should be identical at every slice (same tree), but we take
  ## the union to guard against edge cases.
  all_species <- unique(unlist(lapply(slice_list, function(s) unlist(s$species))))

  ######################################
  ## BUILD SPECIES -> CLADE PER SLICE  ##
  ######################################

  col_names  <- paste0("clade_", slice_ages, "Ma")

  clade_cols <- lapply(seq_along(slice_list), function(i) {
    sp_list    <- slice_list[[i]]$species   # named list: clade_id -> species vector
    assignment <- rep(NA_character_, length(all_species))
    names(assignment) <- all_species
    for (clade_id in names(sp_list)) {
      assignment[sp_list[[clade_id]]] <- clade_id
    }
    assignment
  })

  names(clade_cols) <- col_names

  ############################
  ## ASSEMBLE OUTPUT        ##
  ############################

  out <- as.data.frame(
    c(list(species = all_species), clade_cols),
    stringsAsFactors = FALSE
  )

  ## Sort hierarchically: oldest clade first, progressively younger,
  ## then alphabetically by species within the finest grouping.
  ## This gives a display that reads like a taxonomic table.
  out <- out[do.call(order, out[c(col_names, "species")]), ]
  rownames(out) <- NULL

  out
}

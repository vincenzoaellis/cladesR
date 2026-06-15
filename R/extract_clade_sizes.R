#' Slice an ultrametric tree at a given age and extract the resulting clades
#'
#' Cuts a dated, ultrametric phylogeny at a single age before the present and
#' returns the set of clades (subtrees) whose stems cross that slice, together
#' with their sizes, species membership, and the node numbers that root them in
#' the original tree. This is the workhorse that turns a tree into the vector of
#' same-aged clade sizes consumed by \code{\link{clade_tests}},
#' \code{\link{sd_sim_test}}, \code{\link{geom_expectation}}, and
#' \code{\link{clade_rank_data}} (Ricklefs 2014).
#'
#' The slice is specified as an age before the present. Internally the cut is
#' converted to a height measured from the root
#' (\code{slice_from_root = tree_height - age_before_present}) because
#' \code{phytools::treeSlice()} measures from the root. Lineages that do not
#' branch after the slice are retained as singleton clades
#' (\code{trivial = TRUE}), matching the Nee et al. (1992) convention that a
#' non-branching lineage counts as one progeny (itself).
#'
#' When \code{age_before_present == 0} the slice falls at the present, so every
#' tip is its own singleton clade; this case is built directly without calling
#' \code{phytools::treeSlice()}.
#'
#' @param tree An ultrametric \code{phylo} object with branch lengths.
#' @param age_before_present Single non-negative numeric: the age (in the tree's
#'   time units, measured back from the present) at which to slice.
#' @param ultrametric_tol Tolerance passed to \code{ape::is.ultrametric()}.
#'   Default \code{1e-6}.
#' @param ultrametric_option Option passed to \code{ape::is.ultrametric()},
#'   either \code{1} or \code{2}. Default \code{2} (less strict; works with
#'   trees such as those from \code{clootl}).
#' @return A list with five elements:
#'   \describe{
#'     \item{summary}{Data frame of clade sizes, slice metadata, and the
#'       original-tree node number rooting each clade.}
#'     \item{clade_sizes}{Named integer vector of clade sizes (names = clade
#'       IDs); pass directly to the diagnostic functions.}
#'     \item{clades}{A \code{multiPhylo} object of the extracted clade trees.}
#'     \item{species}{Named list of the tip labels in each clade.}
#'     \item{nodes}{Named integer vector of original-tree node numbers (MRCA for
#'       clades of 2+ tips, the tip node for singletons).}
#'   }
#' @references
#' Nee S, Mooers AO, Harvey PH (1992) Tempo and mode of evolution revealed from
#' molecular phylogenies. Proc. Natl. Acad. Sci. USA 89:8322-8326.
#' \doi{10.1073/pnas.89.17.8322}
#'
#' Ricklefs RE (2014) Reconciling diversification: random pulse models of
#' speciation and extinction. Am. Nat. 184:268-276. \doi{10.1086/676642}
#' @seealso \code{\link{clade_tests}}, \code{\link{build_clade_table}},
#'   \code{\link{ancestor_descendant_stats}}
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(40)
#' sliced <- extract_clade_sizes(tree, age_before_present = 0.5)
#' sliced$clade_sizes
#' clade_tests(sliced$clade_sizes)
#' @export
extract_clade_sizes <- function(tree, age_before_present,
                                ultrametric_tol = 1e-6,
                                ultrametric_option = 2) {

  ####################
  ## CHECK PACKAGES ##
  ####################

  ## Check that ape is available
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.")
  }

  ## Check that phytools is available
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is required.")
  }


  ######################
  ## CHECK INPUT DATA ##
  ######################

  ## Check that tree is a phylo object
  if (!inherits(tree, "phylo")) {
    stop("tree must be an object of class 'phylo'.")
  }

  ## Check that tree has branch lengths
  if (is.null(tree$edge.length)) {
    stop("tree must have branch lengths.")
  }

  ## Check that age_before_present is a single numeric value
  if (!is.numeric(age_before_present) || length(age_before_present) != 1) {
    stop("age_before_present must be a single numeric value.")
  }

  ## Check for missing age
  if (is.na(age_before_present)) {
    stop("age_before_present must not be NA.")
  }

  ## Check that age is not negative
  if (age_before_present < 0) {
    stop("age_before_present must be non-negative.")
  }

  ## Check that ultrametric_tol is a single positive numeric value
  if (!is.numeric(ultrametric_tol) || length(ultrametric_tol) != 1 ||
      is.na(ultrametric_tol) || ultrametric_tol <= 0) {
    stop("ultrametric_tol must be a single positive numeric value.")
  } # I don't think people will mess this one up, so a single error seems ok.

  ## Check that ultrametric_option is either 1 or 2
  if (!ultrametric_option %in% c(1, 2)) {
    stop("ultrametric_option must be either 1 or 2.")
  } # I'm going with option 2 because option 1 is too strict with some of my trees (e.g., clootl)


  ############################
  ## CHECK TREE PROPERTIES  ##
  ############################

  ## Check if tree is ultrametric
  is_ultrametric_tree <- ape::is.ultrametric(tree,
    tol = ultrametric_tol,
    option = ultrametric_option)

  if (!is_ultrametric_tree) {
    warning(
      "Tree does not appear ultrametric using ape::is.ultrametric() with ",
      "tol = ", ultrametric_tol, " and option = ", ultrametric_option,
      ". age_before_present probably won't make sense...be careful."
    )
  }


  ############################
  ## CALCULATE SLICE HEIGHT ##
  ############################

  ## get max height (depth) of the tree
  tree_height <- max(phytools::nodeHeights(tree))

  ## phytools::treeSlice() expects the slice height measured from the root.
  ## However, age_before_present is measured backward from the present. So,
  ## slice_from_root = tree_height - age_before_present
  ## For example, if you want 15Ma clades and your tree goes back 50Ma,
  ## you want to supply treeSlice() with 50 - 15 = 35
  slice_from_root <- tree_height - age_before_present

  ## Check that the requested slice falls within the tree.
  if (slice_from_root < 0 || slice_from_root > tree_height) {
    stop("age_before_present is outside the age range of the tree.")
  }


  ################################
  ## SPECIAL CASE: AGE == 0    ##
  ################################

  ## When slicing at the present, every tip is its own singleton clade.
  ## Skip phytools::treeSlice() entirely and build the output directly.
  if (age_before_present == 0) {

    n_clades  <- ape::Ntip(tree)
    n_digits  <- max(2, nchar(as.character(n_clades)))
    clade_names <- paste0("clade_",
      sprintf(paste0("%0", n_digits, "d"), seq_len(n_clades)))

    clade_sizes <- stats::setNames(rep(1L, n_clades), clade_names)

    species <- stats::setNames(
      as.list(tree$tip.label),
      clade_names
    )

    nodes <- stats::setNames(seq_len(n_clades), clade_names)

    clade_trees <- stats::setNames(
      lapply(tree$tip.label, function(sp) ape::keep.tip(tree, sp)),
      clade_names
    )
    class(clade_trees) <- "multiPhylo"

    summary_table <- data.frame(
      clade              = clade_names,
      age_before_present = age_before_present,
      tree_height        = tree_height,
      slice_from_root    = slice_from_root,
      ultrametric_tol    = ultrametric_tol,
      ultrametric_option = ultrametric_option,
      is_ultrametric_tree = is_ultrametric_tree,
      clade_sizes        = clade_sizes,
      node               = nodes,
      row.names          = NULL
    )

    return(list(
      summary     = summary_table,
      clade_sizes = clade_sizes,
      clades      = clade_trees,
      species     = species,
      nodes       = nodes
    ))
  }


  #############################
  ## EXTRACT SUBTREES/CLADES ##
  #############################

  ## Slice the tree at the specified height and return the descendant clades
  ## Note that treeSlice() argument trivial = TRUE to keep singleton lineages
  subtrees <- phytools::treeSlice(tree,
    slice = slice_from_root,
    trivial = TRUE,
    orientation = "tipwards")


  ############################
  ## STANDARDIZE TREE OUTPUT ##
  ############################

  ## Ensure that the sliced clades are stored as a list of phylo objects.
  ## For example, if only one clade is extracted, treeSlice() may return a single
  ## phylo object rather than a list. Wrapping it in list() ensures that the
  ## rest of the function can treat the result as a list element.
  if (inherits(subtrees, "phylo")) {
    subtree_list <- list(subtrees)
  } else {
    subtree_list <- subtrees
  }



  ############################
  ## NAME EXTRACTED CLADES  ##
  ############################

  ## Give each extracted clade a name.
  ## The number of digits is based on the total number of extracted clades,
  ## but names always use at least two digits.
  n_clades <- length(subtree_list)
  n_digits <- max(2, nchar(as.character(n_clades)))

  clade_names <- paste0("clade_",
    sprintf(paste0("%0", n_digits, "d"), seq_len(n_clades)))

  names(subtree_list) <- clade_names # add the names back into the list of trees

  ## Convert the named list of subtree objects to a multiPhylo object.
  clade_trees <- subtree_list
  class(clade_trees) <- "multiPhylo"


  ############################
  ## EXTRACT CLADE CONTENTS ##
  ############################

  ## Count the number of descendant tips in each extracted subtree/clade.
  clade_sizes <- vapply(clade_trees, ape::Ntip, integer(1))

  ## Extract the species/tip labels in each clade.
  clade_species <- lapply(clade_trees, function(x) x$tip.label)

  ## Find the node number in the ORIGINAL tree that roots each extracted clade.
  ## For clades of 2+ species this is the MRCA node (an internal node number
  ## greater than Ntip); for singletons it is the tip node number (1:Ntip).
  ## These node numbers can be used directly to highlight clades on plots of
  ## the original tree (e.g., ape, ggtree, phytools).
  clade_nodes <- vapply(clade_species, function(sp) {
    if (length(sp) > 1) {
      ape::getMRCA(tree, sp)
    } else {
      which(tree$tip.label == sp)
    }
  }, integer(1))


  ##########################
  ## CREATE SUMMARY TABLE ##
  ##########################

  summary_table <- data.frame(
    clade = clade_names,
    age_before_present = age_before_present,
    tree_height = tree_height,
    slice_from_root = slice_from_root,
    ultrametric_tol = ultrametric_tol,
    ultrametric_option = ultrametric_option,
    is_ultrametric_tree = is_ultrametric_tree,
    clade_sizes = clade_sizes,
    node = clade_nodes,
    row.names = NULL
  )


  ############################
  ## RETURN RESULTS         ##
  ############################

  ## Return a list with five outputs:
  ##
  ## 1. summary:
  ##    Data frame of clade sizes, slice metadata, and original-tree node numbers.
  ##
  ## 2. clade_sizes:
  ##    Named integer vector of clade sizes (names = clade IDs).
  ##    Pass directly to clade_tests(), sd_sim_test(), geom_expectation(),
  ##    or clade_rank_data() to preserve clade identity in results.
  ##
  ## 3. clades:
  ##    multiPhylo object containing each extracted clade tree.
  ##
  ## 4. species:
  ##    Named list of species/tip labels in each clade.
  ##
  ## 5. nodes:
  ##    Named integer vector of original-tree node numbers (one per clade).
  ##    Use with ape::extract.clade(), ggtree geom_hilight(), or
  ##    phytools functions to annotate or subset the original tree.
  out <- list(
    summary    = summary_table,
    clade_sizes = clade_sizes,
    clades     = clade_trees,
    species    = clade_species,
    nodes      = clade_nodes
  )

  return(out)
}

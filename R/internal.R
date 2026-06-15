#' Validate a clade_sizes vector
#'
#' Internal helper called at the top of every public function that takes a
#' \code{clade_sizes} argument, so the checks are not repeated inline. A valid
#' vector contains no NAs and only positive whole numbers (clade sizes are
#' counts of species, 1 or greater; a non-branching lineage counts as a clade
#' of size 1, i.e. itself, following Nee et al. 1992).
#'
#' @param clade_sizes A vector of clade sizes.
#' @return Invisibly \code{NULL}; called for its side effect of stopping on
#'   invalid input.
#' @keywords internal
#' @noRd
check_clade_sizes <- function(clade_sizes) {
  if (any(is.na(clade_sizes))) {
    stop("clade_sizes must not contain NAs")
  }
  if (any(clade_sizes < 1)) {
    stop("clade_sizes must be positive (size 1 or greater)")
  }
  if (any(clade_sizes != round(clade_sizes))) {
    stop("clade_sizes must be whole integers (no decimals)")
  }
  invisible(NULL)
}

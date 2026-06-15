#' Summarize ancestor-descendant clade dynamics across adjacent time slices
#'
#' For a clade membership table produced by \code{\link{build_clade_table}},
#' summarizes the relationship between every parent clade at an older slice and
#' its descendant clades at the next (younger) slice. This supports the
#' Ricklefs (2014) test of whether the size of an ancestral clade predicts the
#' sizes of its descendants: under a homogeneous birth-death process ancestor
#' and descendant clade sizes should be effectively uncorrelated, whereas
#' diversification pulses induce a positive association.
#'
#' @param clade_table A data frame from \code{\link{build_clade_table}} (a
#'   \code{species} column plus two or more \code{clade_<age>Ma} columns ordered
#'   oldest to youngest).
#' @return A data frame with one row per parent clade per adjacent slice pair:
#'   \describe{
#'     \item{old_slice, new_slice}{Names of the two slice columns compared.}
#'     \item{parent_clade}{Clade ID at the older slice.}
#'     \item{parent_size}{Number of species in the parent clade.}
#'     \item{n_desc}{Number of descendant clades at the younger slice.}
#'     \item{mean_desc_size}{Mean size of those descendant clades.}
#'     \item{sd_desc_size}{SD of descendant sizes (\code{NA} when \code{n_desc}
#'       is 1).}
#'     \item{max_desc_size}{Size of the largest descendant clade.}
#'   }
#'   Ready to pass directly to \code{stats::cor.test()}, \code{stats::lm()}, or
#'   \code{ggplot2::ggplot()}.
#' @references
#' Ricklefs RE (2014) Reconciling diversification: random pulse models of
#' speciation and extinction. Am. Nat. 184:268-276. \doi{10.1086/676642}
#' @seealso \code{\link{build_clade_table}}, \code{\link{extract_clade_sizes}}
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(40)
#' ages <- c(1.0, 0.6, 0.3)
#' slices <- lapply(ages, function(a) extract_clade_sizes(tree, a))
#' tab <- build_clade_table(slices, ages)
#' ad <- ancestor_descendant_stats(tab)
#' head(ad)
#' @export
ancestor_descendant_stats <- function(clade_table) {

  ######################
  ## CHECK INPUT DATA ##
  ######################

  if (!is.data.frame(clade_table)) {
    stop("clade_table must be a data frame produced by build_clade_table().")
  }

  if (!"species" %in% names(clade_table)) {
    stop("clade_table must have a 'species' column.")
  }

  slice_cols <- setdiff(names(clade_table), "species")

  if (length(slice_cols) < 2) {
    stop("clade_table must contain at least two slice columns.")
  }


  ####################################
  ## COMPUTE STATS PER ADJACENT PAIR ##
  ####################################

  result <- vector("list", length(slice_cols) - 1)

  for (i in seq_len(length(slice_cols) - 1)) {

    old_col <- slice_cols[i]
    new_col <- slice_cols[i + 1]

    parent_ids <- unique(clade_table[[old_col]])

    pair_rows <- lapply(parent_ids, function(pid) {

      ## Species belonging to this parent clade
      in_parent   <- clade_table[[old_col]] == pid
      child_col   <- clade_table[[new_col]][in_parent]

      ## Size of each descendant clade (number of species)
      desc_sizes  <- tabulate(factor(child_col))

      data.frame(
        old_slice      = old_col,
        new_slice      = new_col,
        parent_clade   = pid,
        parent_size    = sum(in_parent),
        n_desc         = length(desc_sizes),
        mean_desc_size = mean(desc_sizes),
        sd_desc_size   = if (length(desc_sizes) > 1) stats::sd(desc_sizes) else NA_real_,
        max_desc_size  = max(desc_sizes),
        row.names      = NULL
      )
    })

    result[[i]] <- do.call(rbind, pair_rows)
  }

  do.call(rbind, result)
}

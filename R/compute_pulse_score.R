#' Per-tip diversification-pulse score across many time slices
#'
#' Computes a continuous "pulse score" for every tip in a tree, measuring how
#' consistently its lineage has belonged to an anomalously sized clade across a
#' range of time slices. This sidesteps the arbitrary choice of a single slice
#' for clade assignment: rather than asking whether a clade is over- or
#' under-sized at one age, it sweeps many ages and asks, per species, how often
#' (and how significantly) that species sat inside a flagged clade.
#'
#' By default it scores membership in **oversized** clades (\code{too_many};
#' the signature of a diversification pulse, Ricklefs 2014). Set
#' \code{direction = "undersized"} to instead score membership in **too-small**
#' clades (\code{too_few}).
#'
#' This is an exploratory **index, not a hypothesis test**. The per-slice
#' p-values are reused only as weights, and slices are non-independent (nested),
#' so \code{pulse_score} should not be interpreted as a significance value. A
#' species in a lineage that has been part of an oversized clade across most of
#' the phylogeny's history scores near 1 on \code{pulse_fraction} and high on
#' \code{pulse_score}; a species in a steadily diversifying lineage scores near
#' 0. Because high scores mark lineages that have repeatedly out-diversified
#' their same-aged peers, the score is expected to track tip-level speciation
#' rate estimators such as Jetz et al.'s (2012) DR statistic, though it measures
#' clade-size anomaly rather than a per-lineage rate directly.
#'
#' @param tree An ultrametric \code{phylo} object with branch lengths.
#' @param slice_ages Numeric vector of positive ages (before present) to
#'   evaluate. Duplicates are dropped and the values are sorted.
#' @param direction Which tail to score: \code{"oversized"} (default; uses
#'   \code{too_many} / \code{p_upper}) or \code{"undersized"} (uses
#'   \code{too_few} / \code{p_lower}).
#' @param alpha Significance threshold passed to \code{\link{clade_tests}} for
#'   the per-slice flags. Default \code{0.05}.
#' @param adjust Multiple-testing correction passed to \code{\link{clade_tests}}
#'   (one of \code{stats::p.adjust.methods}). Default \code{"none"}. When set,
#'   the adjusted tail p-values drive both the flags and the score weights.
#' @param ultrametric_tol,ultrametric_option Passed to
#'   \code{\link{extract_clade_sizes}} at each slice.
#' @return A data frame with one row per tip (original tip order):
#'   \describe{
#'     \item{species}{Tip label.}
#'     \item{n_slices}{Number of slices evaluated.}
#'     \item{n_flagged}{Slices where the species was in a flagged clade
#'       (oversized or undersized, per \code{direction}).}
#'     \item{pulse_fraction}{\code{n_flagged / n_slices}, in [0, 1].}
#'     \item{pulse_score}{Sum of \code{-log10(p)} (tail p-value for the chosen
#'       direction) over flagged slices only; 0 when never flagged. Higher =
#'       more consistently and more significantly anomalous.}
#'     \item{clade_<age>Ma}{Clade ID at each slice (wide format).}
#'     \item{flagged_<age>Ma}{Logical flag at each slice.}
#'   }
#'   The chosen direction is also recorded in \code{attr(x, "direction")}.
#' @references
#' Jetz W, Thomas GH, Joy JB, Hartmann K, Mooers AO (2012) The global diversity
#' of birds in space and time. Nature 491:444-448. \doi{10.1038/nature11631}
#'
#' Ricklefs RE (2014) Reconciling diversification: random pulse models of
#' speciation and extinction. Am. Nat. 184:268-276. \doi{10.1086/676642}
#' @seealso \code{\link{extract_clade_sizes}}, \code{\link{clade_tests}}
#' @examples
#' set.seed(1)
#' tree <- ape::rcoal(40)
#' ps <- compute_pulse_score(tree, slice_ages = c(0.3, 0.6, 0.9))
#' head(ps[order(-ps$pulse_score), c("species", "pulse_fraction", "pulse_score")])
#' @export
compute_pulse_score <- function(tree, slice_ages,
                                direction = c("oversized", "undersized"),
                                alpha = 0.05, adjust = "none",
                                ultrametric_tol = 1e-6,
                                ultrametric_option = 2) {

  ####################
  ## CHECK PACKAGES ##
  ####################

  if (!requireNamespace("ape", quietly = TRUE)) stop("Package 'ape' is required.")
  if (!requireNamespace("phytools", quietly = TRUE)) stop("Package 'phytools' is required.")


  ######################
  ## CHECK INPUT DATA ##
  ######################

  if (!inherits(tree, "phylo")) stop("tree must be an object of class 'phylo'.")
  if (is.null(tree$edge.length)) stop("tree must have branch lengths.")
  if (!is.numeric(slice_ages) || length(slice_ages) < 1) {
    stop("slice_ages must be a non-empty numeric vector.")
  }
  if (any(is.na(slice_ages) | slice_ages <= 0)) {
    stop("slice_ages must be positive and non-NA.")
  }

  direction <- match.arg(direction)

  slice_ages <- sort(unique(slice_ages))

  ## Pick the tail to score. When adjust != "none" the adjusted tail p-value
  ## drives the weight; otherwise the raw tail p-value does. clade_tests()
  ## already bases too_many/too_few on the (adjusted, if requested) values.
  flag_col <- if (direction == "oversized") "too_many" else "too_few"
  p_col    <- if (direction == "oversized") "p_upper"  else "p_lower"
  if (adjust != "none") p_col <- paste0(p_col, "_adj")


  ####################################
  ## EXTRACT CLADES AT EACH SLICE   ##
  ####################################

  all_species <- tree$tip.label

  slice_dfs <- lapply(slice_ages, function(age) {

    sl    <- extract_clade_sizes(tree, age_before_present = age,
                                 ultrametric_tol    = ultrametric_tol,
                                 ultrametric_option = ultrametric_option)
    tests <- clade_tests(sl$clade_sizes, alpha = alpha, adjust = adjust)

    ## invert species list to species -> clade
    sp_to_clade <- data.frame(
      species  = unlist(sl$species),
      clade_id = rep(names(sl$species), lengths(sl$species)),
      row.names = NULL
    )

    ## attach the chosen tail p-value and flag for this direction
    merged <- merge(sp_to_clade,
                    tests[, c("clade", p_col, flag_col)],
                    by.x = "clade_id", by.y = "clade", all.x = TRUE)
    names(merged)[names(merged) == p_col]    <- "tail_p"
    names(merged)[names(merged) == flag_col] <- "flagged"
    merged$slice_age <- age

    ## Weight = -log10(tail p) for flagged clades, else 0. Floor the p-value at
    ## the smallest positive double so an underflow to 0 cannot produce Inf.
    safe_p <- pmax(merged$tail_p, .Machine$double.xmin)
    merged$pulse_weight <- ifelse(merged$flagged, -log10(safe_p), 0)
    merged
  })

  long_df <- do.call(rbind, slice_dfs)


  ####################################
  ## SUMMARIZE PER SPECIES          ##
  ####################################

  summary_df <- do.call(rbind, lapply(split(long_df, long_df$species), function(d) {
    data.frame(
      species        = d$species[1],
      n_slices       = nrow(d),
      n_flagged      = sum(d$flagged, na.rm = TRUE),
      pulse_fraction = mean(d$flagged, na.rm = TRUE),
      pulse_score    = sum(d$pulse_weight, na.rm = TRUE),
      row.names      = NULL
    )
  }))

  ## wide-format per-slice columns (clade ID and flag per age)
  for (age in slice_ages) {
    age_label <- paste0(age, "Ma")
    sub       <- long_df[long_df$slice_age == age, c("species", "clade_id", "flagged")]
    names(sub)[names(sub) == "clade_id"] <- paste0("clade_",   age_label)
    names(sub)[names(sub) == "flagged"]  <- paste0("flagged_", age_label)
    summary_df <- merge(summary_df, sub, by = "species", all.x = TRUE)
  }

  ## restore original tip order and record the scoring direction
  out <- summary_df[match(all_species, summary_df$species), ]
  rownames(out) <- NULL
  attr(out, "direction") <- direction
  out
}

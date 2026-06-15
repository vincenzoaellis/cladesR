#' Test individual clades against the birth-death (geometric) expectation
#'
#' For a set of same-aged clades, tests whether each clade contains more or
#' fewer species than expected under a homogeneous birth-death
#' (speciation-extinction) process. Under such a process the sizes of equal-aged
#' clades follow a geometric distribution whose mean equals the observed mean
#' clade size (Kendall 1948; Nee et al. 1992; Ricklefs 2003). The geometric
#' parameter is estimated as \eqn{p = 1 / \bar{K}}, the maximum-likelihood
#' estimate of Nee et al. (1992) ("the inverse of the average number of progeny
#' lineages").
#'
#' Only one-sided p-values are reported. A combined two-sided p-value is not
#' computed because the geometric PMF is strictly monotone decreasing (its mode
#' is always at \eqn{K = 1}): any "exact" two-sided sum reduces to the upper
#' tail, and the common "double the smaller tail" approximation misleadingly
#' flags singletons as significant when the mean is large. Use
#' \code{p_lower}/\code{too_few} (clade smaller than expected) and
#' \code{p_upper}/\code{too_many} (clade larger than expected) directly.
#'
#' Two nulls are available through \code{method}:
#' \describe{
#'   \item{\code{"geometric"} (default)}{p-values from the \emph{marginal}
#'     geometric distribution with \eqn{p = 1 / \bar{K}}. This is the approach of
#'     Ricklefs (2003, 2014). It treats clades as marginally geometric and
#'     therefore ignores the mild dependence induced by fixing the total
#'     richness, but is fast and closed-form.}
#'   \item{\code{"conditional"}}{the exact conditional (broken-stick) null of Nee
#'     et al. (1992), which conditions on \emph{both} the number of clades
#'     \eqn{n} and the total richness \eqn{S}. Under a geometric process,
#'     conditioning on the sum makes every composition of \eqn{S} into \eqn{n}
#'     positive parts equally likely (the parameter \eqn{p} drops out), so the
#'     marginal size of a single clade has the exact tail
#'     \eqn{P(K \ge k) = \binom{S-k}{\,n-1} / \binom{S-1}{\,n-1}}. This is the
#'     distribution Nee et al. (1992) used to flag per-clade outliers; no
#'     parameter is estimated.}
#' }
#' The two methods agree well when there are many clades; the conditional null
#' matters most when \eqn{n} is small, so that fixing \eqn{S} constrains the
#' sizes appreciably. Either way, with many clades consider \code{adjust} to
#' control the false-discovery rate.
#'
#' @param clade_sizes A vector of clade sizes (positive whole numbers). If named
#'   (for example the \code{clade_sizes} element returned by
#'   \code{extract_clade_sizes()}), the names are carried through as clade IDs.
#' @param alpha Significance threshold for the \code{too_few}/\code{too_many}
#'   flags. Default \code{0.05}.
#' @param method The null distribution for the per-clade tail probabilities. One
#'   of \code{"geometric"} (the default; marginal geometric tail, Ricklefs 2003,
#'   2014) or \code{"conditional"} (the exact conditional broken-stick tail of
#'   Nee et al. 1992, conditioning on both \eqn{n} and \eqn{S}). See Details.
#' @param adjust Multiple-testing correction applied across clades, passed to
#'   \code{\link[stats]{p.adjust}}. One of \code{stats::p.adjust.methods}
#'   (e.g. \code{"BH"}, \code{"bonferroni"}); default \code{"none"}. When not
#'   \code{"none"}, adjusted columns \code{p_lower_adj} and \code{p_upper_adj}
#'   are added and the \code{too_few}/\code{too_many} flags are based on the
#'   adjusted values. The lower-tail and upper-tail families are corrected
#'   separately.
#' @return A \code{data.frame} with one row per clade, sorted from smallest to
#'   largest, with columns: \code{clade} (ID or index), \code{clade_sizes},
#'   \code{mean_clade_size}, \code{p_lower} (P(K <= observed)), \code{p_upper}
#'   (P(K >= observed)), optionally \code{p_lower_adj}/\code{p_upper_adj} (when
#'   \code{adjust} is not \code{"none"}), \code{too_few}, and \code{too_many}.
#'   The null used is recorded in \code{attr(result, "method")}.
#' @references
#' Nee S, Mooers AO, Harvey PH (1992) Tempo and mode of evolution revealed from
#' molecular phylogenies. Proc. Natl. Acad. Sci. USA 89:8322-8326.
#' \doi{10.1073/pnas.89.17.8322}
#'
#' Ricklefs RE (2003) Global diversification rates of passerine birds. Proc. R.
#' Soc. B 270:2285-2291.
#'
#' Ricklefs RE (2014) Reconciling diversification: random pulse models of
#' speciation and extinction. Am. Nat. 184:268-276. \doi{10.1086/676642}
#' @seealso \code{\link{sd_sim_test}}, \code{\link{geom_expectation}},
#'   \code{\link{clade_rank_data}}
#' @examples
#' clade_sizes <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3,
#'                  4, 4, 5, 6, 7, 8, 10, 12, 20, 27)
#' clade_tests(clade_sizes)
#'
#' ## with a Benjamini-Hochberg correction across clades
#' clade_tests(clade_sizes, adjust = "BH")
#'
#' ## the exact conditional (broken-stick) null of Nee et al. 1992
#' clade_tests(clade_sizes, method = "conditional")
#' @export
clade_tests <- function(clade_sizes, alpha = 0.05, adjust = "none",
                        method = c("geometric", "conditional")) {

  ######################
  ## CHECK INPUT DATA ##
  ######################

  check_clade_sizes(clade_sizes)

  method <- match.arg(method)

  if (length(adjust) != 1 || !adjust %in% stats::p.adjust.methods) {
    stop("adjust must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "))
  }

  ##########################################################
  ## MEAN CLADE SIZE AND GEOMETRIC DISTRIBUTION PARAMETER ##
  ##########################################################

  mean_size <- mean(clade_sizes) # observed mean clade size
  p <- 1 / mean_size # geometric parameter; the Nee et al. (1992) MLE. Small p
                     # corresponds to larger expected clade sizes; expected
                     # K = 1/p, so p = 1 / mean.

  ######################################
  ## SINGLE-TAILED P VALUES PER CLADE ##
  ######################################

  if (method == "geometric") {

    # Marginal geometric tails (Ricklefs 2003, 2014). R's geometric
    # distribution counts failures before the first success: X = K - 1, where
    # K is clade size.
    lower_p <- stats::pgeom(clade_sizes - 1, prob = p, lower.tail = TRUE)  # P(K <= observed)
    upper_p <- stats::pgeom(clade_sizes - 2, prob = p, lower.tail = FALSE) # P(K >= observed)

  } else {

    # Exact conditional (broken-stick) tails (Nee et al. 1992). Conditioning a
    # geometric process on the total richness S makes every composition of S
    # into n positive parts equally likely, so the marginal size K of a single
    # clade has the exact upper tail
    #     P(K >= k) = choose(S - k,     n - 1) / choose(S - 1, n - 1)
    # and lower tail
    #     P(K <= k) = 1 - P(K >= k + 1)
    #               = 1 - choose(S - k - 1, n - 1) / choose(S - 1, n - 1).
    # Computed in log space via lchoose for numerical stability; impossible
    # sizes (k > S - n + 1) give lchoose = -Inf, hence a tail of exactly 0.
    n_clades      <- length(clade_sizes)
    total_species <- sum(clade_sizes)

    if (n_clades == 1L) {
      # A single clade is the whole tree: the conditional null is degenerate
      # (only one possible composition), so no clade can be an outlier.
      lower_p <- 1
      upper_p <- 1
    } else {
      log_denom <- lchoose(total_species - 1L, n_clades - 1L)
      upper_p <- exp(lchoose(total_species - clade_sizes,      n_clades - 1L) - log_denom)
      lower_p <- 1 - exp(lchoose(total_species - clade_sizes - 1L, n_clades - 1L) - log_denom)
      # Guard against tiny floating-point excursions outside [0, 1].
      upper_p <- pmin(pmax(upper_p, 0), 1)
      lower_p <- pmin(pmax(lower_p, 0), 1)
    }
  }

  ##############################
  ## OPTIONAL MULTIPLE TESTING ##
  ##############################

  ## When adjust != "none", correct the lower- and upper-tail families
  ## separately and base the significance flags on the adjusted values.
  if (adjust != "none") {
    lower_p_adj <- stats::p.adjust(lower_p, method = adjust)
    upper_p_adj <- stats::p.adjust(upper_p, method = adjust)
    too_few  <- lower_p_adj < alpha
    too_many <- upper_p_adj < alpha
  } else {
    too_few  <- lower_p < alpha
    too_many <- upper_p < alpha
  }

  ####################
  ## OUTPUT RESULTS ##
  ####################

  ## Use clade names if present (e.g. from extract_clade_sizes()), else index.
  clade_ids <- if (!is.null(names(clade_sizes))) names(clade_sizes) else seq_along(clade_sizes)

  result <- data.frame(
    clade           = clade_ids,
    clade_sizes     = clade_sizes,
    mean_clade_size = mean_size,
    p_lower         = lower_p,
    p_upper         = upper_p,
    row.names       = NULL
  )

  if (adjust != "none") {
    result$p_lower_adj <- lower_p_adj
    result$p_upper_adj <- upper_p_adj
  }

  result$too_few  <- too_few
  result$too_many <- too_many

  ## Sort from smallest to largest so unusual clades are easy to find.
  result <- result[order(result$clade_sizes), ]

  ## Record which null was used so downstream code/reporting can tell.
  attr(result, "method") <- method
  result
}

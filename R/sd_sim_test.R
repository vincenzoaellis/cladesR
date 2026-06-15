#' Broken-stick test for over- or under-dispersed clade sizes
#'
#' Tests whether the observed sample standard deviation of clade sizes is
#' unusually large or small relative to the conditional null of Nee et al.
#' (1992). The null fixes both the number of clades (\eqn{n}) and the total
#' species richness (\eqn{S}) --- both are observed facts about the phylogeny,
#' not random quantities --- and asks: given exactly these constraints, is the
#' spread of clade sizes more uneven than expected by chance?
#'
#' Null datasets are generated with the **broken-stick** method of Nee et al.
#' (1992): a stick of length \eqn{S} is broken into \eqn{n} fragments, with
#' breaks allowed only at unit boundaries (here, draw \eqn{n - 1} distinct
#' breakpoints from \eqn{\{1, \dots, S - 1\}}, sort them, and take successive
#' differences). This samples uniformly from all compositions of \eqn{S} into
#' \eqn{n} positive integers, which is mathematically equivalent to conditioning
#' independent Geometric draws on their sum equalling \eqn{S} (the conditional
#' distribution is uniform over compositions regardless of the geometric
#' parameter \eqn{p}). No parameters are estimated. As Nee et al. (1992) put it,
#' "if progeny numbers are geometrically distributed, then all vectors of
#' progeny number are equally probable, as long as the elements sum to the
#' observed total number of progeny."
#'
#' A standard deviation larger than expected indicates among-clade heterogeneity
#' in diversification rate (the signature of diversification pulses; Ricklefs
#' 2014).
#'
#' @param clade_sizes A vector of clade sizes (positive whole numbers).
#' @param nsim Number of broken-stick null datasets to simulate. Default
#'   \code{10000}.
#' @param seed Optional integer seed for reproducibility. Default \code{NULL}
#'   (does not touch the random state).
#' @return A one-row \code{data.frame} with: \code{n_clades}, \code{total_species},
#'   \code{mean_size}, \code{observed_sd}, \code{mean_simulated_sd},
#'   \code{p_sd_too_large}, \code{p_sd_too_small}, and \code{p_two_sided}.
#'   p-values use the \code{(count + 1) / (nsim + 1)} convention so they are
#'   never exactly zero.
#' @references
#' Nee S, Mooers AO, Harvey PH (1992) Tempo and mode of evolution revealed from
#' molecular phylogenies. Proc. Natl. Acad. Sci. USA 89:8322-8326.
#' \doi{10.1073/pnas.89.17.8322}
#'
#' Ricklefs RE (2014) Reconciling diversification: random pulse models of
#' speciation and extinction. Am. Nat. 184:268-276. \doi{10.1086/676642}
#' @seealso \code{\link{clade_tests}}, \code{\link{geom_expectation}}
#' @examples
#' clade_sizes <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3,
#'                  4, 4, 5, 6, 7, 8, 10, 12, 20, 27)
#' sd_sim_test(clade_sizes, nsim = 2000, seed = 1)
#' @export
sd_sim_test <- function(clade_sizes, nsim = 10000, seed = NULL) {

  ######################
  ## CHECK INPUT DATA ##
  ######################

  check_clade_sizes(clade_sizes)

  n_clades      <- length(clade_sizes)
  total_species <- sum(clade_sizes)
  mean_size     <- mean(clade_sizes)
  observed_sd   <- stats::sd(clade_sizes)

  if (!is.null(seed)) set.seed(seed)

  ############################
  ## SIMULATE NULL DATASETS  ##
  ############################

  ## Edge cases: n_clades == 1 or total_species == n_clades means only one
  ## possible composition exists; the broken-stick null is degenerate.
  if (n_clades == 1 || total_species == n_clades) {
    warning("Null is degenerate (only one possible composition). ",
            "p-values are not meaningful.")
    simulated_sd <- rep(0, nsim)
  } else {
    simulated_sd <- replicate(nsim, {
      bp    <- sort(sample.int(total_species - 1L, n_clades - 1L, replace = FALSE))
      sizes <- diff(c(0L, bp, total_species))
      stats::sd(sizes)
    })
  }

  ################
  ## P-VALUES   ##
  ################

  ## +1 correction to numerator and denominator prevents p-values of exactly 0.
  p_sd_too_large <- (sum(simulated_sd >= observed_sd) + 1) / (nsim + 1)
  p_sd_too_small <- (sum(simulated_sd <= observed_sd) + 1) / (nsim + 1)
  p_two_sided    <- min(1, 2 * min(p_sd_too_large, p_sd_too_small))

  ####################
  ## OUTPUT RESULTS ##
  ####################

  data.frame(
    n_clades          = n_clades,
    total_species     = total_species,
    mean_size         = mean_size,
    observed_sd       = observed_sd,
    mean_simulated_sd = mean(simulated_sd),
    p_sd_too_large    = p_sd_too_large,
    p_sd_too_small    = p_sd_too_small,
    p_two_sided       = p_two_sided
  )
}

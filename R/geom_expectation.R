#' Summarize observed clade sizes against the geometric expectation
#'
#' Summarizes the observed distribution of clade sizes and compares it to a
#' homogeneous birth-death (geometric) null fitted with the same mean clade
#' size. Under that null the geometric parameter is \eqn{p = 1 / \bar{K}} (the
#' MLE of Nee et al. 1992) and the population standard deviation of clade sizes
#' is \eqn{\sqrt{\bar{K}(\bar{K} - 1)}} (Ricklefs 2014). A ratio of observed to
#' expected SD/mean above 1 indicates more spread (among-clade rate
#' heterogeneity) than the geometric null predicts.
#'
#' @param clade_sizes A vector of clade sizes (positive whole numbers).
#' @return A one-row \code{data.frame} with: \code{n_clades},
#'   \code{total_species}, \code{mean_size}, \code{geom_p}, \code{observed_sd},
#'   \code{expected_sd}, \code{observed_sd_over_mean}, \code{expected_sd_over_mean},
#'   \code{ratio_obs_over_exp}, \code{simpson_div} (Gini-Simpson index, the
#'   probability that two randomly chosen species belong to different clades),
#'   \code{prop_singleton}, and \code{expected_prop_singleton} (= \code{geom_p}).
#' @references
#' Nee S, Mooers AO, Harvey PH (1992) Tempo and mode of evolution revealed from
#' molecular phylogenies. Proc. Natl. Acad. Sci. USA 89:8322-8326.
#' \doi{10.1073/pnas.89.17.8322}
#'
#' Ricklefs RE (2014) Reconciling diversification: random pulse models of
#' speciation and extinction. Am. Nat. 184:268-276. \doi{10.1086/676642}
#' @seealso \code{\link{clade_tests}}, \code{\link{sd_sim_test}}
#' @examples
#' clade_sizes <- c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3,
#'                  4, 4, 5, 6, 7, 8, 10, 12, 20, 27)
#' geom_expectation(clade_sizes)
#' @export
geom_expectation <- function(clade_sizes) {

  ######################
  ## CHECK INPUT DATA ##
  ######################

  check_clade_sizes(clade_sizes)

  ########################
  ## FIT GEOMETRIC NULL ##
  ########################

  n_clades      <- length(clade_sizes)
  total_species <- sum(clade_sizes)
  mean_size     <- mean(clade_sizes)

  ## Geometric parameter. Clade sizes K are 1, 2, 3, ...; with E[K] = mean_size,
  ## the null's p = 1 / mean_size so the expected clade size equals the observed
  ## mean clade size (Nee et al. 1992 MLE).
  geom_p <- 1 / mean_size

  ####################
  ## OBSERVED STATS ##
  ####################

  observed_sd <- stats::sd(clade_sizes)

  ## Observed ratio of SD to mean; values > 1 indicate SD exceeds the mean.
  observed_sd_over_mean <- observed_sd / mean_size

  ############################
  ## GEOMETRIC EXPECTATIONS ##
  ############################

  ## Theoretical population SD of the fitted geometric distribution.
  ## For K = 1, 2, 3, ... with E[K] = mean_size:
  ## Var(K) = mean_size * (mean_size - 1); SD(K) = sqrt(mean_size * (mean_size - 1)).
  expected_sd <- sqrt(mean_size * (mean_size - 1))

  expected_sd_over_mean <- expected_sd / mean_size

  ############################
  ## ADDITIONAL DESCRIPTORS ##
  ############################

  ## Gini-Simpson index: probability that two randomly chosen species belong to
  ## different clades. Computed inline as 1 - sum(p^2) where p is each clade's
  ## share of total species. 0 = all species in one clade; approaches 1 as
  ## species spread evenly across many clades.
  clade_shares <- clade_sizes / total_species
  simpson_div  <- 1 - sum(clade_shares^2)

  ## Proportion of singleton clades (size = 1).
  ## Under the geometric null, E[prop_singleton] = geom_p = 1 / mean_size.
  prop_singleton <- sum(clade_sizes == 1) / n_clades

  ####################
  ## RETURN RESULTS ##
  ####################

  data.frame(
    n_clades                = n_clades,
    total_species           = total_species,
    mean_size               = mean_size,
    geom_p                  = geom_p,
    observed_sd             = observed_sd,
    expected_sd             = expected_sd,
    observed_sd_over_mean   = observed_sd_over_mean,
    expected_sd_over_mean   = expected_sd_over_mean,
    ratio_obs_over_exp      = observed_sd_over_mean / expected_sd_over_mean,
    simpson_div             = simpson_div,
    prop_singleton          = prop_singleton,
    expected_prop_singleton = geom_p
  )
}

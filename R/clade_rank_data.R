#' Clade rank data for ln(rank)-versus-size plots
#'
#' Ranks clades from largest to smallest and returns the raw and log-transformed
#' clade ranks for making ln(rank)-versus-clade-size plots, together with the
#' rank predicted under the geometric null. Optionally draws the plot.
#'
#' Under a homogeneous birth-death process the log of a clade's rank (rank 1 =
#' largest clade) is a linear function of clade size (Nee et al. 1992; Ricklefs
#' 2014). When clades are ranked largest to smallest, the rank of a clade of
#' size \eqn{k} is the expected number of clades of size \eqn{\geq k}, namely
#' \eqn{n \cdot P(K \geq k) = n (1 - p)^{k - 1}}, so
#' \eqn{\ln(\mathrm{rank}(k)) = \ln(n) + (k - 1)\ln(1 - p)} with
#' \eqn{p = 1 / \bar{K}}. A change in slope of the observed relationship signals
#' departure from the geometric null (e.g. a diversification pulse).
#'
#' @param clade_sizes A vector of clade sizes (positive whole numbers). If named,
#'   the names are carried through as clade IDs.
#' @param plot Logical; if \code{TRUE}, draw a ggplot of ln(rank) versus clade
#'   size with the geometric prediction line. Requires the \pkg{ggplot2} package.
#'   Default \code{FALSE}.
#' @return A \code{data.frame} (sorted largest to smallest) with columns:
#'   \code{clade}, \code{clade_sizes}, \code{clade_rank}, \code{ln_clade_rank}
#'   (natural log of rank), and \code{ln_rank_predicted} (geometric-null
#'   prediction).
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
#' clade_rank_data(clade_sizes)
#' \dontrun{
#' clade_rank_data(clade_sizes, plot = TRUE)
#' }
#' @export
clade_rank_data <- function(clade_sizes, plot = FALSE) {

  ######################
  ## CHECK INPUT DATA ##
  ######################

  check_clade_sizes(clade_sizes)

  ######################
  ## RANK CLADE SIZES ##
  ######################

  ## Sort clades from largest to smallest. Rank 1 = largest clade.
  sort_order   <- order(clade_sizes, decreasing = TRUE)
  sorted_sizes <- clade_sizes[sort_order]
  clade_rank   <- seq_along(sorted_sizes)

  ############
  ## OUTPUT ##
  ############

  ## Geometric null parameters for the predicted ranks (see Details).
  n_clades <- length(clade_sizes)
  geom_p   <- 1 / mean(clade_sizes)

  ## Use clade names if present, otherwise the original numeric index.
  clade_ids <- if (!is.null(names(clade_sizes))) names(clade_sizes)[sort_order] else sort_order

  out <- data.frame(
    clade             = clade_ids,
    clade_sizes       = sorted_sizes,
    clade_rank        = clade_rank,
    ln_clade_rank     = log(clade_rank),
    ln_rank_predicted = log(n_clades) + (sorted_sizes - 1) * log(1 - geom_p),
    row.names         = NULL
  )

  ###################
  ## OPTIONAL PLOT ##
  ###################

  if (plot) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for plotting. Install it or use plot = FALSE.")
    }

    ## Smooth prediction line across all integer clade sizes from 1 to max.
    k_seq        <- seq(1, max(clade_sizes))
    ln_rank_pred <- log(n_clades) + (k_seq - 1) * log(1 - geom_p)

    ## Keep only the portion where predicted rank >= 1 (ln(rank) >= 0).
    pred_df <- data.frame(clade_sizes  = k_seq,
                          ln_rank_pred = ln_rank_pred)
    pred_df <- pred_df[pred_df$ln_rank_pred >= 0, ]

    p <- ggplot2::ggplot(out, ggplot2::aes(x = .data$clade_sizes, y = .data$ln_clade_rank)) +
      ggplot2::geom_line(data = pred_df,
                         ggplot2::aes(x = .data$clade_sizes, y = .data$ln_rank_pred),
                         linetype = "dashed") +
      ggplot2::geom_point(shape = 1) +
      ggplot2::labs(x = "Clade size", y = "ln(Clade rank)") +
      ggplot2::theme_bw(base_size = 12)

    print(p)
  }

  out
}

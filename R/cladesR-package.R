#' @keywords internal
"_PACKAGE"

## The `.data` pronoun is used inside ggplot2::aes() mappings (e.g. in
## clade_rank_data() and plot_cladetracker()) so that R CMD check does not flag
## bare column names as undefined globals. ggplot2 lives in Suggests, so we
## cannot importFrom it; instead we register `.data` as a known global here.
utils::globalVariables(".data")

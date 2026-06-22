#' Tidy table of sensitivity indices
#'
#' @param result An [ospgsa_result].
#' @param output,method,index Optional filters.
#' @param scalar_only Drop time-resolved rows (default `TRUE`).
#' @param digits Rounding for the numeric columns.
#' @return A `data.table` ordered by output / method / index / rank.
#' @export
gsa_table <- function(
  result,
  output = NULL,
  method = NULL,
  index = NULL,
  scalar_only = TRUE,
  digits = 3L
) {
  d <- data.table::copy(result$indices)
  # Build the filter mask in base R: inside d[...] the arguments output/method/
  # index would be shadowed by the like-named columns (data.table NSE).
  keep <- rep(TRUE, nrow(d))
  if (scalar_only) {
    keep <- keep & is.na(d$time)
  }
  if (!is.null(output)) {
    keep <- keep & d$output %in% output
  }
  if (!is.null(method)) {
    keep <- keep & d$method %in% method
  }
  if (!is.null(index)) {
    keep <- keep & d$index %in% index
  }
  d <- d[keep]
  num <- c("estimate", "bias", "std_error", "conf_low", "conf_high")
  for (cn in num) {
    if (cn %in% names(d)) d[[cn]] <- round(d[[cn]], digits)
  }
  data.table::setorder(d, output, method, index, rank)
  d[]
}

#' ospgsa: Global Sensitivity Analysis for PK-Sim / OSP models
#'
#' `ospgsa` performs global sensitivity analysis (GSA) of physiologically based
#' pharmacokinetic (PBPK) models built in PK-Sim and run through the
#' [ospsuite](https://www.open-systems-pharmacology.org/OSPSuite-R/) package.
#'
#' The workflow is:
#' \enumerate{
#'   \item Describe the uncertain model parameters and their distributions with
#'         [gsa_parameter()] / [gsa_parameters()], optionally with a correlation
#'         matrix ([gsa_correlation()]).
#'   \item Build a model evaluator. For PK-Sim use [ospsuite_evaluator()]; for
#'         testing or non-PK-Sim GSA use [function_evaluator()].
#'   \item Run one or more GSA methods with [gsa()] (or the building blocks
#'         [gsa_delta()], [gsa_sobol()], [gsa_morris()], [gsa_regression()]).
#'   \item Inspect results with [print()]/[summary()], [gsa_table()] and the
#'         `plot_*()` family (write the standard set with [gsa_save_plots()]).
#' }
#'
#' The main method is the Borgonovo delta moment-independent index with a
#' two-stage decomposition for correlated inputs (after De Carlo et al. 2023 and
#' Cuquerella-Gilabert et al. 2026). Sobol, Morris and regression methods are
#' provided for cross-checking.
#'
#' @section Design:
#' The GSA engine (sampling, designs, estimators, plots, reports) is
#' model-agnostic: it only needs an evaluator that maps an input matrix to an
#' output matrix. PK-Sim is the first-class evaluator, but any R function can be
#' used, which makes the estimators testable on analytic benchmarks and reusable
#' for other GSA tasks.
#'
#' @keywords internal
#' @importFrom stats density qnorm pnorm qunif qlnorm runif rnorm sd var cor quantile approx coef lm complete.cases setNames cov2cor
#' @importFrom utils head tail modifyList packageVersion
#' @importFrom data.table data.table as.data.table rbindlist setcolorder setorder fcase copy :=
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom mvtnorm rmvnorm
"_PACKAGE"

utils::globalVariables(c(
  ".",
  ".N",
  ".SD",
  ".data",
  "..keep",
  "parameter",
  "output",
  "method",
  "index",
  "value",
  "stage",
  "time",
  "estimate",
  "bias",
  "std_error",
  "conf_low",
  "conf_high",
  "conf_level",
  "significant",
  "rank",
  "mu",
  "mu_star",
  "sigma",
  "S1",
  "ST",
  "delta",
  "QuantityPath",
  "Parameter",
  "Value",
  "Time",
  "est_ind",
  "est_full",
  "d1",
  "d2",
  "sig1",
  "sig2",
  "class"
))

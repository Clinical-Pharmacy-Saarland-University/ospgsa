.sobol_design <- function(parameters, n, method = c("lhs", "random")) {
  method <- match.arg(method)
  k <- length(parameters)
  draw <- function() if (method == "lhs") .lhs(n, k) else matrix(stats::runif(n * k), n, k)
  UA <- draw()
  UB <- draw()
  XA <- .map_marginals(parameters, UA)
  XB <- .map_marginals(parameters, UB)
  blocks <- list(XA, XB)
  for (i in seq_len(k)) {
    Ui <- UA
    Ui[, i] <- UB[, i]
    blocks[[2L + i]] <- .map_marginals(parameters, Ui)
  }
  X <- do.call(rbind, blocks)
  list(X = X, n = n, k = k)
}

.sobol_core <- function(yA, yB, yAB, idx = NULL) {
  if (!is.null(idx)) {
    yA <- yA[idx]
    yB <- yB[idx]
    yAB <- lapply(yAB, function(v) v[idx])
  }
  both <- c(yA, yB)
  VY <- stats::var(both)
  k <- length(yAB)
  Si <- numeric(k)
  STi <- numeric(k)
  if (!is.finite(VY) || VY == 0) {
    return(list(Si = rep(NA_real_, k), STi = rep(NA_real_, k)))
  }
  for (i in seq_len(k)) {
    yABi <- yAB[[i]]
    Si[i] <- mean(yB * (yABi - yA)) / VY # Saltelli (2010)
    STi[i] <- 0.5 * mean((yA - yABi)^2) / VY # Jansen (1999)
  }
  list(Si = Si, STi = STi)
}

#' Sobol first- and total-order sensitivity indices
#'
#' Variance-based indices via the pick-freeze design (Saltelli 2010 first-order,
#' Jansen 1999 total-order) with bootstrap confidence intervals. Total model
#' evaluations: `n * (k + 2)`.
#'
#' **Independence assumption.** Sobol indices are only interpretable for
#' independent inputs; under correlation `sum(Si)` need not be <= 1 and a
#' non-influential parameter can appear important. For correlated PBPK parameters
#' use [gsa_delta_two_stage()]. If `correlation` is supplied here it is ignored
#' with a warning (the design is built from the marginals).
#'
#' @param parameters A [gsa_parameters()] object.
#' @param evaluator A model evaluator: a function mapping a natural-unit input
#'   matrix (columns named by parameter) to an output matrix (see
#'   [function_evaluator()], [ospsuite_evaluator()]).
#' @param n Base sample size (a power of two is conventional).
#' @param boot,conf Bootstrap resamples and confidence level.
#' @param method Marginal sampling: `"lhs"` (default) or `"random"`.
#' @param log_output Analyse `log(Y)` (recommended for AUC / Cmax).
#' @param dummy Optional phantom-parameter name used as the significance floor.
#' @param correlation Ignored (with a warning) if non-`NULL`.
#'
#' @return An [ospgsa_result] with `index` values `"S1"` and `"ST"`.
#' @references Saltelli et al. (2010) *Comput. Phys. Commun.* 181:259;
#'   Jansen (1999) *Comput. Phys. Commun.* 117:35.
#' @examples
#' p <- gsa_parameters(
#'   gsa_parameter("x1", "uniform", min = -pi, max = pi),
#'   gsa_parameter("x2", "uniform", min = -pi, max = pi),
#'   gsa_parameter("x3", "uniform", min = -pi, max = pi))
#' ishigami <- function(M) sin(M[, 1]) + 7 * sin(M[, 2])^2 + 0.1 * M[, 3]^4 * sin(M[, 1])
#' gsa_sobol(p, function_evaluator(ishigami), n = 2048, boot = 50)
#' @export
gsa_sobol <- function(
  parameters,
  evaluator,
  n = 1024L,
  boot = 100L,
  conf = 0.95,
  method = c("lhs", "random"),
  log_output = FALSE,
  dummy = NULL,
  correlation = NULL
) {
  .check_parameters(parameters)
  .check_evaluator(evaluator)
  if (!is.null(correlation)) {
    .warn(c(
      "{.fn gsa_sobol} ignores the supplied {.arg correlation}.",
      "i" = "The Sobol pick-freeze design assumes independent inputs; indices are computed from the marginals.",
      "i" = "For correlated inputs use {.fn gsa_delta_two_stage}."
    ))
  }
  method <- match.arg(method)
  n <- .assert_count(n, min = 4L)
  boot <- as.integer(boot)
  k <- length(parameters)
  pnames <- .param_names(parameters)

  des <- .sobol_design(parameters, n, method)
  Y <- .coerce_outputs(evaluator(des$X))
  if (nrow(Y) != nrow(des$X)) {
    .abort_eval_rows(nrow(Y), nrow(des$X))
  }
  onames <- colnames(Y)

  rows <- list()
  n_failed_total <- 0L
  for (oc in seq_along(onames)) {
    yv <- Y[, oc]
    if (log_output) {
      if (any(yv[is.finite(yv)] <= 0)) {
        .stop(c(
          "{.arg log_output = TRUE} requires strictly positive outputs.",
          "x" = "Output {.val {onames[oc]}} contains non-positive values."
        ))
      }
      yv <- log(yv)
    }
    yA <- yv[seq_len(n)]
    yB <- yv[n + seq_len(n)]
    yAB <- lapply(seq_len(k), function(i) yv[(1L + i) * n + seq_len(n)])
    keep <- is.finite(yA) & is.finite(yB) & Reduce(`&`, lapply(yAB, is.finite))
    n_failed_total <- n_failed_total + sum(!keep)
    if (sum(keep) < 4L) {
      .warn(c(
        "Skipping output {.val {onames[oc]}}.",
        "i" = "Too few valid Sobol rows ({sum(keep)}) after removing failed runs."
      ))
      next
    }
    if (any(!keep)) {
      .warn(c(
        "Output {.val {onames[oc]}}: {sum(!keep)}/{n} pick-freeze row{?s} failed.",
        "!" = "Sobol indices use complete rows only and may be biased; consider re-running or using {.fn gsa_delta}."
      ))
      yA <- yA[keep]
      yB <- yB[keep]
      yAB <- lapply(yAB, function(v) v[keep])
    }
    est <- .sobol_core(yA, yB, yAB)
    Nk <- length(yA)
    bS <- matrix(NA_real_, boot, k)
    bT <- matrix(NA_real_, boot, k)
    if (boot > 0L) {
      for (b in seq_len(boot)) {
        r <- sample.int(Nk, Nk, replace = TRUE)
        cc <- .sobol_core(yA, yB, yAB, idx = r)
        bS[b, ] <- cc$Si
        bT[b, ] <- cc$STi
      }
    }
    for (i in seq_len(k)) {
      ciS <- if (boot > 0L) .boot_ci(bS[, i], conf) else c(NA, NA)
      ciT <- if (boot > 0L) .boot_ci(bT[, i], conf) else c(NA, NA)
      rows[[length(rows) + 1L]] <- .indices_dt(
        output = onames[oc],
        parameter = pnames[i],
        method = "sobol",
        index = "S1",
        estimate = est$Si[i],
        std_error = if (boot > 0L) stats::sd(bS[, i]) else NA_real_,
        conf_low = ciS[[1]],
        conf_high = ciS[[2]],
        conf_level = conf
      )
      rows[[length(rows) + 1L]] <- .indices_dt(
        output = onames[oc],
        parameter = pnames[i],
        method = "sobol",
        index = "ST",
        estimate = est$STi[i],
        std_error = if (boot > 0L) stats::sd(bT[, i]) else NA_real_,
        conf_low = ciT[[1]],
        conf_high = ciT[[2]],
        conf_level = conf
      )
    }
  }
  if (!length(rows)) {
    .stop(c(
      "No Sobol indices could be computed.",
      "i" = "All outputs were skipped (too few valid runs). Increase {.arg n} or check the model."
    ))
  }
  d <- .finalize_indices(data.table::rbindlist(rows, fill = TRUE), dummy = dummy)
  new_ospgsa_result(
    d,
    design = des$X,
    Y = Y,
    meta = list(
      method = "sobol",
      n = n,
      boot = boot,
      conf = conf,
      n_runs = nrow(des$X),
      n_failed = n_failed_total,
      log_output = log_output
    )
  )
}

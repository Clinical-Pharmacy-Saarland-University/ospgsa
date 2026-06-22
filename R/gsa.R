.with_dummy <- function(parameters, dummy_name) {
  gsa_parameters(c(unclass(parameters), list(gsa_dummy(dummy_name))))
}

#' Run a global sensitivity analysis
#'
#' One entry point that builds the appropriate design(s), evaluates the model and
#' computes the requested sensitivity indices. For correlated inputs with
#' `method = "delta"` it automatically runs the **two-stage** analysis
#' (independent + correlated designs). A dummy noise-floor parameter is added by
#' default so significance is calibrated empirically.
#'
#' @param parameters A [gsa_parameters()] object.
#' @param evaluator A model evaluator ([ospsuite_evaluator()] for PK-Sim, or
#'   [function_evaluator()]).
#' @param method One or more of `"delta"` (default, recommended for PBPK),
#'   `"sobol"`, `"morris"`, `"regression"`.
#' @param n Base sample size for `delta`/`regression`/`sobol`. Defaults to 2000
#'   (given-data) or 1024 (sobol).
#' @param correlation Optional correlation matrix ([gsa_correlation()]). Triggers
#'   two-stage delta; ignored by sobol/morris (with a warning).
#' @param corr_type,correlation_method Passed to [gsa_sample()].
#' @param boot,conf,ci Bootstrap settings.
#' @param log_output Analyse `log(Y)` (recommended for AUC / Cmax).
#' @param add_dummy,dummy_name Add a phantom parameter as the significance floor.
#' @param seed Optional RNG seed.
#' @param r,levels Morris trajectory count and grid levels.
#' @param reg_methods Regression indices to compute (see [gsa_regression()]).
#' @param n_cores Workers for the delta estimator's per-task bootstrap (passed to
#'   [gsa_delta()] / [gsa_delta_two_stage()]); `1` (default) is serial. With a
#'   source package, parallel workers must load it (see [gsa_delta()] `cl`).
#'
#' @return An [ospgsa_result] (combined across methods if several are requested).
#' @seealso [gsa_delta_two_stage()], [gsa_sobol()], [gsa_morris()],
#'   [gsa_regression()], [gsa_convergence()]
#' @examples
#' p <- gsa_parameters(
#'   gsa_parameter("x1", "uniform", min = 0, max = 1),
#'   gsa_parameter("x2", "uniform", min = 0, max = 1),
#'   gsa_parameter("x3", "uniform", min = 0, max = 1))
#' ev <- function_evaluator(function(M) 3 * M[, 1] + M[, 2], "y")
#' gsa(p, ev, method = "delta", n = 1500, boot = 40, seed = 1)
#' @export
gsa <- function(
  parameters,
  evaluator,
  method = c("delta", "sobol", "morris", "regression"),
  n = NULL,
  correlation = NULL,
  corr_type = "spearman",
  correlation_method = "iman-conover",
  boot = 100L,
  conf = 0.95,
  ci = "percentile",
  log_output = FALSE,
  add_dummy = TRUE,
  dummy_name = "dummy",
  seed = NULL,
  r = 10L,
  levels = 4L,
  reg_methods = c("SRC", "PRCC"),
  n_cores = 1L
) {
  method <- match.arg(
    tolower(method),
    c("delta", "sobol", "morris", "regression"),
    several.ok = TRUE
  )
  .check_parameters(parameters)
  .check_evaluator(evaluator)

  params <- parameters
  dummy <- .dummy_name(params)
  if (add_dummy && is.null(dummy)) {
    params <- .with_dummy(parameters, dummy_name)
    dummy <- dummy_name
  }
  if (is.null(n)) {
    n <- if ("sobol" %in% method) 1024L else 2000L
  }
  if (!is.null(correlation)) {
    nm <- .param_names(params)
    if (!all(nm %in% rownames(correlation))) {
      R <- diag(length(nm))
      dimnames(R) <- list(nm, nm)
      cn <- intersect(rownames(correlation), nm)
      R[cn, cn] <- correlation[cn, cn]
      correlation <- R
    }
  }
  res <- list()

  if (any(c("delta", "regression") %in% method)) {
    if (!is.null(correlation)) {
      Sf <- gsa_sample(
        params,
        n,
        correlation = correlation,
        corr_type = corr_type,
        correlation_method = correlation_method,
        seed = seed
      )
      Yf <- .coerce_outputs(evaluator(Sf$X))
      if ("delta" %in% method) {
        Si <- gsa_sample(params, n, seed = if (is.null(seed)) NULL else seed + 1L)
        Yi <- .coerce_outputs(evaluator(Si$X))
        res$delta <- gsa_delta_two_stage(
          Sf$X,
          Yf,
          Si$X,
          Yi,
          boot = boot,
          conf = conf,
          ci = ci,
          log_output = log_output,
          dummy = dummy,
          seed = seed,
          n_cores = n_cores
        )
      }
      if ("regression" %in% method) {
        res$reg <- gsa_regression(
          Sf$X,
          Yf,
          methods = reg_methods,
          parameters = params,
          boot = boot,
          conf = conf,
          log_output = log_output,
          dummy = dummy
        )
      }
    } else {
      S <- gsa_sample(params, n, seed = seed)
      Y <- .coerce_outputs(evaluator(S$X))
      if ("delta" %in% method) {
        res$delta <- gsa_delta(
          S$X,
          Y,
          boot = boot,
          conf = conf,
          ci = ci,
          log_output = log_output,
          dummy = dummy,
          seed = seed,
          n_cores = n_cores
        )
      }
      if ("regression" %in% method) {
        res$reg <- gsa_regression(
          S$X,
          Y,
          methods = reg_methods,
          parameters = params,
          boot = boot,
          conf = conf,
          log_output = log_output,
          dummy = dummy
        )
      }
    }
  }
  if ("sobol" %in% method) {
    res$sobol <- gsa_sobol(
      params,
      evaluator,
      n = n,
      boot = boot,
      conf = conf,
      log_output = log_output,
      dummy = dummy,
      correlation = correlation
    )
  }
  if ("morris" %in% method) {
    res$morris <- gsa_morris(
      params,
      evaluator,
      r = r,
      levels = levels,
      boot = boot,
      conf = conf,
      log_output = log_output,
      dummy = dummy
    )
  }

  out <- if (length(res) == 1L) res[[1L]] else .combine_results(res)
  out$meta$dummy <- dummy
  out$meta$methods <- method
  out$meta$correlated <- !is.null(correlation)
  out
}

#' Convergence of GSA indices with sample size
#'
#' Re-estimates the indices at increasing sample sizes to check Monte Carlo
#' convergence (stable rankings, tight CIs).
#'
#' @param parameters,evaluator As in [gsa()].
#' @param method One of `"delta"`, `"sobol"`, `"regression"` (given-data /
#'   variance methods).
#' @param n_seq Increasing sample sizes.
#' @param boot,conf,log_output,seed As in [gsa()].
#' @param index Which index to track (default depends on method: `"delta"`,
#'   `"ST"` or `"PRCC"`).
#' @return A `data.table` (class `ospgsa_convergence`) with columns `n`,
#'   `output`, `parameter`, `index`, `estimate`, `conf_low`, `conf_high`.
#' @export
gsa_convergence <- function(
  parameters,
  evaluator,
  method = c("sobol", "delta", "regression"),
  n_seq = c(256L, 512L, 1024L, 2048L),
  boot = 50L,
  conf = 0.95,
  log_output = FALSE,
  seed = NULL,
  index = NULL
) {
  method <- match.arg(method)
  idx_name <- if (is.null(index)) {
    switch(method, delta = "delta", sobol = "ST", regression = "PRCC")
  } else {
    index
  }
  out <- vector("list", length(n_seq))
  for (i in seq_along(n_seq)) {
    ni <- n_seq[i]
    sd_i <- if (is.null(seed)) NULL else seed + i
    r <- switch(
      method,
      sobol = gsa_sobol(
        parameters,
        evaluator,
        n = ni,
        boot = boot,
        conf = conf,
        log_output = log_output
      ),
      {
        S <- gsa_sample(parameters, ni, seed = sd_i)
        Y <- .coerce_outputs(evaluator(S$X))
        if (method == "delta") {
          gsa_delta(S$X, Y, boot = boot, conf = conf, log_output = log_output)
        } else {
          gsa_regression(
            S$X,
            Y,
            methods = idx_name,
            parameters = parameters,
            boot = boot,
            conf = conf,
            log_output = log_output
          )
        }
      }
    )
    ix <- as.data.frame(r$indices)
    d <- ix[ix$index == idx_name, , drop = FALSE]
    d$n <- ni
    out[[i]] <- d[, c("n", "output", "parameter", "index", "estimate", "conf_low", "conf_high")]
  }
  res <- data.table::rbindlist(out, fill = TRUE)
  data.table::setattr(res, "class", c("ospgsa_convergence", class(res)))
  res[]
}

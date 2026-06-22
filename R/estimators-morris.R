.morris_trajectory <- function(k, p) {
  Delta <- p / (2 * (p - 1))
  grid <- seq(0, 1, by = 1 / (p - 1))
  Dstar <- sample(c(-1, 1), k, replace = TRUE)
  base <- numeric(k)
  for (i in seq_len(k)) {
    cand <- if (Dstar[i] > 0) grid[grid <= 1 - Delta + 1e-9] else grid[grid >= Delta - 1e-9]
    base[i] <- cand[sample.int(length(cand), 1L)]
  }
  perm <- sample.int(k)
  pts <- matrix(0, k + 1L, k)
  pts[1L, ] <- base
  cur <- base
  changed <- integer(k)
  sgn <- numeric(k)
  for (t in seq_len(k)) {
    i <- perm[t]
    cur[i] <- cur[i] + Dstar[i] * Delta
    pts[t + 1L, ] <- cur
    changed[t] <- i
    sgn[t] <- Dstar[i]
  }
  list(pts = pts, changed = changed, sgn = sgn, Delta = Delta)
}

#' Morris elementary-effects screening
#'
#' Computes the Morris screening measures mu (mean effect), mu* (mean absolute
#' effect, a robust importance ranking / proxy for the total effect) and sigma
#' (spread of effects, indicating nonlinearity and interactions), with bootstrap
#' CIs for mu* obtained by resampling trajectories.
#'
#' Elementary effects are computed in the normalized unit-hypercube design space,
#' so the marginal distributions (including log-scaled parameters) are accounted
#' for through the inverse-CDF mapping and effects are comparable across
#' parameters. Use Morris as a cheap first pass to drop clearly inactive factors
#' before a quantitative analysis.
#'
#' @param parameters A [gsa_parameters()] object.
#' @param evaluator A model evaluator (see [function_evaluator()],
#'   [ospsuite_evaluator()]).
#' @param r Number of trajectories (10-50 typical).
#' @param levels Number of grid levels `p` (even; default 4).
#' @param boot,conf Bootstrap resamples / confidence level for mu*.
#' @param log_output Analyse `log(Y)`.
#' @param dummy Optional phantom-parameter name for the significance floor.
#'
#' @return An [ospgsa_result] with `index` values `"mu_star"`, `"mu"`, `"sigma"`
#'   and `"mu_star_norm"` (mu* normalized to sum to 1 within each output).
#' @references Morris (1991) *Technometrics* 33:161; Campolongo, Cariboni &
#'   Saltelli (2007) *Environ. Model. Softw.* 22:1509.
#' @examples
#' p <- gsa_parameters(
#'   gsa_parameter("x1", "uniform", min = 0, max = 1),
#'   gsa_parameter("x2", "uniform", min = 0, max = 1),
#'   gsa_parameter("x3", "uniform", min = 0, max = 1))
#' f <- function(M) 3 * M[, 1] + 2 * M[, 2]               # x3 inert
#' gsa_morris(p, function_evaluator(f), r = 20)
#' @export
gsa_morris <- function(
  parameters,
  evaluator,
  r = 10L,
  levels = 4L,
  boot = 100L,
  conf = 0.95,
  log_output = FALSE,
  dummy = NULL
) {
  .check_parameters(parameters)
  .check_evaluator(evaluator)
  r <- .assert_count(r, min = 2L)
  levels <- .assert_count(levels, min = 2L)
  if (levels %% 2L != 0L) {
    .warn(c(
      "Morris {.arg levels} (the grid resolution {.field p}) is usually even.",
      "i" = "You supplied {.val {levels}}; an even value keeps the step on the grid."
    ))
  }
  boot <- as.integer(boot)
  k <- length(parameters)
  pnames <- .param_names(parameters)

  trajs <- lapply(seq_len(r), function(i) .morris_trajectory(k, levels))
  bigU <- do.call(rbind, lapply(trajs, `[[`, "pts"))
  bigX <- .map_marginals(parameters, bigU)
  Y <- .coerce_outputs(evaluator(bigX))
  if (nrow(Y) != nrow(bigX)) {
    .abort_eval_rows(nrow(Y), nrow(bigX))
  }
  onames <- colnames(Y)
  step <- k + 1L

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
    EE <- matrix(NA_real_, r, k)
    for (ti in seq_len(r)) {
      tr <- trajs[[ti]]
      base_row <- (ti - 1L) * step
      yt <- yv[base_row + seq_len(step)]
      for (t in seq_len(k)) {
        i <- tr$changed[t]
        EE[ti, i] <- (yt[t + 1L] - yt[t]) / (tr$sgn[t] * tr$Delta)
      }
    }
    n_failed_total <- n_failed_total + sum(!is.finite(yv))
    mu <- apply(EE, 2L, function(e) mean(e, na.rm = TRUE))
    mu_star <- apply(EE, 2L, function(e) mean(abs(e), na.rm = TRUE))
    sigma <- apply(EE, 2L, function(e) stats::sd(e, na.rm = TRUE))
    denom <- sum(mu_star, na.rm = TRUE)
    mu_star_norm <- if (denom > 0) mu_star / denom else mu_star * NA

    bM <- matrix(NA_real_, boot, k)
    if (boot > 0L) {
      for (b in seq_len(boot)) {
        rr <- sample.int(r, r, replace = TRUE)
        bM[b, ] <- apply(EE[rr, , drop = FALSE], 2L, function(e) mean(abs(e), na.rm = TRUE))
      }
    }
    for (i in seq_len(k)) {
      ci <- if (boot > 0L) .boot_ci(bM[, i], conf) else c(NA, NA)
      rows[[length(rows) + 1L]] <- .indices_dt(
        output = onames[oc],
        parameter = pnames[i],
        method = "morris",
        index = "mu_star",
        estimate = mu_star[i],
        std_error = if (boot > 0L) stats::sd(bM[, i]) else NA_real_,
        conf_low = ci[[1]],
        conf_high = ci[[2]],
        conf_level = conf
      )
      rows[[length(rows) + 1L]] <- .indices_dt(
        output = onames[oc],
        parameter = pnames[i],
        method = "morris",
        index = "mu",
        estimate = mu[i],
        conf_level = conf
      )
      rows[[length(rows) + 1L]] <- .indices_dt(
        output = onames[oc],
        parameter = pnames[i],
        method = "morris",
        index = "sigma",
        estimate = sigma[i],
        conf_level = conf
      )
      rows[[length(rows) + 1L]] <- .indices_dt(
        output = onames[oc],
        parameter = pnames[i],
        method = "morris",
        index = "mu_star_norm",
        estimate = mu_star_norm[i],
        conf_level = conf
      )
    }
  }
  if (!length(rows)) {
    .stop(c(
      "No Morris indices could be computed.",
      "i" = "All outputs were skipped. Increase {.arg r} or check the model."
    ))
  }
  d <- .finalize_indices(data.table::rbindlist(rows, fill = TRUE), dummy = dummy)
  new_ospgsa_result(
    d,
    design = bigX,
    Y = Y,
    meta = list(
      method = "morris",
      r = r,
      levels = levels,
      boot = boot,
      conf = conf,
      n_runs = nrow(bigX),
      n_failed = n_failed_total,
      log_output = log_output
    )
  )
}

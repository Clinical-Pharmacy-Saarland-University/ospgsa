.src_one <- function(X, y, rankit) {
  if (rankit) {
    X <- apply(X, 2L, rank)
    y <- rank(y)
  }
  df <- as.data.frame(X)
  df$.y <- y
  fit <- stats::lm(.y ~ ., data = df)
  b <- stats::coef(fit)[-1L]
  sdx <- apply(X, 2L, stats::sd)
  sdy <- stats::sd(y)
  src <- b * sdx / sdy
  names(src) <- colnames(X)
  list(src = src, r2 = suppressWarnings(summary(fit)$r.squared))
}

.pcc_all <- function(X, y, rankit) {
  M <- cbind(X, .y = y)
  C <- if (rankit) stats::cor(M, method = "spearman") else stats::cor(M)
  P <- tryCatch(solve(C), error = function(e) MASS_ginv(C))
  yi <- ncol(M)
  pcc <- vapply(seq_len(ncol(X)), function(i) -P[i, yi] / sqrt(P[i, i] * P[yi, yi]), numeric(1))
  names(pcc) <- colnames(X)
  pcc
}

# Moore-Penrose fallback to avoid a hard MASS dependency.
MASS_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  d <- s$d
  ok <- d > max(tol * d[1L], 0)
  s$v[, ok, drop = FALSE] %*% (t(s$u[, ok, drop = FALSE]) / d[ok])
}

#' Regression / partial-correlation sensitivity indices
#'
#' Computes standardized regression coefficients (SRC / SRRC) and partial
#' correlation coefficients (PCC / PRCC) from a single input-output sample, with
#' bootstrap confidence intervals and the model `R^2` (SRC / SRRC are only
#' reliable when `R^2` is high; escalate to [gsa_delta()] / [gsa_sobol()] for
#' non-(rank-)linear models).
#'
#' @param X Input design ([gsa_sample()], matrix or data.frame). For SRC, columns
#'   are transformed to the per-parameter estimator scale (e.g. log) when a
#'   [gsa_sample()] / `parameters` is supplied.
#' @param Y Output(s) (vector, matrix, data.frame or evaluation result).
#' @param methods Which indices to compute: any of `"SRC"`, `"SRRC"`, `"PCC"`,
#'   `"PRCC"`. Default `c("SRC", "PRCC")`.
#' @param parameters Optional [gsa_parameters()] used to apply per-parameter
#'   scale transforms for SRC (ignored for rank methods).
#' @param boot,conf Bootstrap resamples / confidence level.
#' @param log_output Analyse `log(Y)`.
#' @param dummy Optional phantom-parameter name for the significance floor.
#'
#' @return An [ospgsa_result]; `index` is one of the requested methods plus an
#'   `"R2"` row (parameter `"__model__"`) per output and regression method.
#' @references Marino et al. (2008) *J. Theor. Biol.* 254:178 (PRCC); Saltelli et
#'   al. (2008) *Global Sensitivity Analysis: The Primer*.
#' @examples
#' set.seed(1)
#' X <- matrix(runif(1000 * 3), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
#' y <- 2 * X[, 1] - X[, 2]
#' gsa_regression(X, y, boot = 50)
#' @export
gsa_regression <- function(
  X,
  Y,
  methods = c("SRC", "PRCC"),
  parameters = NULL,
  boot = 100L,
  conf = 0.95,
  log_output = FALSE,
  dummy = NULL
) {
  methods <- toupper(methods)
  methods <- match.arg(methods, c("SRC", "SRRC", "PCC", "PRCC"), several.ok = TRUE)
  des <- .coerce_design(X)
  Xm <- des$X
  if (is.null(parameters)) {
    parameters <- des$parameters
  }
  Ym <- .coerce_outputs(Y)
  if (nrow(Xm) != nrow(Ym)) {
    .stop(c(
      "{.arg X} and {.arg Y} must have the same number of rows.",
      "x" = "{.arg X} has {nrow(Xm)} row{?s} but {.arg Y} has {nrow(Ym)}."
    ))
  }
  pnames <- colnames(Xm)
  onames <- colnames(Ym)
  boot <- as.integer(boot)

  Xsrc <- Xm
  if (!is.null(parameters) && inherits(parameters, "gsa_parameters")) {
    Xsrc <- .to_estimator_inputs(Xm, parameters)
  }

  one <- function(Xuse, y, idx, mth) {
    if (!is.null(idx)) {
      Xuse <- Xuse[idx, , drop = FALSE]
      y <- y[idx]
    }
    switch(
      mth,
      SRC = .src_one(Xuse, y, FALSE)$src,
      SRRC = .src_one(Xuse, y, TRUE)$src,
      PCC = .pcc_all(Xuse, y, FALSE),
      PRCC = .pcc_all(Xuse, y, TRUE)
    )
  }

  rows <- list()
  n_failed_total <- 0L
  for (oc in seq_along(onames)) {
    y0 <- Ym[, oc]
    keep <- is.finite(y0) & stats::complete.cases(Xm)
    n_failed_total <- n_failed_total + sum(!is.finite(y0))
    y <- y0[keep]
    Xk <- Xm[keep, , drop = FALSE]
    Xks <- Xsrc[keep, , drop = FALSE]
    if (log_output) {
      if (any(y <= 0)) {
        .stop(c(
          "{.arg log_output = TRUE} requires strictly positive outputs.",
          "x" = "Output {.val {onames[oc]}} contains non-positive values."
        ))
      }
      y <- log(y)
    }
    if (length(y) < ncol(Xk) + 2L || stats::sd(y) == 0) {
      .warn(c(
        "Skipping output {.val {onames[oc]}}.",
        "i" = "Too few valid runs ({length(y)}) or zero variance for a regression with {ncol(Xk)} predictor{?s}."
      ))
      next
    }
    Nk <- length(y)
    for (mth in methods) {
      Xuse <- if (mth %in% c("SRC", "SRRC")) Xks else Xk
      est <- one(Xuse, y, NULL, mth)
      bM <- matrix(NA_real_, boot, ncol(Xk))
      if (boot > 0L) {
        for (b in seq_len(boot)) {
          r <- sample.int(Nk, Nk, replace = TRUE)
          bM[b, ] <- tryCatch(one(Xuse, y, r, mth), error = function(e) rep(NA_real_, ncol(Xk)))
        }
      }
      for (i in seq_along(pnames)) {
        ci <- if (boot > 0L) .boot_ci(bM[, i], conf) else c(NA, NA)
        rows[[length(rows) + 1L]] <- .indices_dt(
          output = onames[oc],
          parameter = pnames[i],
          method = "regression",
          index = mth,
          estimate = est[i],
          std_error = if (boot > 0L) stats::sd(bM[, i], na.rm = TRUE) else NA_real_,
          conf_low = ci[[1]],
          conf_high = ci[[2]],
          conf_level = conf
        )
      }
      if (mth %in% c("SRC", "SRRC")) {
        r2 <- .src_one(Xuse, y, mth == "SRRC")$r2
        rows[[length(rows) + 1L]] <- .indices_dt(
          output = onames[oc],
          parameter = "__model__",
          method = "regression",
          index = paste0("R2_", mth),
          estimate = r2,
          conf_level = conf
        )
      }
    }
  }
  if (!length(rows)) {
    .stop(c(
      "No regression indices could be computed.",
      "i" = "All outputs were skipped (too few valid runs or zero variance)."
    ))
  }
  d <- .finalize_indices(data.table::rbindlist(rows, fill = TRUE), dummy = dummy)
  new_ospgsa_result(
    d,
    design = Xm,
    Y = Ym,
    meta = list(
      method = "regression",
      methods = methods,
      boot = boot,
      conf = conf,
      n_runs = nrow(Xm),
      n_failed = n_failed_total,
      log_output = log_output
    )
  )
}

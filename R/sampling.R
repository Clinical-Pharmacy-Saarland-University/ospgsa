.lhs <- function(n, k) {
  M <- matrix(0, n, k)
  for (j in seq_len(k)) {
    M[, j] <- (sample.int(n) - stats::runif(n)) / n
  }
  M
}

.map_marginals <- function(parameters, U) {
  k <- length(parameters)
  X <- matrix(0, nrow(U), k, dimnames = list(NULL, .param_names(parameters)))
  for (j in seq_len(k)) {
    X[, j] <- parameters[[j]]$qfun(U[, j])
  }
  X
}

# Iman & Conover (1982): induce target rank correlation, preserving marginals.
.iman_conover <- function(S, R) {
  n <- nrow(S)
  k <- ncol(S)
  if (k == 1L) {
    return(S)
  }
  scores <- stats::qnorm(seq_len(n) / (n + 1)) # van der Waerden scores
  M <- vapply(seq_len(k), function(j) sample(scores), numeric(n))
  cM <- stats::cor(M)
  Q <- chol(cM)
  Mstar <- M %*% solve(Q) # decorrelate: cor(Mstar) ~ I
  P <- chol(R)
  Tt <- Mstar %*% P # impose target: cor(Tt) ~ R
  out <- matrix(0, n, k, dimnames = dimnames(S))
  for (j in seq_len(k)) {
    out[, j] <- sort(S[, j])[rank(Tt[, j], ties.method = "first")]
  }
  out
}

#' Draw a correlated / uncorrelated input sample
#'
#' Generates a Monte Carlo design in natural parameter units, honouring each
#' marginal distribution and an optional correlation structure. This sample is
#' the basis for the moment-independent delta index and the regression methods.
#'
#' @param parameters A [gsa_parameters()] object.
#' @param n Sample size.
#' @param method `"lhs"` (Latin hypercube, default) or `"random"`.
#' @param correlation Optional correlation matrix from [gsa_correlation()]. If
#'   `NULL`, inputs are sampled independently.
#' @param correlation_method How a non-`NULL` correlation is induced:
#'   `"iman-conover"` (default; targets the **rank** correlation and preserves
#'   marginals/LHS exactly) or `"gaussian"` (Gaussian copula / NORTA via the
#'   normal scale).
#' @param corr_type Interpretation of `correlation`: `"spearman"` (rank, default)
#'   or `"pearson"`. For the Gaussian copula a Spearman target is mapped to the
#'   latent normal correlation via `rho_Z = 2 * sin(pi * rho_S / 6)`; Iman-Conover
#'   always induces a rank correlation.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `gsa_sample`: a list with `X` (an `n x k` numeric
#'   matrix in natural units, columns named by parameter), plus the originating
#'   `parameters`, `correlation` and method metadata.
#' @examples
#' p <- gsa_parameters(
#'   gsa_parameter("a", "uniform", min = 0, max = 1),
#'   gsa_parameter("b", "lognormal", median = 2, gsd = 1.4)
#' )
#' s <- gsa_sample(p, n = 200, correlation = gsa_correlation(p, c("a", "b", 0.6)))
#' stats::cor(s$X, method = "spearman")
#' @export
gsa_sample <- function(
  parameters,
  n,
  method = c("lhs", "random"),
  correlation = NULL,
  correlation_method = c("iman-conover", "gaussian"),
  corr_type = c("spearman", "pearson"),
  seed = NULL
) {
  .check_parameters(parameters)
  n <- .assert_count(n, min = 2L)
  method <- match.arg(method)
  correlation_method <- match.arg(correlation_method)
  corr_type <- match.arg(corr_type)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  k <- length(parameters)
  nms <- .param_names(parameters)

  if (is.null(correlation)) {
    U <- if (method == "lhs") .lhs(n, k) else matrix(stats::runif(n * k), n, k)
    X <- .map_marginals(parameters, U)
  } else {
    correlation <- correlation[nms, nms, drop = FALSE]
    if (correlation_method == "gaussian") {
      Rz <- correlation
      if (corr_type == "spearman") {
        Rz <- 2 * sin(pi * correlation / 6)
        diag(Rz) <- 1
        if (!.is_valid_corr(Rz)) Rz <- .nearest_pd_corr(Rz)
      }
      Z <- mvtnorm::rmvnorm(n, mean = rep(0, k), sigma = Rz)
      U <- stats::pnorm(Z)
      X <- .map_marginals(parameters, U)
    } else {
      U <- if (method == "lhs") .lhs(n, k) else matrix(stats::runif(n * k), n, k)
      S <- .map_marginals(parameters, U)
      X <- .iman_conover(S, correlation)
    }
  }
  structure(
    list(
      X = X,
      parameters = parameters,
      n = n,
      method = method,
      correlation = correlation,
      correlation_method = correlation_method,
      corr_type = corr_type
    ),
    class = "gsa_sample"
  )
}

#' @export
print.gsa_sample <- function(x, ...) {
  cat(sprintf(
    "<gsa_sample> n=%d, k=%d, method=%s%s\n",
    x$n,
    ncol(x$X),
    x$method,
    if (is.null(x$correlation)) {
      ""
    } else {
      sprintf(", correlated (%s)", x$correlation_method)
    }
  ))
  invisible(x)
}

#' Build the design pair for a two-stage delta analysis
#'
#' Draws the two designs that [gsa_delta_two_stage()] compares: an **independent**
#' (Stage 1) design from the marginals alone, and a **correlated** (Stage 2)
#' design from the joint distribution. The two use *distinct* seeds (`seed` for
#' the correlated design, `seed + 1` for the independent one) so their Monte-Carlo
#' errors are independent -- the convention required for the `delta_2 - delta_1`
#' contrast and the dummy noise floor. The returned designs are ready to feed to
#' [gsa_evaluate()] (e.g. with `checkpoint_dir =` for long runs).
#'
#' @param parameters A [gsa_parameters()] object.
#' @param n Design size per stage.
#' @param correlation A correlation matrix from [gsa_correlation()] for the
#'   correlated (Stage 2) design.
#' @param seed Optional integer seed for the correlated design; the independent
#'   design uses `seed + 1`.
#' @param method,correlation_method,corr_type Passed to [gsa_sample()].
#' @return A list with two [gsa_sample()] objects: `full` (correlated, Stage 2)
#'   and `ind` (independent, Stage 1).
#' @seealso [gsa_evaluate()], [gsa_delta_two_stage()]
#' @examples
#' p <- gsa_parameters(
#'   gsa_parameter("a", "uniform", min = 0, max = 1),
#'   gsa_parameter("b", "uniform", min = 0, max = 1)
#' )
#' R <- gsa_correlation(p, c("a", "b", 0.6))
#' d <- gsa_two_stage_designs(p, 200, correlation = R, seed = 1)
#' dim(d$full$X)
#' @export
gsa_two_stage_designs <- function(
  parameters,
  n,
  correlation,
  seed = NULL,
  method = c("lhs", "random"),
  correlation_method = c("iman-conover", "gaussian"),
  corr_type = c("spearman", "pearson")
) {
  .check_parameters(parameters)
  method <- match.arg(method)
  correlation_method <- match.arg(correlation_method)
  corr_type <- match.arg(corr_type)
  seed_ind <- if (is.null(seed)) NULL else as.integer(seed) + 1L
  full <- gsa_sample(
    parameters,
    n,
    method = method,
    correlation = correlation,
    correlation_method = correlation_method,
    corr_type = corr_type,
    seed = seed
  )
  ind <- gsa_sample(parameters, n, method = method, seed = seed_ind)
  list(full = full, ind = ind)
}

# Apply per-parameter scale (identity/log) for Sobol/Morris/regression; delta uses natural X.
.to_estimator_inputs <- function(X, parameters) {
  Z <- X
  for (j in seq_len(ncol(X))) {
    Z[, j] <- parameters[[j]]$tx(X[, j])
  }
  Z
}

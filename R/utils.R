`%||%` <- function(x, y) if (is.null(x)) y else x

.stop <- function(message, ..., .envir = parent.frame(), call = .envir) {
  cli::cli_abort(message, ..., .envir = .envir, call = call)
}

.warn <- function(message, ..., .envir = parent.frame()) {
  cli::cli_warn(message, ..., .envir = .envir)
}

.obj_desc <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  if (is.function(x)) {
    return("a function")
  }
  if (is.data.frame(x)) {
    return(sprintf("a data frame (%d x %d)", nrow(x), ncol(x)))
  }
  if (is.matrix(x)) {
    return(sprintf("a %s matrix (%d x %d)", typeof(x), nrow(x), ncol(x)))
  }
  if (is.atomic(x) && length(x) == 1L && !is.na(x)) {
    if (is.character(x)) {
      safe <- gsub("([{}])", "\\1\\1", x) # escape braces for glue/cli
      return(sprintf('the string "%s"', safe))
    }
    if (is.logical(x)) {
      return(if (isTRUE(x)) "TRUE" else "FALSE")
    }
    return(sprintf("the number %s", format(x)))
  }
  sprintf("a %s of length %d", class(x)[1L], length(x))
}

.assert_flag <- function(x, name = deparse(substitute(x)), call = parent.frame()) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    desc <- .obj_desc(x)
    .stop(
      c(
        "{.arg {name}} must be a single {.code TRUE} or {.code FALSE}.",
        "x" = "You supplied {desc}."
      ),
      call = call
    )
  }
  invisible(x)
}

.assert_count <- function(x, name = deparse(substitute(x)), min = 1L, call = parent.frame()) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != round(x)) {
    desc <- .obj_desc(x)
    .stop(
      c("{.arg {name}} must be a single whole number >= {min}.", "x" = "You supplied {desc}."),
      call = call
    )
  }
  invisible(as.integer(x))
}

.assert_number <- function(x, name = deparse(substitute(x)), call = parent.frame()) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    desc <- .obj_desc(x)
    .stop(
      c("{.arg {name}} must be a single finite number.", "x" = "You supplied {desc}."),
      call = call
    )
  }
  invisible(as.numeric(x))
}

.assert_string <- function(x, name = deparse(substitute(x)), call = parent.frame()) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    desc <- .obj_desc(x)
    .stop(c("{.arg {name}} must be a single string.", "x" = "You supplied {desc}."), call = call)
  }
  invisible(x)
}

.check_parameters <- function(x, arg = "parameters", call = parent.frame()) {
  if (!inherits(x, "gsa_parameters")) {
    desc <- .obj_desc(x)
    .stop(
      c(
        "{.arg {arg}} must be a {.cls gsa_parameters} object.",
        "i" = "Build one with {.fn gsa_parameters}.",
        "x" = "You supplied {desc}."
      ),
      call = call
    )
  }
  invisible(x)
}

.check_evaluator <- function(x, arg = "evaluator", call = parent.frame()) {
  if (!is.function(x)) {
    desc <- .obj_desc(x)
    .stop(
      c(
        "{.arg {arg}} must be a model evaluator (a function).",
        "i" = "Build one with {.fn ospsuite_evaluator} (OSP model) or {.fn function_evaluator}.",
        "x" = "You supplied {desc}."
      ),
      call = call
    )
  }
  invisible(x)
}

.abort_eval_rows <- function(got, expected, call = parent.frame()) {
  .stop(
    c(
      "The evaluator returned the wrong number of rows.",
      "x" = "Expected {expected} row{?s} (one per design point) but got {got}.",
      "i" = "An evaluator must return one output row per input row, with {.code NA} for failed runs."
    ),
    call = call
  )
}

# x assumed sorted increasing
.trapz <- function(x, y) {
  n <- length(x)
  if (n < 2L) {
    return(0)
  }
  sum((x[-1L] - x[-n]) * (y[-1L] + y[-n])) / 2
}

.boot_ci <- function(reps, conf = 0.95) {
  reps <- reps[is.finite(reps)]
  if (length(reps) < 2L) {
    return(c(low = NA_real_, high = NA_real_))
  }
  a <- (1 - conf) / 2
  q <- stats::quantile(reps, probs = c(a, 1 - a), names = FALSE, na.rm = TRUE)
  c(low = q[1L], high = q[2L])
}

.is_valid_corr <- function(R, tol = 1e-8) {
  if (!is.matrix(R) || nrow(R) != ncol(R)) {
    return(FALSE)
  }
  if (any(!is.finite(R))) {
    return(FALSE)
  }
  if (max(abs(R - t(R))) > 1e-6) {
    return(FALSE)
  }
  if (any(abs(diag(R) - 1) > 1e-6)) {
    return(FALSE)
  }
  ev <- tryCatch(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), error = function(e) {
    -Inf
  })
  ev > -tol
}

# Higham nearest positive-definite correlation matrix
.nearest_pd_corr <- function(R, eig_tol = 1e-8, maxit = 100L) {
  R <- (R + t(R)) / 2
  X <- R
  DS <- matrix(0, nrow(R), ncol(R))
  for (i in seq_len(maxit)) {
    Rk <- X - DS
    e <- eigen(Rk, symmetric = TRUE)
    V <- e$vectors
    d <- pmax(e$values, eig_tol)
    Xnew <- V %*% diag(d, length(d)) %*% t(V)
    DS <- Xnew - Rk
    X <- Xnew
    diag(X) <- 1
    if (max(abs(Xnew - X)) < eig_tol) break
  }
  X <- (X + t(X)) / 2
  diag(X) <- 1
  stats::cov2cor(X)
}

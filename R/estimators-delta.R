# Estimator: Plischke, Borgonovo and Smith (2013), Eur J Oper Res 226(3):536-550.
.coerce_design <- function(X, call = parent.frame()) {
  if (inherits(X, "gsa_sample")) {
    return(list(X = X$X, parameters = X$parameters))
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X) || !is.numeric(X)) {
    desc <- .obj_desc(X)
    .stop(
      c(
        "{.arg X} must be a numeric matrix, data frame or {.cls gsa_sample}.",
        "x" = "You supplied {desc}."
      ),
      call = call
    )
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  list(X = X, parameters = NULL)
}

.coerce_outputs <- function(Y) {
  if (is.list(Y) && !is.data.frame(Y) && !is.null(Y$Y)) {
    Y <- Y$Y
  }
  if (is.numeric(Y) && is.null(dim(Y))) {
    Y <- matrix(Y, ncol = 1L, dimnames = list(NULL, "Y"))
  } else if (is.data.frame(Y)) {
    nm <- names(Y)
    Y <- as.matrix(Y)
    colnames(Y) <- nm
  }
  if (!is.matrix(Y)) {
    desc <- .obj_desc(Y)
    .stop(c(
      "{.arg Y} must be a numeric vector, matrix, data frame or evaluation result.",
      "x" = "You supplied {desc}."
    ))
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("Y", seq_len(ncol(Y)))
  }
  Y
}

# Number of equal-frequency rank classes (SALib calibration formula).
.delta_n_classes <- function(N, min_per_class = 10L, cap = 48L) {
  ex <- 2 / (7 + tanh((1500 - N) / 500))
  M <- min(ceiling(N^ex), cap)
  M <- min(M, max(1L, floor(N / min_per_class)))
  as.integer(max(M, 1L))
}

.kde_eval <- function(y, grid) {
  n <- length(y)
  if (n < 2L) {
    return(rep(0, length(grid)))
  }
  s <- stats::sd(y)
  if (!is.finite(s) || s == 0) {
    bw <- diff(range(grid)) / length(grid)
    return(stats::dnorm(grid, mean = y[1L], sd = max(bw, 1e-9)))
  }
  d <- stats::density(y, bw = "nrd0", n = 512L, from = grid[1L], to = grid[length(grid)])
  stats::approx(d$x, d$y, xout = grid, rule = 2)$y
}

.calc_delta_one <- function(Xi, Y, grid, M, min_per_class = 10L) {
  N <- length(Y)
  fY <- .kde_eval(Y, grid)
  rX <- rank(Xi, ties.method = "first")
  bounds <- seq(0, N, length.out = M + 1L)
  acc <- 0
  for (j in seq_len(M)) {
    idx <- which(rX > bounds[j] & rX <= bounds[j + 1L])
    nm <- length(idx)
    if (nm < min_per_class) {
      next
    }
    fYc <- .kde_eval(Y[idx], grid)
    acc <- acc + (nm / (2 * N)) * .trapz(grid, abs(fY - fYc))
  }
  acc
}

.calc_s1_partition <- function(Xi, Y, M) {
  N <- length(Y)
  mY <- mean(Y)
  vY <- stats::var(Y)
  if (!is.finite(vY) || vY == 0) {
    return(0)
  }
  rX <- rank(Xi, ties.method = "first")
  bounds <- seq(0, N, length.out = M + 1L)
  num <- 0
  for (j in seq_len(M)) {
    idx <- which(rX > bounds[j] & rX <= bounds[j + 1L])
    nm <- length(idx)
    if (nm < 1L) {
      next
    }
    num <- num + (nm / N) * (mean(Y[idx]) - mY)^2
  }
  num / vY
}

# Run expr under rng_seed's L'Ecuyer-CMRG stream, restoring caller's RNG after; NULL = ambient stream.
.with_rng_stream <- function(rng_seed, expr) {
  if (is.null(rng_seed)) {
    return(force(expr))
  }
  old_kind <- RNGkind()
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(
    {
      RNGkind(old_kind[1], old_kind[2], old_kind[3])
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", rng_seed, envir = .GlobalEnv)
  force(expr)
}

.delta_param <- function(Xi, Y, grid, M, boot, conf, ci, min_per_class, rng_seed = NULL) {
  dhat <- .calc_delta_one(Xi, Y, grid, M, min_per_class)
  if (boot > 0L) {
    N <- length(Y)
    db <- .with_rng_stream(rng_seed, {
      out <- numeric(boot)
      for (b in seq_len(boot)) {
        r <- sample.int(N, N, replace = TRUE)
        out[b] <- .calc_delta_one(Xi[r], Y[r], grid, M, min_per_class)
      }
      out
    })
    dc <- 2 * dhat - db
    est <- mean(dc)
    bias <- mean(db) - dhat
    se <- stats::sd(db)
    if (ci == "percentile") {
      q <- .boot_ci(dc, conf)
      lo <- q[[1]]
      hi <- q[[2]]
    } else {
      z <- stats::qnorm(0.5 + conf / 2)
      lo <- est - z * se
      hi <- est + z * se
    }
  } else {
    est <- dhat
    bias <- NA_real_
    se <- NA_real_
    lo <- NA_real_
    hi <- NA_real_
  }
  list(estimate = est, bias = bias, std_error = se, conf_low = lo, conf_high = hi)
}

# Build n independent L'Ecuyer-CMRG substreams without disturbing caller's RNG; NULL seed draws+advances the ambient stream.
.delta_make_seeds <- function(n, seed = NULL) {
  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  root <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as.integer(seed)
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  RNGkind("L'Ecuyer-CMRG")
  set.seed(root)
  s <- get(".Random.seed", envir = .GlobalEnv)
  seeds <- vector("list", n)
  for (i in seq_len(n)) {
    seeds[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  if (is.null(saved)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", saved, envir = .GlobalEnv)
  }
  seeds
}

#' Borgonovo delta moment-independent sensitivity index
#'
#' Estimates the delta index from a single input-output sample using the
#' given-data estimator (Plischke, Borgonovo and Smith, 2013): equal-frequency
#' rank classes, Gaussian KDE of the unconditional and class-conditional output
#' densities, half-L1 shift weighted by class probability, bootstrap bias
#' correction and bootstrap confidence intervals.
#'
#' The delta index is **moment-independent** (it responds to changes in any
#' feature of the output distribution, not only its variance) and is **well
#' defined under correlated inputs**, which makes it the recommended index for
#' PBPK models. See [gsa_delta_two_stage()] to separate a parameter's structural
#' effect from the part it inherits through correlation.
#'
#' @param X Input design: a [gsa_sample()], a numeric matrix, or a data.frame
#'   (one column per parameter, in natural units). Column names are the parameter
#'   labels.
#' @param Y Model output(s): a numeric vector (single output), a matrix /
#'   data.frame (one column per output), or an evaluation result from
#'   [gsa_evaluate()]. Failed runs encoded as `NA` are dropped per output.
#' @param boot Number of bootstrap resamples for bias correction and CIs (0 to
#'   skip). Use >= 500 for publication.
#' @param conf Confidence level for the bootstrap interval.
#' @param ci Bootstrap interval type: `"percentile"` (default, recommended for
#'   this bounded, near-zero-skewed index) or `"normal"`.
#' @param min_per_class Minimum number of points required in a rank class for its
#'   conditional KDE to be used.
#' @param grid_n Number of output grid points for the KDE integration.
#' @param include_s1 Also report the first-order Sobol index obtained for free
#'   from the same rank partition (a cross-check; no CI).
#' @param log_output Analyse `log(Y)` instead of `Y` (recommended for strongly
#'   right-skewed metrics such as AUC and Cmax). Requires positive outputs.
#' @param classes Optional fixed number of rank classes (overrides the automatic
#'   rule).
#' @param dummy Optional name of a phantom parameter used as the significance
#'   noise floor (see [gsa_parameters()] and [gsa_dummy()]).
#' @param seed Optional integer. When supplied (or whenever `n_cores > 1`), each
#'   (output, parameter) task's bootstrap is run on an independent
#'   L'Ecuyer-CMRG substream derived from `seed`, so results are reproducible
#'   **and identical whether computed serially or in parallel**. When `NULL` and
#'   `n_cores == 1`, the ambient RNG stream is used (unchanged legacy behaviour).
#' @param n_cores Number of workers for evaluating the (output x parameter)
#'   tasks in parallel. `1` (default) runs serially. Parallelism uses a PSOCK
#'   cluster (works on Windows; no `fork`).
#' @param cl Optional pre-built `parallel` cluster to reuse (its workers must
#'   have ospgsa available, e.g. an installed package or one loaded via
#'   `pkgload::load_all` on each worker). Overrides `n_cores`; not stopped by
#'   this function. Useful to share one cluster across both delta stages.
#'
#' @return An [ospgsa_result].
#' @references
#' Borgonovo E (2007) _Reliab Eng Syst Saf_ 92(6):771-784. \doi{10.1016/j.ress.2006.04.015}.
#' Plischke E, Borgonovo E, Smith CL (2013) _Eur J Oper Res_ 226(3):536-550. \doi{10.1016/j.ejor.2012.11.047}.
#' @seealso [gsa_delta_two_stage()], [gsa_sobol()]
#' @examples
#' set.seed(1)
#' X <- matrix(runif(2000 * 3), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
#' y <- 3 * X[, 1] + X[, 2]                 # c is inert
#' gsa_delta(X, y, boot = 50)
#' @export
gsa_delta <- function(
  X,
  Y,
  boot = 100L,
  conf = 0.95,
  ci = c("percentile", "normal"),
  min_per_class = 10L,
  grid_n = 110L,
  include_s1 = TRUE,
  log_output = FALSE,
  classes = NULL,
  dummy = NULL,
  seed = NULL,
  n_cores = 1L,
  cl = NULL
) {
  ci <- match.arg(ci)
  des <- .coerce_design(X)
  Xm <- des$X
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

  prep <- vector("list", 0L)
  n_failed_total <- 0L
  for (oc in seq_along(onames)) {
    y0 <- Ym[, oc]
    keep <- is.finite(y0) & stats::complete.cases(Xm)
    n_failed_total <- n_failed_total + sum(!is.finite(y0))
    y <- y0[keep]
    Xk <- Xm[keep, , drop = FALSE]
    if (log_output) {
      if (any(y <= 0)) {
        .stop(c(
          "{.arg log_output = TRUE} requires strictly positive outputs.",
          "x" = "Output {.val {onames[oc]}} contains non-positive values."
        ))
      }
      y <- log(y)
    }
    N <- length(y)
    if (N < 2L * min_per_class || stats::sd(y) == 0) {
      .warn(c(
        "Skipping output {.val {onames[oc]}}.",
        "i" = "It has too few valid runs ({N}) or zero variance for the delta estimator."
      ))
      next
    }
    M <- if (is.null(classes)) .delta_n_classes(N, min_per_class) else as.integer(classes)
    rng <- range(y)
    pad <- 0.05 * diff(rng)
    if (pad == 0) {
      pad <- 1e-6
    }
    grid <- seq(rng[1] - pad, rng[2] + pad, length.out = grid_n)
    prep[[length(prep) + 1L]] <- list(output = onames[oc], y = y, Xk = Xk, M = M, grid = grid)
  }
  if (!length(prep)) {
    .stop(c(
      "No delta indices could be computed.",
      "i" = "All outputs were skipped (too few valid runs or zero variance). Increase {.arg n} or check the model."
    ))
  }

  tasks <- list()
  for (pi in seq_along(prep)) {
    for (j in seq_along(pnames)) {
      tasks[[length(tasks) + 1L]] <- c(pi, j)
    }
  }
  n_tasks <- length(tasks)

  resolved_cores <- if (!is.null(cl)) length(cl) else max(1L, as.integer(n_cores))
  use_streams <- resolved_cores > 1L || !is.null(seed)
  task_seeds <- if (use_streams) .delta_make_seeds(n_tasks, seed) else vector("list", n_tasks)

  run_task <- function(t, rng_seed) {
    pp <- prep[[t[1]]]
    xi <- pp$Xk[, t[2]]
    r <- .delta_param(xi, pp$y, pp$grid, pp$M, boot, conf, ci, min_per_class, rng_seed = rng_seed)
    s1 <- if (include_s1) .calc_s1_partition(xi, pp$y, pp$M) else NA_real_
    list(r = r, s1 = s1)
  }

  if (resolved_cores > 1L) {
    if (is.null(cl)) {
      cl <- parallel::makePSOCKcluster(resolved_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      # For a SOURCE package pass a pre-initialised `cl` (workers pkgload::load_all'd).
      try(
        parallel::clusterCall(cl, function() {
          suppressWarnings(requireNamespace("ospgsa", quietly = TRUE))
        }),
        silent = TRUE
      )
    }
    ns <- parent.env(environment())
    parallel::clusterExport(
      cl,
      c("prep", "boot", "conf", "ci", "min_per_class", "include_s1"),
      envir = environment()
    )
    parallel::clusterExport(
      cl,
      c(
        ".delta_param",
        ".with_rng_stream",
        ".calc_delta_one",
        ".kde_eval",
        ".trapz",
        ".calc_s1_partition",
        ".boot_ci"
      ),
      envir = ns
    )
    worker <- run_task
    environment(worker) <- globalenv()
    res <- parallel::clusterMap(
      cl,
      worker,
      tasks,
      task_seeds,
      SIMPLIFY = FALSE,
      .scheduling = "dynamic"
    )
  } else {
    res <- Map(run_task, tasks, task_seeds)
  }

  rows <- vector("list", 0L)
  for (ti in seq_along(tasks)) {
    out <- prep[[tasks[[ti]][1]]]$output
    pn <- pnames[tasks[[ti]][2]]
    rr <- res[[ti]]$r
    rows[[length(rows) + 1L]] <- .indices_dt(
      output = out,
      parameter = pn,
      method = "delta",
      index = "delta",
      stage = NA_character_,
      estimate = rr$estimate,
      bias = rr$bias,
      std_error = rr$std_error,
      conf_low = rr$conf_low,
      conf_high = rr$conf_high,
      conf_level = conf
    )
    if (include_s1) {
      rows[[length(rows) + 1L]] <- .indices_dt(
        output = out,
        parameter = pn,
        method = "delta",
        index = "S1",
        stage = NA_character_,
        estimate = res[[ti]]$s1,
        conf_level = conf
      )
    }
  }
  d <- data.table::rbindlist(rows, fill = TRUE)
  d <- .finalize_indices(d, dummy = dummy)
  new_ospgsa_result(
    d,
    design = Xm,
    Y = Ym,
    meta = list(
      method = "delta",
      boot = boot,
      conf = conf,
      ci = ci,
      log_output = log_output,
      n_runs = nrow(Xm),
      n_failed = n_failed_total,
      n_cores = resolved_cores
    )
  )
}

#' Two-stage delta analysis for correlated inputs
#'
#' Computes the delta index on two designs and compares them to separate each
#' parameter's **structural / causal** effect from the effect it inherits through
#' **correlation** with other parameters (De Carlo et al. 2023; adapted to PBPK
#' by Cuquerella-Gilabert et al. 2026):
#'
#' * **Stage 1 (independent):** parameters sampled from their marginals with all
#'   correlations switched off -> `delta_1`. A non-zero `delta_1` certifies a
#'   genuine direct effect (delta = 0 iff Y is independent of the parameter *and*
#'   the parameter is uncorrelated with the others).
#' * **Stage 2 (full):** parameters sampled from the realistic correlated joint
#'   distribution -> `delta_2` (importance under the true output distribution).
#'
#' Each parameter is then classified as *causal*, *indirect only* (influential
#' purely via correlation) or *both* (see [delta_classification()]).
#'
#' @param X_ind,Y_ind Stage-1 independent design and its outputs.
#' @param X_full,Y_full Stage-2 correlated design and its outputs.
#' @param boot,conf,ci,min_per_class,grid_n,log_output,classes,dummy Passed to
#'   [gsa_delta()].
#' @param seed,n_cores,cl Bootstrap parallelism, passed to [gsa_delta()] for each
#'   stage. The two stages use distinct (offset) seeds so their bootstraps are
#'   independent; results are reproducible and identical to the serial path. A
#'   single cluster is built (or the supplied `cl` reused) across both stages.
#'
#' @return An [ospgsa_result] whose `indices` contains `stage` values
#'   `"independent"` (delta_1), `"full"` (delta_2) and `"correlative"`
#'   (delta_2 - delta_1, an attribution heuristic; delta is not additive).
#' @seealso [delta_classification()], [gsa()]
#' @examples
#' set.seed(1)
#' p <- gsa_parameters(
#'   gsa_parameter("k2",   "uniform", min = 0, max = 1),
#'   gsa_parameter("IC50", "uniform", min = 0, max = 1)  # inert but correlated to k2
#' )
#' R  <- gsa_correlation(p, c("k2", "IC50", 0.9))
#' Sf <- gsa_sample(p, 2000, correlation = R, seed = 1)
#' Si <- gsa_sample(p, 2000, seed = 1)
#' f  <- function(M) 2 * M[, "k2"]                       # only k2 acts
#' gsa_delta_two_stage(Sf$X, f(Sf$X), Si$X, f(Si$X), boot = 50)
#' @export
gsa_delta_two_stage <- function(
  X_full,
  Y_full,
  X_ind,
  Y_ind,
  boot = 100L,
  conf = 0.95,
  ci = c("percentile", "normal"),
  min_per_class = 10L,
  grid_n = 110L,
  log_output = FALSE,
  classes = NULL,
  dummy = NULL,
  seed = NULL,
  n_cores = 1L,
  cl = NULL
) {
  ci <- match.arg(ci)
  if (is.null(cl) && as.integer(n_cores) > 1L) {
    cl <- parallel::makePSOCKcluster(as.integer(n_cores))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    try(
      parallel::clusterCall(cl, function() {
        suppressWarnings(requireNamespace("ospgsa", quietly = TRUE))
      }),
      silent = TRUE
    )
  }
  # Distinct seeds per stage: shared resamples would correlate the delta_2 - delta_1 attribution.
  seed_ind <- if (is.null(seed)) NULL else as.integer(seed) + 1L
  seed_full <- if (is.null(seed)) NULL else as.integer(seed) + 2L
  r_ind <- gsa_delta(
    X_ind,
    Y_ind,
    boot = boot,
    conf = conf,
    ci = ci,
    min_per_class = min_per_class,
    grid_n = grid_n,
    include_s1 = FALSE,
    log_output = log_output,
    classes = classes,
    dummy = dummy,
    seed = seed_ind,
    n_cores = n_cores,
    cl = cl
  )
  r_full <- gsa_delta(
    X_full,
    Y_full,
    boot = boot,
    conf = conf,
    ci = ci,
    min_per_class = min_per_class,
    grid_n = grid_n,
    include_s1 = FALSE,
    log_output = log_output,
    classes = classes,
    dummy = dummy,
    seed = seed_full,
    n_cores = n_cores,
    cl = cl
  )
  di <- r_ind$indices[index == "delta"][, stage := "independent"]
  df <- r_full$indices[index == "delta"][, stage := "full"]

  m <- merge(
    di[, .(output, parameter, est_ind = estimate)],
    df[, .(output, parameter, est_full = estimate)],
    by = c("output", "parameter")
  )
  dc <- m[, .indices_dt(
    output = output,
    parameter = parameter,
    method = "delta",
    index = "delta",
    stage = "correlative",
    estimate = est_full - est_ind,
    conf_level = conf
  )]
  d <- data.table::rbindlist(list(di, df, dc), fill = TRUE)
  d <- .finalize_indices(d, dummy = dummy)
  new_ospgsa_result(
    d,
    meta = list(
      method = "delta_two_stage",
      methods = "delta_two_stage",
      boot = boot,
      conf = conf,
      ci = ci,
      log_output = log_output,
      correlated = TRUE,
      dummy = dummy,
      n_runs = r_full$meta$n_runs,
      n_failed = r_full$meta$n_failed
    )
  )
}

#' Classify parameters from a two-stage delta analysis
#'
#' @param result An [ospgsa_result] from [gsa_delta_two_stage()] (or [gsa()] with
#'   `method = "delta"` and a correlation matrix).
#' @return A `data.table` with one row per output / parameter giving `delta_1`,
#'   `delta_2`, their significance, and a `class` in
#'   `{causal, indirect-only, both, non-influential}`.
#' @export
delta_classification <- function(result) {
  d <- result$indices
  ind <- d[
    index == "delta" & stage == "independent",
    .(output, parameter, d1 = estimate, sig1 = significant)
  ]
  full <- d[
    index == "delta" & stage == "full",
    .(output, parameter, d2 = estimate, sig2 = significant)
  ]
  if (!nrow(ind) || !nrow(full)) {
    .stop(c(
      "This result does not contain two-stage delta indices.",
      "i" = "Produce one with {.fn gsa_delta_two_stage}, or {.fn gsa} with {.code method = \"delta\"} and a {.arg correlation} matrix."
    ))
  }
  m <- merge(ind, full, by = c("output", "parameter"))
  m[,
    class := data.table::fcase(
      sig1 & sig2  , "both"          ,
      !sig1 & sig2 , "indirect-only" ,
      sig1 & !sig2 , "causal"        ,
      default = "non-influential"
    )
  ]
  data.table::setorder(m, output, -d2)
  m[]
}

#' Create a parallel cluster with ospgsa loaded on the workers
#'
#' Convenience wrapper for parallel GSA: builds a PSOCK cluster (works on
#' Windows; no `fork`) and **silently** loads ospgsa on every worker, so it can
#' be passed as `cl` to [gsa_delta()], [gsa_delta_two_stage()] or [gsa()]. Use it
#' instead of a hand-rolled `parallel::clusterCall()`, which prints its
#' per-worker return value at the top level of a script.
#'
#' @param n_cores Number of workers.
#' @param pkg_dir Optional path to the ospgsa **source** tree; when given, workers
#'   load it with `pkgload::load_all()` (for a not-installed / development
#'   package). When `NULL`, workers run `library(ospgsa)` (installed package).
#' @param libpaths Library paths to prepend on each worker (e.g. a non-default
#'   user library). Defaults to the current `.libPaths()`.
#' @return A running PSOCK cluster. Stop it with [parallel::stopCluster()] when
#'   done (e.g. via `on.exit()`).
#' @examples
#' \dontrun{
#' cl <- ospgsa_cluster(8, pkg_dir = "path/to/ospgsa")
#' on.exit(parallel::stopCluster(cl))
#' res <- gsa_delta_two_stage(Xf, Yf, Xi, Yi, seed = 1, cl = cl)
#' }
#' @export
ospgsa_cluster <- function(n_cores, pkg_dir = NULL, libpaths = .libPaths()) {
  n_cores <- max(1L, as.integer(n_cores))
  cl <- parallel::makePSOCKcluster(n_cores)
  tryCatch(
    # invisible() so the per-worker return value is not auto-printed.
    invisible(parallel::clusterCall(
      cl,
      function(lp, dir) {
        if (length(lp)) {
          .libPaths(lp)
        }
        loaded <- FALSE
        if (!is.null(dir) && requireNamespace("pkgload", quietly = TRUE)) {
          loaded <- tryCatch(
            {
              suppressWarnings(suppressMessages(pkgload::load_all(dir, quiet = TRUE)))
              TRUE
            },
            error = function(e) FALSE
          )
        }
        if (!loaded) {
          suppressWarnings(suppressMessages(requireNamespace("ospgsa", quietly = TRUE)))
        }
        invisible(TRUE)
      },
      libpaths,
      pkg_dir
    )),
    error = function(e) {
      parallel::stopCluster(cl)
      stop(e)
    }
  )
  cl
}

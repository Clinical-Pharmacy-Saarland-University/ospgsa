.dist_choices <- c("uniform", "loguniform", "normal", "lognormal", "truncnorm")

.make_qfun <- function(dist, p) {
  switch(
    dist,
    uniform = function(u) stats::qunif(u, min = p$min, max = p$max),
    loguniform = function(u) exp(stats::qunif(u, min = log(p$min), max = log(p$max))),
    normal = function(u) stats::qnorm(u, mean = p$mean, sd = p$sd),
    lognormal = function(u) stats::qlnorm(u, meanlog = p$meanlog, sdlog = p$sdlog),
    truncnorm = function(u) {
      lo <- stats::pnorm(p$lower, p$mean, p$sd)
      hi <- stats::pnorm(p$upper, p$mean, p$sd)
      stats::qnorm(lo + u * (hi - lo), mean = p$mean, sd = p$sd)
    }
  )
}

.resolve_dist <- function(dist, args, call = parent.frame()) {
  need <- function(...) {
    miss <- setdiff(c(...), names(args))
    if (length(miss)) {
      .stop(
        c(
          "Distribution {.val {dist}} is missing required argument{?s}: {.arg {miss}}.",
          "i" = "See {.fn gsa_parameter} for the arguments each distribution needs."
        ),
        call = call
      )
    }
  }
  p <- list()
  if (dist == "uniform") {
    need("min", "max")
    p$min <- args$min
    p$max <- args$max
    if (p$min >= p$max) {
      .stop(
        c(
          "A {.val uniform} distribution needs {.arg min} < {.arg max}.",
          "x" = "You supplied {.arg min} = {.val {p$min}} and {.arg max} = {.val {p$max}}."
        ),
        call = call
      )
    }
  } else if (dist == "loguniform") {
    need("min", "max")
    p$min <- args$min
    p$max <- args$max
    if (p$min <= 0) {
      .stop(
        c(
          "A {.val loguniform} distribution needs {.arg min} > 0.",
          "x" = "You supplied {.arg min} = {.val {p$min}}."
        ),
        call = call
      )
    }
    if (p$min >= p$max) {
      .stop(
        c(
          "A {.val loguniform} distribution needs {.arg min} < {.arg max}.",
          "x" = "You supplied {.arg min} = {.val {p$min}} and {.arg max} = {.val {p$max}}."
        ),
        call = call
      )
    }
  } else if (dist == "normal") {
    need("mean", "sd")
    p$mean <- args$mean
    p$sd <- args$sd
    if (p$sd <= 0) {
      .stop(
        c(
          "A {.val normal} distribution needs {.arg sd} > 0.",
          "x" = "You supplied {.arg sd} = {.val {p$sd}}."
        ),
        call = call
      )
    }
  } else if (dist == "lognormal") {
    if (all(c("meanlog", "sdlog") %in% names(args))) {
      p$meanlog <- args$meanlog
      p$sdlog <- args$sdlog
    } else if (all(c("median", "gsd") %in% names(args))) {
      if (args$median <= 0 || args$gsd <= 1) {
        .stop(
          c(
            "A {.val lognormal} via {.arg median}/{.arg gsd} needs {.arg median} > 0 and {.arg gsd} > 1.",
            "x" = "You supplied {.arg median} = {.val {args$median}} and {.arg gsd} = {.val {args$gsd}}."
          ),
          call = call
        )
      }
      p$meanlog <- log(args$median)
      p$sdlog <- log(args$gsd)
    } else if (all(c("mean", "cv") %in% names(args))) {
      if (args$mean <= 0 || args$cv <= 0) {
        .stop(
          c(
            "A {.val lognormal} via {.arg mean}/{.arg cv} needs {.arg mean} > 0 and {.arg cv} > 0.",
            "x" = "You supplied {.arg mean} = {.val {args$mean}} and {.arg cv} = {.val {args$cv}}."
          ),
          call = call
        )
      }
      p$sdlog <- sqrt(log(1 + args$cv^2))
      p$meanlog <- log(args$mean) - p$sdlog^2 / 2
    } else {
      .stop(
        c(
          "A {.val lognormal} distribution needs one of these argument pairs:",
          "*" = "{.arg meanlog} and {.arg sdlog}",
          "*" = "{.arg median} and {.arg gsd} (geometric SD)",
          "*" = "{.arg mean} and {.arg cv} (mean and coefficient of variation)"
        ),
        call = call
      )
    }
    if (p$sdlog <= 0) {
      .stop(
        c(
          "A {.val lognormal} distribution needs {.arg sdlog} > 0.",
          "x" = "The supplied parameters give {.arg sdlog} = {.val {p$sdlog}}."
        ),
        call = call
      )
    }
  } else if (dist == "truncnorm") {
    need("mean", "sd", "lower", "upper")
    p$mean <- args$mean
    p$sd <- args$sd
    p$lower <- args$lower
    p$upper <- args$upper
    if (p$sd <= 0) {
      .stop(
        c(
          "A {.val truncnorm} distribution needs {.arg sd} > 0.",
          "x" = "You supplied {.arg sd} = {.val {p$sd}}."
        ),
        call = call
      )
    }
    if (p$lower >= p$upper) {
      .stop(
        c(
          "A {.val truncnorm} distribution needs {.arg lower} < {.arg upper}.",
          "x" = "You supplied {.arg lower} = {.val {p$lower}} and {.arg upper} = {.val {p$upper}}."
        ),
        call = call
      )
    }
  }
  p
}

#' Describe one uncertain model parameter
#'
#' Defines the marginal uncertainty distribution of a single model parameter for
#' global sensitivity analysis.
#'
#' @param name Short, unique label used in tables and plots.
#' @param dist Distribution: one of `"uniform"`, `"loguniform"`, `"normal"`,
#'   `"lognormal"`, `"truncnorm"`.
#' @param ... Distribution arguments (in **natural** parameter units):
#'   * `uniform`, `loguniform`: `min`, `max` (loguniform needs `min > 0`).
#'   * `normal`: `mean`, `sd`.
#'   * `lognormal`: one of `(meanlog, sdlog)`, `(median, gsd)` (geometric SD), or
#'     `(mean, cv)` (mean and coefficient of variation on the natural scale).
#'   * `truncnorm`: `mean`, `sd`, `lower`, `upper`.
#' @param path PK-Sim parameter path used by [ospsuite_evaluator()]. Defaults to
#'   `name`. Ignored by non-PK-Sim evaluators.
#' @param unit Optional display unit of the PK-Sim parameter. If supplied, sampled
#'   values are converted to base units before being set (see [ospsuite_evaluator()]).
#'   Leave `NULL` to provide values already in base units.
#' @param scale Representation used for **scale-dependent** estimators (Sobol,
#'   Morris, regression): `"natural"` or `"log"`. Defaults to `"log"` for
#'   `loguniform`/`lognormal` parameters and `"natural"` otherwise. The
#'   moment-independent delta index is invariant to this choice. See the GSA
#'   primer in the package repository for the log- vs linear-scale discussion.
#' @param mol_weight Optional molecular weight (used only for unit conversions of
#'   amount/concentration parameters).
#' @param mol_weight_unit Unit of `mol_weight` (default `"g/mol"`).
#'
#' @return An object of class `gsa_parameter`.
#' @seealso [gsa_parameters()], [gsa_correlation()]
#' @examples
#' gsa_parameter("Lipophilicity", "normal", mean = -0.1, sd = 0.3,
#'               path = "Aciclovir|Lipophilicity")
#' gsa_parameter("CL", "lognormal", median = 5, gsd = 1.5)
#' @export
gsa_parameter <- function(
  name,
  dist,
  ...,
  path = name,
  unit = NULL,
  scale = NULL,
  mol_weight = NULL,
  mol_weight_unit = "g/mol"
) {
  .assert_string(name)
  dist <- match.arg(dist, .dist_choices)
  args <- list(...)
  p <- .resolve_dist(dist, args, call = environment())
  if (is.null(scale)) {
    scale <- if (dist %in% c("loguniform", "lognormal")) "log" else "natural"
  }
  scale <- match.arg(scale, c("natural", "log"))

  qfun <- .make_qfun(dist, p)
  support <- c(qfun(1e-6), qfun(1 - 1e-6))

  tx <- switch(scale, natural = function(x) x, log = function(x) log(x))
  itx <- switch(scale, natural = function(x) x, log = function(x) exp(x))
  if (scale == "log" && any(support <= 0)) {
    .stop(c(
      "{.arg scale} = {.val log} requires strictly positive support.",
      "x" = "Parameter {.val {name}} can take non-positive values under its {.val {dist}} distribution.",
      "i" = "Use a positive distribution (e.g. {.val lognormal}/{.val loguniform}) or set {.code scale = \"natural\"}."
    ))
  }

  structure(
    list(
      name = name,
      dist = dist,
      dist_args = p,
      path = path,
      unit = unit,
      scale = scale,
      qfun = qfun,
      tx = tx,
      itx = itx,
      support = support,
      nominal = qfun(0.5),
      mol_weight = mol_weight,
      mol_weight_unit = mol_weight_unit,
      is_dummy = FALSE
    ),
    class = "gsa_parameter"
  )
}

#' Add a phantom ("dummy") parameter as a significance noise floor
#'
#' A dummy parameter is sampled like any other input but is **not connected to
#' the model** (PK-Sim evaluators skip it). Its estimated sensitivity reflects
#' pure Monte Carlo / estimator noise, giving an empirical significance floor:
#' any real parameter whose index is indistinguishable from the dummy's is not
#' meaningfully influential. Pass the dummy's `name` as the `dummy` argument to
#' the estimators / [gsa()].
#'
#' @param name Label for the dummy parameter (default `"dummy"`).
#' @return A `gsa_parameter` flagged as a dummy (uniform on `[0, 1]`).
#' @examples
#' gsa_parameters(
#'   gsa_parameter("a", "uniform", min = 0, max = 1),
#'   gsa_dummy()
#' )
#' @export
gsa_dummy <- function(name = "dummy") {
  p <- gsa_parameter(name, "uniform", min = 0, max = 1, path = NA_character_)
  p$is_dummy <- TRUE
  p
}

#' @export
print.gsa_parameter <- function(x, ...) {
  da <- paste(
    sprintf("%s=%s", names(x$dist_args), formatC(unlist(x$dist_args), digits = 4)),
    collapse = ", "
  )
  cat(sprintf(
    "<gsa_parameter> %s  [%s: %s]  scale=%s%s\n",
    x$name,
    x$dist,
    da,
    x$scale,
    if (!is.null(x$unit)) sprintf("  unit=%s", x$unit) else ""
  ))
  cat(sprintf("  path: %s\n", x$path))
  invisible(x)
}

#' Collect parameter descriptions
#'
#' @param ... One or more [gsa_parameter()] objects, or a single list of them.
#' @return An object of class `gsa_parameters` (an ordered, named collection).
#' @seealso [gsa_parameter()], [gsa_correlation()]
#' @examples
#' gsa_parameters(
#'   gsa_parameter("a", "uniform", min = 0, max = 1),
#'   gsa_parameter("b", "lognormal", median = 2, gsd = 1.5)
#' )
#' @export
gsa_parameters <- function(...) {
  items <- list(...)
  if (length(items) == 1L && is.list(items[[1]]) && !inherits(items[[1]], "gsa_parameter")) {
    items <- items[[1]]
  }
  if (!length(items)) {
    .stop(c(
      "At least one parameter is required.",
      "i" = "Pass {.fn gsa_parameter} objects, e.g. {.code gsa_parameters(gsa_parameter(\"CL\", \"lognormal\", median = 5, gsd = 1.5))}."
    ))
  }
  ok <- vapply(items, inherits, logical(1), "gsa_parameter")
  if (!all(ok)) {
    bad <- which(!ok)
    .stop(c(
      "Every argument to {.fn gsa_parameters} must be a {.cls gsa_parameter}.",
      "x" = "Invalid at argument position(s): {bad}.",
      "i" = "Create them with {.fn gsa_parameter} / {.fn gsa_dummy}."
    ))
  }
  nms <- vapply(items, function(p) p$name, character(1))
  if (anyDuplicated(nms)) {
    dup <- unique(nms[duplicated(nms)])
    .stop(c("Parameter names must be unique.", "x" = "Duplicated name{?s}: {.val {dup}}."))
  }
  names(items) <- nms
  structure(items, class = "gsa_parameters")
}

#' @export
print.gsa_parameters <- function(x, ...) {
  cat(sprintf("<gsa_parameters> %d parameter(s)\n", length(x)))
  for (p in x) {
    print(p)
  }
  invisible(x)
}

#' @export
length.gsa_parameters <- function(x) length(unclass(x))

.param_names <- function(params) vapply(params, function(p) p$name, character(1))
.param_paths <- function(params) vapply(params, function(p) p$path, character(1))
.param_scales <- function(params) vapply(params, function(p) p$scale, character(1))
.param_nominal <- function(params) vapply(params, function(p) p$nominal, numeric(1))
.param_is_dummy <- function(params) vapply(params, function(p) isTRUE(p$is_dummy), logical(1))
.dummy_name <- function(params) {
  d <- .param_is_dummy(params)
  if (any(d)) unname(.param_names(params)[which(d)[1L]]) else NULL
}

# `width` per dist: gsd (lognormal), CV-fraction of nominal (truncnorm), absolute SD (normal).
.gsa_anchor_parameter <- function(
  name,
  dist,
  nominal,
  width,
  lower = NA_real_,
  upper = NA_real_,
  path = name,
  unit = NULL,
  call = parent.frame()
) {
  nominal <- as.numeric(nominal)
  width <- as.numeric(width)
  switch(
    dist,
    lognormal = gsa_parameter(
      name,
      "lognormal",
      median = nominal,
      gsd = width,
      path = path,
      unit = unit
    ),
    truncnorm = gsa_parameter(
      name,
      "truncnorm",
      mean = nominal,
      sd = max(width * nominal, 1e-6),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      path = path,
      unit = unit
    ),
    normal = gsa_parameter(name, "normal", mean = nominal, sd = width, path = path, unit = unit),
    .stop(
      c(
        "Unsupported {.arg dist} {.val {dist}} for an anchored parameter.",
        "i" = "Use {.val lognormal}, {.val truncnorm} or {.val normal}, or build it directly with {.fn gsa_parameter}."
      ),
      call = call
    )
  )
}

#' Build / validate an input correlation matrix
#'
#' Constructs a (rank) correlation matrix over the parameters. Supply pairwise
#' correlations as named triplets, or pass a full matrix to validate it. The
#' correlation is interpreted as a Spearman/rank target and induced with the
#' Iman-Conover method, or as a Gaussian-copula correlation, depending on the
#' sampler (see [gsa_sample()]).
#'
#' @param parameters A [gsa_parameters()] object (defines size and ordering).
#' @param ... Pairwise correlations as `c("a", "b", 0.6)` triplets (character
#'   names, numeric value), a single **list** of such triplets, or a single full
#'   numeric matrix.
#' @param repair If `TRUE` (default) a non-positive-definite matrix is projected
#'   to the nearest valid correlation matrix (with a warning).
#' @return A named correlation matrix.
#' @examples
#' p <- gsa_parameters(
#'   gsa_parameter("a", "uniform", min = 0, max = 1),
#'   gsa_parameter("b", "uniform", min = 0, max = 1),
#'   gsa_parameter("c", "uniform", min = 0, max = 1)
#' )
#' gsa_correlation(p, c("a", "b", 0.7), c("a", "c", -0.3))
#' # equivalently, pass the triplets as one list:
#' gsa_correlation(p, list(c("a", "b", 0.7), c("a", "c", -0.3)))
#' @export
gsa_correlation <- function(parameters, ..., repair = TRUE) {
  .check_parameters(parameters)
  nms <- .param_names(parameters)
  k <- length(nms)
  spec <- list(...)

  if (
    length(spec) == 1L &&
      is.list(spec[[1]]) &&
      !is.data.frame(spec[[1]]) &&
      !is.matrix(spec[[1]]) &&
      all(vapply(spec[[1]], function(e) length(e) == 3L, logical(1)))
  ) {
    spec <- spec[[1]]
  }

  if (length(spec) == 1L && is.matrix(spec[[1]])) {
    R <- spec[[1]]
    if (is.null(dimnames(R)) || is.null(rownames(R))) {
      if (nrow(R) != k) {
        .stop(c(
          "The correlation matrix has the wrong dimension.",
          "x" = "It is {nrow(R)}x{ncol(R)} but there are {k} parameter{?s}.",
          "i" = "Provide a {k}x{k} matrix, or name its rows/columns after the parameters."
        ))
      }
      dimnames(R) <- list(nms, nms)
    }
    R <- R[nms, nms, drop = FALSE]
  } else {
    R <- diag(k)
    dimnames(R) <- list(nms, nms)
    for (tr in spec) {
      if (length(tr) != 3L) {
        .stop(c(
          "Each pairwise correlation must be a length-3 triplet.",
          "i" = "Use {.code c(\"name1\", \"name2\", value)}, e.g. {.code c(\"CL\", \"V\", 0.6)}.",
          "x" = "Got a value of length {length(tr)}."
        ))
      }
      a <- as.character(tr[[1]])
      b <- as.character(tr[[2]])
      v <- as.numeric(tr[[3]])
      unknown <- setdiff(c(a, b), nms)
      if (length(unknown)) {
        .stop(c(
          "Unknown parameter{?s} in a correlation triplet: {.val {unknown}}.",
          "i" = "Known parameters are {.val {nms}}."
        ))
      }
      if (is.na(v) || abs(v) > 1) {
        .stop(c(
          "A correlation must be a number in {.field [-1, 1]}.",
          "x" = "You supplied {.val {v}} for {.val {a}}--{.val {b}}."
        ))
      }
      R[a, b] <- v
      R[b, a] <- v
    }
  }
  diag(R) <- 1
  if (!.is_valid_corr(R)) {
    if (!repair) {
      .stop(c(
        "The correlation matrix is not positive semi-definite.",
        "i" = "Set {.code repair = TRUE} to project it to the nearest valid matrix."
      ))
    }
    .warn(c(
      "The correlation matrix is not positive (semi-)definite.",
      "i" = "Projecting it to the nearest valid correlation matrix."
    ))
    R <- .nearest_pd_corr(R)
    dimnames(R) <- list(nms, nms)
  }
  R
}

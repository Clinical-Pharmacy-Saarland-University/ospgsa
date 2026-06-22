.need_ospsuite <- function(call = parent.frame()) {
  if (!requireNamespace("ospsuite", quietly = TRUE)) {
    .stop(
      c(
        "The {.pkg ospsuite} package is required for the PK-Sim evaluator.",
        "i" = "Install it from {.url https://www.open-systems-pharmacology.org/}."
      ),
      call = call
    )
  }
}

.short_output <- function(path) {
  seg <- strsplit(path, "|", fixed = TRUE)[[1L]]
  utils::tail(seg, 1L)
}

.extract_outputs <- function(sr, out_paths, out_labels, pk_params, time_points, col_names) {
  vals <- stats::setNames(rep(NA_real_, length(col_names)), col_names)
  if (is.null(sr)) {
    return(vals)
  }

  if (length(pk_params)) {
    pk <- tryCatch(
      {
        df <- ospsuite::pkAnalysesToDataFrame(ospsuite::calculatePKAnalyses(results = sr))
        data.table::as.data.table(df)
      },
      error = function(e) NULL
    )
    if (!is.null(pk) && nrow(pk)) {
      for (oi in seq_along(out_paths)) {
        for (pp in pk_params) {
          v <- pk[QuantityPath == out_paths[oi] & Parameter == pp, Value]
          cn <- paste0(out_labels[oi], "__", pp)
          if (length(v)) vals[cn] <- suppressWarnings(as.numeric(v[1L]))
        }
      }
    }
  }

  if (length(time_points)) {
    for (oi in seq_along(out_paths)) {
      ov <- tryCatch(
        ospsuite::getOutputValues(sr, quantitiesOrPaths = out_paths[oi]),
        error = function(e) NULL
      )
      if (!is.null(ov) && !is.null(ov$data)) {
        d <- ov$data
        tt <- d[["Time"]]
        yy <- d[[out_paths[oi]]]
        if (!is.null(tt) && !is.null(yy)) {
          yi <- stats::approx(tt, yy, xout = time_points, rule = 2)$y
          for (ti in seq_along(time_points)) {
            vals[paste0(out_labels[oi], "@t", ti)] <- yi[ti]
          }
        }
      }
    }
  }
  vals
}

# One simulated run as a named list of data.frame(time, value), one per output.
.sr_profiles <- function(sr, out_paths, out_labels) {
  stats::setNames(
    lapply(seq_along(out_paths), function(oi) {
      ov <- tryCatch(
        ospsuite::getOutputValues(sr, quantitiesOrPaths = out_paths[oi]),
        error = function(e) NULL
      )
      if (is.null(ov) || is.null(ov$data)) {
        return(data.frame(time = numeric(0), value = numeric(0)))
      }
      d <- ov$data
      data.frame(time = d[["Time"]], value = d[[out_paths[oi]]])
    }),
    out_labels
  )
}

# Align a summary_fn's (named) return onto the fixed output columns, NA-filling
# missing names and dropping extras. Unnamed returns map positionally if lengths
# match. Keeps the evaluator's "matrix in -> fixed matrix out" contract.
.align_summary <- function(v, col_names) {
  vals <- stats::setNames(rep(NA_real_, length(col_names)), col_names)
  if (is.null(v)) {
    return(vals)
  }
  v <- unlist(v)
  nm <- names(v)
  v <- suppressWarnings(as.numeric(v))
  if (is.null(nm)) {
    if (length(v) == length(col_names)) {
      vals[] <- v
    }
    return(vals)
  }
  names(v) <- nm
  keep <- intersect(nm, col_names)
  vals[keep] <- v[keep]
  vals
}

#' Build a PK-Sim (ospsuite) model evaluator
#'
#' Returns an evaluator that, given a matrix of parameter values, runs the
#' PK-Sim simulation once per row and returns scalar PK parameters and/or
#' time-resolved concentrations. The result plugs directly into [gsa()],
#' [gsa_sobol()], [gsa_morris()] and [gsa_evaluate()].
#'
#' The evaluator loads **one** compiled model and a **single** `SimulationBatch`,
#' and reuses them for every call: all rows are queued on that batch and solved
#' concurrently across `n_cores` (a single batch parallelises its queued runs
#' internally). Memory is therefore ~one model copy regardless of `n_cores`, and
#' there is no per-call reloading.
#'
#' Values for parameters declared with a display `unit` are converted to base
#' units (with `ospsuite::toBaseUnit()`) before each run. Dummy parameters
#' ([gsa_dummy()]) are ignored by the model. Failed runs are returned as `NA`
#' rows, which the estimators drop.
#'
#' @param simulation A loaded `ospsuite` `Simulation`, or a path to a `.pkml`
#'   file.
#' @param parameters A [gsa_parameters()] object. Each non-dummy parameter's
#'   `path` must resolve in the simulation.
#' @param outputs Output quantity path(s) to record. If `NULL`, the simulation's
#'   existing output selection(s) are used.
#' @param pk_parameters Scalar PK parameters to extract (subset of
#'   `ospsuite::allPKParameterNames()`), e.g. `c("AUC_tEnd", "C_max", "t_max")`.
#'   Set to `character(0)` for time-resolved only.
#' @param time_points Optional numeric vector of times (in the simulation's base
#'   time unit, minutes) at which to record concentrations (time-resolved GSA).
#' @param n_cores Number of cores used to solve the queued runs of the single
#'   batch concurrently. Scaling plateaus near the number of physical cores
#'   (hyperthreads add little for ODE solving).
#' @param output_names Optional named character vector mapping output paths to
#'   short labels used in result column names.
#' @param silent Suppress the ospsuite run progress/among messages.
#' @param chunk Maximum number of rows queued and solved per inner sub-batch.
#'   Bounds the number of `SimulationResults` held at once for very large single
#'   calls; does not change which model is used (always the one warm model).
#' @param summary_fn Optional custom reducer, a function `function(profiles, sr)`
#'   that turns one simulated run into a named numeric vector of outputs. When
#'   supplied it REPLACES the built-in `pk_parameters` / `time_points` extraction,
#'   so you can compute any metric (`t_max`, a partial AUC, a slope, a ratio, time
#'   above a threshold, ...). `profiles` is a named list (one per output) of
#'   `data.frame(time, value)` in base units. `sr` is the raw `SimulationResults`
#'   (use `ospsuite::calculatePKAnalyses(sr)` for PK-Sim's own PK parameters). The
#'   returned names become the result columns and must be the same on every run.
#'   Develop and debug it with [ospsuite_test_summary()].
#'
#' @return A function with the evaluator contract (matrix in -> matrix out).
#' @seealso [gsa()], [function_evaluator()], [ospsuite_test_summary()]
#' @examples
#' \dontrun{
#' sim <- ospsuite::loadSimulation(
#'   system.file("extdata", "Aciclovir.pkml", package = "ospsuite"))
#' p <- gsa_parameters(
#'   gsa_parameter("Lipophilicity", "normal", mean = -0.1, sd = 0.3,
#'                 path = "Aciclovir|Lipophilicity"))
#' ev <- ospsuite_evaluator(sim, p, pk_parameters = c("AUC_tEnd", "C_max"))
#' ev(gsa_sample(p, 8)$X)
#' }
#' @export
ospsuite_evaluator <- function(
  simulation,
  parameters,
  outputs = NULL,
  pk_parameters = c("AUC_tEnd", "C_max", "t_max"),
  time_points = NULL,
  n_cores = 1L,
  output_names = NULL,
  silent = TRUE,
  chunk = 2000L,
  summary_fn = NULL
) {
  .need_ospsuite()
  .check_parameters(parameters)
  sim_path <- NULL
  if (is.character(simulation)) {
    sim_path <- simulation
    if (!file.exists(sim_path)) {
      .stop(c("Simulation file not found.", "x" = "No file at {.path {sim_path}}."))
    }
  } else if (inherits(simulation, "Simulation")) {
    sim_path <- simulation$sourceFile
    if (is.null(sim_path) || !nzchar(sim_path) || !file.exists(sim_path)) {
      .stop(c(
        "Cannot locate the simulation's source {.file .pkml} file.",
        "i" = "Pass the {.file .pkml} path directly so independent copies can be loaded for batching."
      ))
    }
  } else {
    desc <- .obj_desc(simulation)
    .stop(c(
      "{.arg simulation} must be an {.pkg ospsuite} {.cls Simulation} or a path to a {.file .pkml} file.",
      "x" = "You supplied {desc}."
    ))
  }
  n_cores <- .assert_count(n_cores, min = 1L)
  chunk <- .assert_count(chunk, min = 1L)
  pk_parameters <- as.character(pk_parameters)
  use_summary <- !is.null(summary_fn)
  if (use_summary && !is.function(summary_fn)) {
    desc <- .obj_desc(summary_fn)
    .stop(c("{.arg summary_fn} must be a function.", "x" = "You supplied {desc}."))
  }
  if (!use_summary && !length(pk_parameters) && !length(time_points)) {
    .stop(c(
      "Nothing to record from the model.",
      "i" = "Provide {.arg pk_parameters} (e.g. {.val AUC_tEnd}), {.arg time_points}, or a {.arg summary_fn}."
    ))
  }

  is_dum <- .param_is_dummy(parameters)
  var_idx <- which(!is_dum)
  if (!length(var_idx)) {
    .stop(c(
      "There are no model parameters to vary.",
      "x" = "All supplied parameters are dummies ({.fn gsa_dummy})."
    ))
  }
  var_names <- .param_names(parameters)[var_idx]
  var_paths <- .param_paths(parameters)[var_idx]
  var_units <- lapply(parameters[var_idx], function(p) p$unit)
  var_mw <- lapply(parameters[var_idx], function(p) p$mol_weight)
  var_mwu <- lapply(parameters[var_idx], function(p) p$mol_weight_unit)

  sim0 <- ospsuite::loadSimulation(sim_path, loadFromCache = FALSE, addToCache = FALSE)
  if (!is.null(outputs)) {
    ospsuite::clearOutputs(sim0)
    ospsuite::addOutputs(quantitiesOrPaths = outputs, simulation = sim0)
    out_paths <- as.character(outputs)
  } else {
    out_paths <- NULL
  }
  ref_params <- lapply(var_paths, function(pp) ospsuite::getParameter(pp, sim0))

  to_base_vec <- function(disp_vals) {
    vapply(
      seq_along(disp_vals),
      function(j) {
        u <- var_units[[j]]
        if (is.null(u)) {
          return(disp_vals[j])
        }
        ospsuite::toBaseUnit(
          ref_params[[j]],
          disp_vals[j],
          u,
          molWeight = var_mw[[j]],
          molWeightUnit = var_mwu[[j]]
        )
      },
      numeric(1)
    )
  }

  b0 <- ospsuite::createSimulationBatch(simulation = sim0, parametersOrPaths = var_paths)
  nom <- to_base_vec(.param_nominal(parameters)[var_idx])
  b0$addRunValues(parameterValues = nom)
  res0 <- ospsuite::runSimulationBatches(simulationBatches = list(b0), silentMode = silent)
  sr0 <- res0[[1L]][[1L]]
  if (is.null(sr0)) {
    .stop(c(
      "The nominal PK-Sim run failed.",
      "i" = "Check that the simulation solves at the parameters' median values and that the parameter paths are correct.",
      "i" = "Parameter path{?s}: {.val {var_paths}}."
    ))
  }
  if (is.null(out_paths)) {
    out_paths <- sr0$allQuantityPaths
  }
  out_labels <- if (!is.null(output_names)) {
    unname(ifelse(
      out_paths %in% names(output_names),
      output_names[out_paths],
      vapply(out_paths, .short_output, "")
    ))
  } else {
    vapply(out_paths, .short_output, "")
  }
  if (use_summary) {
    probe <- tryCatch(
      unlist(summary_fn(.sr_profiles(sr0, out_paths, out_labels), sr0)),
      error = function(e) {
        .stop(c(
          "{.arg summary_fn} failed on the nominal run.",
          "x" = conditionMessage(e),
          "i" = "Debug it with {.fn ospsuite_test_summary}."
        ))
      }
    )
    if (!length(probe) || is.null(names(probe)) || any(!nzchar(names(probe)))) {
      .stop(c(
        "{.arg summary_fn} must return a named numeric vector.",
        "i" = "Give every metric a name, e.g. {.code c(tmax = ..., partial_auc = ...)}."
      ))
    }
    col_names <- names(probe)
  } else {
    col_names <- character(0)
    if (length(pk_parameters)) {
      col_names <- c(col_names, as.vector(t(outer(out_labels, pk_parameters, paste, sep = "__"))))
    }
    if (length(time_points)) {
      col_names <- c(
        col_names,
        as.vector(t(outer(out_labels, paste0("t", seq_along(time_points)), paste, sep = "@")))
      )
    }
  }

  # sim0/b0 are reused across every call; b0's queue clears after each run.
  run_opts <- ospsuite::SimulationRunOptions$new(numberOfCores = n_cores, showProgress = !silent)
  function(X) {
    X <- as.matrix(X)
    if (is.null(colnames(X))) {
      colnames(X) <- .param_names(parameters)
    }
    miss <- setdiff(var_names, colnames(X))
    if (length(miss)) {
      .stop(c(
        "The design matrix is missing parameter column{?s}: {.val {miss}}.",
        "i" = "Its columns must be named after the parameters; {.fn gsa_sample} does this automatically."
      ))
    }
    n <- nrow(X)
    Y <- matrix(NA_real_, n, length(col_names), dimnames = list(NULL, col_names))
    for (s in seq.int(1L, n, by = chunk)) {
      idx <- s:min(s + chunk - 1L, n)
      ids <- character(length(idx))
      for (k in seq_along(idx)) {
        ids[k] <- b0$addRunValues(parameterValues = to_base_vec(as.numeric(X[idx[k], var_names])))
      }
      res <- ospsuite::runSimulationBatches(
        simulationBatches = list(b0),
        simulationRunOptions = run_opts,
        silentMode = silent,
        stopIfFails = FALSE
      )
      bres <- res[[b0$id]]
      for (k in seq_along(idx)) {
        sr <- tryCatch(bres[[ids[k]]], error = function(e) NULL)
        Y[idx[k], ] <- if (is.null(sr)) {
          rep(NA_real_, length(col_names))
        } else if (use_summary) {
          .align_summary(
            tryCatch(summary_fn(.sr_profiles(sr, out_paths, out_labels), sr), error = function(e) {
              NULL
            }),
            col_names
          )
        } else {
          .extract_outputs(sr, out_paths, out_labels, pk_parameters, time_points, col_names)
        }
      }
      rm(res, bres)
    }
    Y
  }
}

#' Test and debug a custom `summary_fn`
#'
#' Runs the simulation once (at the parameters' nominal values, or at `at`) and
#' shows exactly what an [ospsuite_evaluator()] `summary_fn` receives and returns,
#' so you can develop a custom output reducer without launching a full GSA. Unlike
#' inside a GSA run, the `summary_fn` is called here WITHOUT error trapping, so a
#' mistake surfaces directly (with its traceback).
#'
#' @param simulation A loaded `ospsuite` `Simulation`, or a path to a `.pkml`.
#' @param parameters A [gsa_parameters()] object.
#' @param summary_fn The reducer to test: `function(profiles, sr)` returning a
#'   named numeric vector (see [ospsuite_evaluator()]).
#' @param outputs Output quantity path(s) to record, as in [ospsuite_evaluator()].
#' @param at Optional named numeric vector of parameter values to simulate at;
#'   parameters not named fall back to their nominal value. Default: all nominal.
#' @param n_cores Cores for the single run (default `1`).
#' @return Invisibly, a list with `value` (the `summary_fn` return), `profiles`
#'   (the per-output `data.frame(time, value)` list it was given) and `sr` (the
#'   raw `SimulationResults`). Prints a readable diagnostic.
#' @seealso [ospsuite_evaluator()]
#' @examples
#' \dontrun{
#' sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
#' p <- gsa_parameters(
#'   gsa_parameter("Lipophilicity", "normal", mean = -0.1, sd = 0.3,
#'                 path = "Aciclovir|Lipophilicity"))
#' tmax <- function(profiles, sr) {
#'   d <- profiles[[1]]
#'   c(tmax = d$time[which.max(d$value)])
#' }
#' ospsuite_test_summary(sim, p, tmax)
#' }
#' @export
ospsuite_test_summary <- function(
  simulation,
  parameters,
  summary_fn,
  outputs = NULL,
  at = NULL,
  n_cores = 1L
) {
  .need_ospsuite()
  .check_parameters(parameters)
  if (!is.function(summary_fn)) {
    desc <- .obj_desc(summary_fn)
    .stop(c("{.arg summary_fn} must be a function.", "x" = "You supplied {desc}."))
  }
  capture <- new.env(parent = emptyenv())
  cap_fn <- function(profiles, sr) {
    capture$profiles <- profiles
    capture$sr <- sr
    c(.probe = 0)
  }
  ev <- ospsuite_evaluator(
    simulation,
    parameters,
    outputs = outputs,
    summary_fn = cap_fn,
    n_cores = n_cores,
    silent = TRUE
  )
  pnames <- .param_names(parameters)
  x <- stats::setNames(.param_nominal(parameters), pnames)
  if (!is.null(at)) {
    common <- intersect(names(at), pnames)
    x[common] <- as.numeric(at[common])
  }
  probe_out <- ev(matrix(x, nrow = 1L, dimnames = list(NULL, pnames)))
  if (all(is.na(probe_out))) {
    .stop(c(
      "The simulation run failed at the requested parameter values.",
      "i" = "Check that the model solves there before debugging {.arg summary_fn}."
    ))
  }
  prof <- capture$profiles
  sr <- capture$sr

  cli::cli_h2("Profiles passed to summary_fn (per output)")
  for (nm in names(prof)) {
    p <- prof[[nm]]
    if (nrow(p)) {
      cli::cli_inform(
        "{.field {nm}}: {nrow(p)} point{?s}, t {signif(min(p$time), 3)} to {signif(max(p$time), 3)}, value {signif(min(p$value), 3)} to {signif(max(p$value), 3)}"
      )
    } else {
      cli::cli_alert_warning("{.field {nm}}: no data returned for this output.")
    }
  }

  out <- summary_fn(prof, sr)
  v <- unlist(out)
  cli::cli_h2("summary_fn returned")
  if (is.null(names(v)) || any(!nzchar(names(v)))) {
    cli::cli_alert_warning("Unnamed value(s): a GSA matches columns BY NAME, so name every metric.")
  }
  if (!is.numeric(v)) {
    cli::cli_alert_warning("Not numeric: summary_fn must return numbers.")
  }
  if (any(!is.finite(suppressWarnings(as.numeric(v))))) {
    cli::cli_alert_warning("Some values are non-finite (NA/NaN/Inf) at this point.")
  }
  print(out)
  invisible(list(value = out, profiles = prof, sr = sr))
}

.as_simulation <- function(simulation, call = parent.frame()) {
  .need_ospsuite(call)
  if (inherits(simulation, "Simulation")) {
    return(simulation)
  }
  if (is.character(simulation) && length(simulation) == 1L) {
    if (!file.exists(simulation)) {
      .stop(c("Simulation file not found.", "x" = "No file at {.path {simulation}}."), call = call)
    }
    return(ospsuite::loadSimulation(simulation, loadFromCache = FALSE, addToCache = FALSE))
  }
  desc <- .obj_desc(simulation)
  .stop(
    c(
      "{.arg simulation} must be an {.pkg ospsuite} {.cls Simulation} or a path to a {.file .pkml} file.",
      "x" = "You supplied {desc}."
    ),
    call = call
  )
}

.path_near_matches <- function(pattern, paths, n = 5L) {
  toks <- unlist(strsplit(gsub("\\\\", "", pattern), "[^A-Za-z0-9]+"))
  toks <- unique(tolower(toks[nchar(toks) >= 3L]))
  if (!length(toks)) {
    return(character(0))
  }
  lp <- tolower(paths)
  score <- vapply(
    lp,
    function(p) sum(vapply(toks, grepl, logical(1), x = p, fixed = TRUE)),
    integer(1)
  )
  hit <- which(score > 0L)
  if (!length(hit)) {
    return(character(0))
  }
  utils::head(paths[hit][order(score[hit], decreasing = TRUE)], n)
}

.resolve_unique <- function(name, pattern, paths, call = parent.frame()) {
  hits <- grep(pattern, paths, value = TRUE, perl = TRUE)
  if (length(hits) == 1L) {
    return(hits)
  }
  if (length(hits) == 0L) {
    near <- .path_near_matches(pattern, paths)
    .stop(
      c(
        "Parameter {.val {name}}: pattern {.val {pattern}} matched no parameter path.",
        if (length(near)) c("i" = "Did you mean: {.val {near}}") else NULL,
        "i" = "List or search the model's paths with {.fn ospsuite_parameter_paths}."
      ),
      call = call
    )
  }
  .stop(
    c(
      "Parameter {.val {name}}: pattern {.val {pattern}} matched {length(hits)} paths (need exactly 1).",
      "i" = "Make the pattern more specific.",
      "x" = "Matches: {.val {hits}}."
    ),
    call = call
  )
}

#' Resolve parameter-path patterns against an OSP simulation
#'
#' Matches regular-expression patterns against the absolute parameter paths of a
#' PK-Sim/MoBi simulation, requiring each to resolve to **exactly one** path
#' (erroring with the candidate list otherwise). Handy for turning short, stable
#' regexes into the full paths that [gsa_parameter()] / [ospsuite_evaluator()]
#' need.
#'
#' @param simulation A loaded `ospsuite` `Simulation`, or a path to a `.pkml`.
#' @param patterns A character vector of PCRE patterns. If named, the names are
#'   used in error messages and carried onto the result; otherwise each pattern
#'   labels itself.
#' @return A named character vector mapping each pattern's name to its unique
#'   absolute parameter path.
#' @seealso [ospsuite_parameters()]
#' @examples
#' \dontrun{
#' sim <- ospsuite::loadSimulation(
#'   system.file("extdata", "Aciclovir.pkml", package = "ospsuite"))
#' ospsuite_resolve_paths(sim, c(lipo = "Aciclovir\\|Lipophilicity$"))
#' }
#' @export
ospsuite_resolve_paths <- function(simulation, patterns) {
  sim <- .as_simulation(simulation)
  patterns <- as.character(patterns)
  nms <- names(patterns)
  if (is.null(nms)) {
    nms <- patterns
  }
  nms[!nzchar(nms)] <- patterns[!nzchar(nms)]
  all_paths <- ospsuite::getAllParameterPathsIn(sim)
  out <- vapply(
    seq_along(patterns),
    function(i) {
      .resolve_unique(nms[i], patterns[i], all_paths)
    },
    character(1)
  )
  stats::setNames(out, nms)
}

#' List or search the parameter paths of an OSP simulation
#'
#' Discovery aid for finding a parameter's path in a PK-Sim/MoBi model. Returns
#' the absolute parameter paths, optionally filtered by a pattern and optionally
#' with each path's current value. Use it to find the path (or build the regex)
#' that [ospsuite_parameters()] / [ospsuite_resolve_paths()] need, instead of
#' reaching for raw `ospsuite::getAllParameterPathsIn()` plus `grep()`.
#'
#' @param simulation A loaded `ospsuite` `Simulation`, or a path to a `.pkml`.
#' @param pattern Optional search string. `NULL` (default) returns every path.
#' @param ignore_case Match case-insensitively (default `TRUE`).
#' @param fixed If `TRUE`, treat `pattern` as a literal substring (no regex), so
#'   you do not have to escape the `|` path separator. Default `FALSE` (PCRE).
#' @param value If `TRUE`, return a data frame of `path` and current `value`
#'   (base units) rather than a character vector.
#' @param max Maximum number of paths to return (default `Inf`).
#' @return A sorted character vector of matching paths, or, when `value = TRUE`,
#'   a data frame with columns `path` and `value`.
#' @seealso [ospsuite_resolve_paths()], [ospsuite_parameters()]
#' @examples
#' \dontrun{
#' sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
#' ospsuite_parameter_paths(sim, "lipophilicity")
#' ospsuite_parameter_paths(sim, "Fraction unbound", fixed = TRUE, value = TRUE)
#' }
#' @export
ospsuite_parameter_paths <- function(
  simulation,
  pattern = NULL,
  ignore_case = TRUE,
  fixed = FALSE,
  value = FALSE,
  max = Inf
) {
  sim <- .as_simulation(simulation)
  paths <- sort(ospsuite::getAllParameterPathsIn(sim))
  if (!is.null(pattern)) {
    .assert_string(pattern)
    if (fixed && ignore_case) {
      paths <- paths[grepl(tolower(pattern), tolower(paths), fixed = TRUE)]
    } else {
      paths <- grep(
        pattern,
        paths,
        value = TRUE,
        perl = !fixed,
        fixed = fixed,
        ignore.case = ignore_case
      )
    }
  }
  if (is.finite(max) && length(paths) > max) {
    paths <- utils::head(paths, max)
  }
  if (!isTRUE(value)) {
    return(paths)
  }
  vals <- vapply(
    paths,
    function(p) {
      as.numeric(tryCatch(ospsuite::getParameter(p, sim)$value, error = function(e) NA_real_))
    },
    numeric(1)
  )
  data.frame(path = paths, value = unname(vals), stringsAsFactors = FALSE)
}

#' Build model-anchored GSA parameters from a specification table
#'
#' One-call setup of a [gsa_parameters()] object for a PK-Sim/MoBi model: it
#' resolves each parameter's path from a regex pattern, reads the model's current
#' value at that path (in base units) and **anchors** the marginal there, and
#' (optionally) appends a dummy parameter. This replaces the hand-rolled
#' resolve/anchor/validate boilerplate of a typical driver script.
#'
#' The marginal is anchored at the nominal value `m` per `dist`:
#' * `lognormal` -> `median = m`, `gsd = width` (geometric SD),
#' * `truncnorm` -> `mean = m`, `sd = width * m` (so `width` is a CV), with
#'   `lower`/`upper` from the table,
#' * `normal` -> `mean = m`, `sd = width` (absolute).
#'
#' Values are taken in base units (`unit = NULL`), matching the values read from
#' the model, so no unit conversion is applied. Lognormal anchors must be
#' strictly positive (checked).
#'
#' @param simulation A loaded `ospsuite` `Simulation`, or a path to a `.pkml`.
#' @param specs A data frame (or list coercible to one) with columns `name`,
#'   `pattern`, `dist`, `width`, and optional `lower`, `upper` (used by
#'   `truncnorm`).
#' @param dummy Name for an appended [gsa_dummy()] noise-floor parameter, or
#'   `NULL`/`NA` to add none. Default `"dummy"`.
#' @return A [gsa_parameters()] object, with the resolved table (name, dist,
#'   nominal, path) attached as attribute `"resolved"`.
#' @seealso [ospsuite_resolve_paths()], [gsa_parameter()], [gsa_correlation()]
#' @examples
#' \dontrun{
#' sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
#' specs <- data.frame(
#'   name    = c("Lipophilicity", "FractionUnbound"),
#'   pattern = c("Aciclovir\\|Lipophilicity$",
#'               "Aciclovir\\|Fraction unbound \\(plasma\\)$"),
#'   dist    = c("normal", "truncnorm"),
#'   width   = c(0.3, 0.15),
#'   lower   = c(NA, 1e-4), upper = c(NA, 0.999)
#' )
#' params <- ospsuite_parameters(sim, specs)
#' }
#' @export
ospsuite_parameters <- function(simulation, specs, dummy = "dummy") {
  sim <- .as_simulation(simulation)
  specs <- as.data.frame(specs, stringsAsFactors = FALSE)
  required <- c("name", "pattern", "dist", "width")
  missing_cols <- setdiff(required, names(specs))
  if (length(missing_cols)) {
    .stop(c(
      "{.arg specs} is missing column{?s}: {.val {missing_cols}}.",
      "i" = "Required: {.val {required}}; optional: {.val lower}, {.val upper}."
    ))
  }
  if (!nrow(specs)) {
    .stop(c("{.arg specs} has no rows.", "i" = "Provide at least one parameter specification."))
  }
  if (!"lower" %in% names(specs)) {
    specs$lower <- NA_real_
  }
  if (!"upper" %in% names(specs)) {
    specs$upper <- NA_real_
  }

  all_paths <- ospsuite::getAllParameterPathsIn(sim)
  paths <- vapply(
    seq_len(nrow(specs)),
    function(i) {
      .resolve_unique(specs$name[i], specs$pattern[i], all_paths)
    },
    character(1)
  )
  nominal <- vapply(
    paths,
    function(p) {
      as.numeric(ospsuite::getParameter(p, sim)$value)
    },
    numeric(1)
  )

  bad <- which(specs$dist == "lognormal" & (!is.finite(nominal) | nominal <= 0))
  if (length(bad)) {
    .stop(c(
      "Lognormal parameters need a strictly positive nominal value.",
      "x" = "Non-positive / non-finite nominal for {.val {specs$name[bad]}}.",
      "i" = "Pick a different distribution or check the resolved path."
    ))
  }

  plist <- lapply(seq_len(nrow(specs)), function(i) {
    .gsa_anchor_parameter(
      specs$name[i],
      specs$dist[i],
      nominal[i],
      specs$width[i],
      specs$lower[i],
      specs$upper[i],
      path = paths[i]
    )
  })
  if (!is.null(dummy) && !is.na(dummy) && nzchar(dummy)) {
    plist <- c(plist, list(gsa_dummy(dummy)))
  }
  params <- gsa_parameters(plist)
  attr(params, "resolved") <- data.frame(
    name = specs$name,
    dist = specs$dist,
    nominal = nominal,
    path = paths,
    stringsAsFactors = FALSE
  )
  params
}

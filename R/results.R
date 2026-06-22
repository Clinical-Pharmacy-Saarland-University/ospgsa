.idx_cols <- c(
  "output",
  "time",
  "parameter",
  "method",
  "index",
  "stage",
  "estimate",
  "bias",
  "std_error",
  "conf_low",
  "conf_high",
  "conf_level",
  "significant",
  "rank"
)

.indices_dt <- function(...) {
  d <- data.table::data.table(...)
  for (cn in .idx_cols) {
    if (!cn %in% names(d)) {
      d[[cn]] <- switch(
        cn,
        output = NA_character_,
        parameter = NA_character_,
        method = NA_character_,
        index = NA_character_,
        stage = NA_character_,
        significant = NA,
        rank = NA_integer_,
        NA_real_
      )
    }
  }
  data.table::setcolorder(d, .idx_cols)
  d[]
}

.finalize_indices <- function(d, dummy = NULL) {
  if (!nrow(d)) {
    return(d)
  }
  grp <- c("output", "time", "method", "index", "stage")
  d[,
    rank := {
      o <- order(-estimate, na.last = TRUE)
      r <- integer(.N)
      r[o] <- seq_len(.N)
      r
    },
    by = grp
  ]
  ci_zero <- function(lo, hi, est) {
    s <- (is.finite(lo) & lo > 0) | (is.finite(hi) & hi < 0)
    s[is.na(lo) & is.na(hi)] <- is.finite(est[is.na(lo) & is.na(hi)]) &
      abs(est[is.na(lo) & is.na(hi)]) > 0.01
    s
  }
  if (!is.null(dummy)) {
    d[,
      significant := {
        if (any(parameter == dummy)) {
          # Significant iff this CI does not overlap the dummy's CI band (MC noise floor).
          dest <- estimate[parameter == dummy][1]
          dhi <- conf_high[parameter == dummy][1]
          dlo <- conf_low[parameter == dummy][1]
          hi <- if (is.finite(dhi)) dhi else dest
          lo <- if (is.finite(dlo)) dlo else dest
          ifelse(
            is.finite(conf_low) & is.finite(conf_high),
            conf_low > hi | conf_high < lo,
            is.finite(estimate) & abs(estimate) > abs(if (is.finite(dest)) dest else 0)
          )
        } else {
          ci_zero(conf_low, conf_high, estimate)
        }
      },
      by = grp
    ]
  } else {
    d[, significant := ci_zero(conf_low, conf_high, estimate)]
  }
  d[]
}

#' The `ospgsa_result` object
#'
#' The common return value of the GSA estimators and [gsa()]. It is a list with:
#' \describe{
#'   \item{`indices`}{a tidy `data.table`, one row per
#'     output x parameter x method x index x stage, with columns `estimate`,
#'     `bias`, `std_error`, `conf_low`, `conf_high`, `conf_level`, `significant`
#'     and `rank`.}
#'   \item{`design`}{the evaluated input matrix (when available).}
#'   \item{`Y`}{the output matrix (when available).}
#'   \item{`meta`}{run metadata (method, sample size, failures, ...).}
#' }
#' Available methods: [print()], [summary()], [as.data.frame()] and [plot()]
#' (see [plot.ospgsa_result()]). Use [gsa_table()] for a formatted view.
#'
#' @name ospgsa_result
#' @seealso [gsa()], [gsa_delta()], [gsa_sobol()], [gsa_morris()], [gsa_regression()]
NULL

new_ospgsa_result <- function(indices, design = NULL, Y = NULL, meta = list()) {
  indices <- data.table::as.data.table(indices)
  structure(list(indices = indices[], design = design, Y = Y, meta = meta), class = "ospgsa_result")
}

#' @export
print.ospgsa_result <- function(x, ...) {
  m <- x$meta
  methods <- sort(unique(x$indices$method))
  outs <- sort(unique(stats::na.omit(x$indices$output)))
  cat("<ospgsa_result>\n")
  cat(sprintf("  method(s): %s\n", paste(methods, collapse = ", ")))
  cat(sprintf("  output(s): %s\n", paste(outs, collapse = ", ")))
  if (!is.null(m$n_runs)) {
    cat(
      sprintf("  model runs: %d", m$n_runs),
      if (!is.null(m$n_failed)) sprintf(" (%d failed)", m$n_failed) else "",
      "\n",
      sep = ""
    )
  }
  np <- length(unique(stats::na.omit(x$indices$parameter)))
  cat(sprintf("  parameters: %d\n", np))
  d <- x$indices[is.na(time)]
  if (nrow(d)) {
    show <- d[order(output, method, index, rank)]
    show <- show[, utils::head(.SD, 6), by = .(output, method, index, stage)]
    keep <- c(
      "output",
      "parameter",
      "method",
      "index",
      "stage",
      "estimate",
      "conf_low",
      "conf_high",
      "significant",
      "rank"
    )
    print(show[, ..keep], nrows = 60)
  }
  invisible(x)
}

#' @export
summary.ospgsa_result <- function(object, ...) {
  d <- object$indices
  out <- d[,
    .(
      n_param = length(unique(parameter)),
      n_significant = sum(significant, na.rm = TRUE),
      max_index = max(estimate, na.rm = TRUE)
    ),
    by = .(output, method, index, stage)
  ]
  data.table::setorder(out, output, method, index, stage)
  out[]
}

#' Coerce GSA result indices to a data frame
#' @param x An `ospgsa_result`.
#' @param ... Unused.
#' @return The tidy indices `data.frame`.
#' @export
as.data.frame.ospgsa_result <- function(x, ...) as.data.frame(x$indices)

#' @export
as.data.table.ospgsa_result <- function(x, keep.rownames = FALSE, ...) {
  data.table::copy(x$indices)
}

.combine_results <- function(results) {
  results <- Filter(Negate(is.null), results)
  if (!length(results)) {
    .stop("There are no results to combine.")
  }
  idx <- data.table::rbindlist(lapply(results, function(r) r$indices), fill = TRUE)
  meta <- list(
    combined = TRUE,
    methods = unique(unlist(lapply(results, function(r) r$meta$method)))
  )
  new_ospgsa_result(idx, meta = meta)
}

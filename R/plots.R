.theme_gsa <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92")
    )
}

.filter_scalar <- function(d, output = NULL, method = NULL, index = NULL) {
  d <- d[is.na(d$time), ]
  if (!is.null(output)) {
    d <- d[d$output %in% output, ]
  }
  if (!is.null(method)) {
    d <- d[d$method %in% method, ]
  }
  if (!is.null(index)) {
    d <- d[d$index %in% index, ]
  }
  d[d$parameter != "__model__", ]
}

.order_param <- function(d) {
  ord <- stats::aggregate(estimate ~ parameter, data = d, FUN = function(x) {
    max(x, na.rm = TRUE)
  })
  ord <- ord[order(ord$estimate), ]
  d$parameter <- factor(d$parameter, levels = ord$parameter)
  d
}

#' Bar plot of sensitivity indices with confidence intervals
#'
#' @param result An [ospgsa_result].
#' @param index Which index to show (e.g. `"delta"`, `"ST"`, `"mu_star"`,
#'   `"PRCC"`). Defaults to the first available.
#' @param output Optional output(s) to include (faceted if several).
#' @param method Optional method filter.
#' @param threshold Optional dashed reference line (e.g. an importance cutoff).
#' @param labels Optional named character vector mapping parameter names to
#'   display labels for the axis (parameters absent from the map keep their name).
#' @param show_values If `TRUE`, print each bar's estimate at the bar end
#'   (sign-aware placement; matching a local-SA style figure).
#' @param digits Decimal places for the printed values.
#' @return A `ggplot` object.
#' @export
plot_indices <- function(
  result,
  index = NULL,
  output = NULL,
  method = NULL,
  threshold = NULL,
  labels = NULL,
  show_values = FALSE,
  digits = 3L
) {
  d <- as.data.frame(result$indices)
  if (is.null(index)) {
    index <- d$index[!d$index %in% c("R2_SRC", "R2_SRRC", "S1")][1L]
  }
  d <- .filter_scalar(d, output, method, index)
  if (!nrow(d)) {
    avail <- setdiff(unique(as.data.frame(result$indices)$index), c("R2_SRC", "R2_SRRC"))
    .stop(c(
      "Nothing to plot for index {.val {index}}.",
      "i" = "Available index value{?s}: {.val {avail}}."
    ))
  }
  if (!is.null(labels)) {
    hit <- d$parameter %in% names(labels)
    d$parameter[hit] <- unname(labels[d$parameter[hit]])
  }
  d <- .order_param(d)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$parameter, y = .data$estimate)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$significant), width = 0.7) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#2c7fb8", `FALSE` = "grey70"),
      na.value = "grey70",
      name = "significant"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = index, title = sprintf("%s sensitivity index", index)) +
    .theme_gsa()
  if (all(c("conf_low", "conf_high") %in% names(d)) && any(is.finite(d$conf_low))) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
        width = 0.25
      )
  }
  if (isTRUE(show_values)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = formatC(.data$estimate, format = "f", digits = digits)),
        hjust = ifelse(d$estimate >= 0, -0.15, 1.15),
        size = 2.5
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.18)))
  }
  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = threshold, linetype = 2, colour = "red")
  }
  if (length(unique(d$output)) > 1L) {
    p <- p + ggplot2::facet_wrap(~output, scales = "free_x")
  }
  p
}

#' Grouped Sobol first- vs total-order bar plot
#'
#' @param result An [ospgsa_result] containing Sobol indices.
#' @param output Optional output filter (faceted if several).
#' @return A `ggplot` object (the S1/ST gap indicates interactions).
#' @export
plot_sobol <- function(result, output = NULL) {
  d <- .filter_scalar(
    as.data.frame(result$indices),
    output,
    method = "sobol",
    index = c("S1", "ST")
  )
  if (!nrow(d)) {
    .stop(c(
      "This result contains no Sobol indices.",
      "i" = "Compute them with {.fn gsa_sobol} or {.code gsa(..., method = \"sobol\")}."
    ))
  }
  d <- .order_param(d)
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$parameter, y = .data$estimate, fill = .data$index)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.6) +
    ggplot2::scale_fill_manual(values = c(S1 = "#a6bddb", ST = "#2c7fb8"), name = NULL) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Sobol index", title = "Sobol first- vs total-order") +
    .theme_gsa()
  if (any(is.finite(d$conf_low))) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
        position = ggplot2::position_dodge(0.7),
        width = 0.25
      )
  }
  if (length(unique(d$output)) > 1L) {
    p <- p + ggplot2::facet_wrap(~output, scales = "free_x")
  }
  p
}

#' Two-stage delta plot (independent vs full)
#'
#' @param result An [ospgsa_result] from [gsa_delta_two_stage()] or [gsa()].
#' @param output Optional output filter.
#' @param labels Optional named character vector mapping parameter names (as in
#'   `result$indices$parameter`) to display labels for the axis (e.g. publication
#'   names). Parameters absent from the map keep their original name.
#' @param show_values If `TRUE`, print each bar's delta estimate at the bar end
#'   (matching a local-SA style figure).
#' @param digits Decimal places for the printed values.
#' @return A `ggplot` object comparing delta_1 (structural) and delta_2 (full).
#' @export
plot_delta_two_stage <- function(
  result,
  output = NULL,
  labels = NULL,
  show_values = FALSE,
  digits = 3L
) {
  d <- as.data.frame(result$indices)
  d <- d[is.na(d$time) & d$index == "delta" & d$stage %in% c("independent", "full"), ]
  if (!is.null(output)) {
    d <- d[d$output %in% output, ]
  }
  if (!nrow(d)) {
    .stop(c(
      "This result contains no two-stage delta indices.",
      "i" = "Compute them with {.fn gsa_delta_two_stage} or {.code gsa(..., method = \"delta\", correlation = R)}."
    ))
  }
  if (!is.null(labels)) {
    hit <- d$parameter %in% names(labels)
    d$parameter[hit] <- unname(labels[d$parameter[hit]])
  }
  d$stage <- factor(
    d$stage,
    levels = c("independent", "full"),
    labels = c("delta[1] (structural)", "delta[2] (full)")
  )
  d <- .order_param(d)
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$parameter, y = .data$estimate, fill = .data$stage)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.6) +
    ggplot2::scale_fill_manual(
      values = c("grey70", "#2c7fb8"),
      name = NULL,
      labels = c(expression(delta[1] ~ "(structural)"), expression(delta[2] ~ "(full)"))
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = expression(delta), title = "Two-stage delta analysis") +
    .theme_gsa()
  if (any(is.finite(d$conf_low))) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
        position = ggplot2::position_dodge(0.7),
        width = 0.25
      )
  }
  if (isTRUE(show_values)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          y = ifelse(is.finite(.data$conf_high), .data$conf_high, .data$estimate),
          label = formatC(.data$estimate, format = "f", digits = digits)
        ),
        position = ggplot2::position_dodge(0.7),
        hjust = -0.15,
        size = 2.5
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.20)))
  }
  if (length(unique(d$output)) > 1L) {
    p <- p + ggplot2::facet_wrap(~output, scales = "free_x")
  }
  p
}

#' Morris (mu*, sigma) screening plot
#'
#' @param result An [ospgsa_result] with Morris indices.
#' @param output Optional output filter.
#' @return A `ggplot` object; high mu* = important, high sigma = nonlinear/interacting.
#' @export
plot_morris <- function(result, output = NULL) {
  d <- as.data.frame(result$indices)
  d <- d[
    is.na(d$time) &
      d$method == "morris" &
      d$index %in% c("mu_star", "sigma") &
      d$parameter != "__model__",
  ]
  if (!is.null(output)) {
    d <- d[d$output %in% output, ]
  }
  if (!nrow(d)) {
    .stop(c(
      "This result contains no Morris indices.",
      "i" = "Compute them with {.fn gsa_morris} or {.code gsa(..., method = \"morris\")}."
    ))
  }
  w <- stats::reshape(
    d[, c("output", "parameter", "index", "estimate")],
    idvar = c("output", "parameter"),
    timevar = "index",
    direction = "wide"
  )
  names(w) <- sub("estimate\\.", "", names(w))
  p <- ggplot2::ggplot(
    w,
    ggplot2::aes(x = .data$mu_star, y = .data$sigma, label = .data$parameter)
  ) +
    ggplot2::geom_abline(slope = 1, linetype = 3, colour = "grey50") +
    ggplot2::geom_point(colour = "#2c7fb8", size = 2) +
    ggplot2::geom_text(vjust = -0.6, size = 3) +
    ggplot2::labs(
      x = expression(mu * "*"),
      y = expression(sigma),
      title = "Morris elementary effects"
    ) +
    .theme_gsa()
  if (length(unique(w$output)) > 1L) {
    p <- p + ggplot2::facet_wrap(~output, scales = "free")
  }
  p
}

#' Time-resolved index heatmap
#'
#' @param result An [ospgsa_result] with time-resolved indices (non-`NA` `time`).
#' @param index Index to display.
#' @param output Optional output filter.
#' @return A `ggplot` tile heatmap (parameters x time).
#' @export
plot_time_heatmap <- function(result, index = "delta", output = NULL) {
  d <- as.data.frame(result$indices)
  d <- d[!is.na(d$time) & d$index == index & d$parameter != "__model__", ]
  if (!is.null(output)) {
    d <- d[d$output %in% output, ]
  }
  if (!nrow(d)) {
    .stop(c(
      "No time-resolved {.val {index}} indices in this result.",
      "i" = "Build a time-resolved evaluator with {.fn ospsuite_evaluator} using {.arg time_points}."
    ))
  }
  ggplot2::ggplot(d, ggplot2::aes(x = .data$time, y = .data$parameter, fill = .data$estimate)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = index) +
    ggplot2::labs(x = "time", y = NULL, title = sprintf("Time-resolved %s", index)) +
    .theme_gsa() +
    (if (length(unique(d$output)) > 1L) ggplot2::facet_wrap(~output) else NULL)
}

#' Convergence plot
#'
#' @param conv An `ospgsa_convergence` object from [gsa_convergence()].
#' @param output Optional output filter.
#' @return A `ggplot` of the tracked index vs sample size with CI ribbons.
#' @export
plot_convergence <- function(conv, output = NULL) {
  d <- as.data.frame(conv)
  if (!is.null(output)) {
    d <- d[d$output %in% output, ]
  }
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$n, y = .data$estimate, colour = .data$parameter, fill = .data$parameter)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
      alpha = 0.15,
      colour = NA
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "sample size (N)", y = unique(d$index)[1L], title = "GSA convergence") +
    .theme_gsa() +
    (if (length(unique(d$output)) > 1L) {
      ggplot2::facet_wrap(~output, scales = "free_y")
    } else {
      NULL
    })
}

#' Scatter plot of output vs each input
#'
#' @param X Input design ([gsa_sample()] or matrix).
#' @param y Output vector.
#' @param max_points Subsample for speed.
#' @return A faceted `ggplot` scatter (a raw sanity check before indices).
#' @export
plot_scatter <- function(X, y, max_points = 2000L) {
  des <- .coerce_design(X)
  Xm <- des$X
  if (nrow(Xm) > max_points) {
    s <- sample.int(nrow(Xm), max_points)
    Xm <- Xm[s, , drop = FALSE]
    y <- y[s]
  }
  long <- data.table::rbindlist(lapply(colnames(Xm), function(cn) {
    data.table::data.table(parameter = cn, value = Xm[, cn], y = y)
  }))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$y)) +
    ggplot2::geom_point(alpha = 0.2, size = 0.6, colour = "#2c7fb8") +
    ggplot2::facet_wrap(~parameter, scales = "free_x") +
    ggplot2::labs(x = NULL, y = "output") +
    .theme_gsa()
}

#' Plot method for ospgsa results
#'
#' @param x An [ospgsa_result].
#' @param type One of `"indices"`, `"sobol"`, `"delta_two_stage"`, `"morris"`,
#'   `"time"`.
#' @param ... Passed to the specific plotting function.
#' @return A `ggplot` object.
#' @export
plot.ospgsa_result <- function(
  x,
  type = c("indices", "sobol", "delta_two_stage", "morris", "time"),
  ...
) {
  type <- match.arg(type)
  switch(
    type,
    indices = plot_indices(x, ...),
    sobol = plot_sobol(x, ...),
    delta_two_stage = plot_delta_two_stage(x, ...),
    morris = plot_morris(x, ...),
    time = plot_time_heatmap(x, ...)
  )
}

#' @rdname plot.ospgsa_result
#' @export
gsa_plot <- function(x, type = "indices", ...) {
  plot.ospgsa_result(x, type = type, ...)
}

#' Save the standard GSA plots as PNG files
#'
#' Writes one PNG per output for each index family present in `result`: the
#' two-stage delta split (when a `"full"` stage is present), `PRCC`/`SRC`
#' regression bars, Sobol `S1`/`ST`, and the Morris `mu*` plane. Plots that
#' cannot be produced are skipped (not fatal). This packages the plot-writing
#' loop that GSA driver scripts otherwise repeat by hand.
#'
#' @param result An [ospgsa_result].
#' @param dir Output directory (created if needed).
#' @param prefix File-name prefix.
#' @param width,height,dpi Passed to [ggplot2::ggsave()].
#' @param quiet Suppress the per-file progress messages.
#' @return (Invisibly) the character vector of files written.
#' @seealso [plot_delta_two_stage()], [plot_indices()], [plot_sobol()]
#' @examples
#' p <- gsa_parameters(
#'   gsa_parameter("a", "uniform", min = 0, max = 1),
#'   gsa_parameter("b", "uniform", min = 0, max = 1)
#' )
#' f <- function_evaluator(function(M) 2 * M[, "a"] + M[, "b"], "y")
#' res <- gsa(p, f, method = "sobol", n = 256, boot = 0)
#' gsa_save_plots(res, dir = tempfile("plots"))
#' @export
gsa_save_plots <- function(
  result,
  dir = ".",
  prefix = "plot_",
  width = 8,
  height = 5,
  dpi = 150,
  quiet = FALSE
) {
  if (!inherits(result, "ospgsa_result")) {
    .stop(c(
      "{.arg result} must be an {.cls ospgsa_result}.",
      "x" = "You supplied {.obj_desc {result}}."
    ))
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  idx <- unique(result$indices$index)
  outputs <- sort(unique(stats::na.omit(result$indices$output)))
  if (!length(outputs)) {
    outputs <- NA_character_
  }
  has_two_stage <- "delta" %in% idx && any(result$indices$stage %in% "full", na.rm = TRUE)
  written <- character(0)
  save1 <- function(p, file) {
    if (is.null(p)) {
      return(invisible())
    }
    ok <- tryCatch(
      {
        ggplot2::ggsave(file.path(dir, file), plot = p, width = width, height = height, dpi = dpi)
        TRUE
      },
      error = function(e) FALSE
    )
    if (ok) {
      written <<- c(written, file)
      if (!quiet) cli::cli_inform(c("v" = "wrote {.path {file}}"))
    }
  }
  for (oo in outputs) {
    safe <- if (is.na(oo)) "all" else gsub("[^A-Za-z0-9]+", "_", oo)
    arg_out <- if (is.na(oo)) NULL else oo
    if (has_two_stage) {
      save1(
        tryCatch(plot_delta_two_stage(result, output = arg_out), error = function(e) NULL),
        sprintf("%sdelta_two_stage_%s.png", prefix, safe)
      )
    }
    for (ix in intersect(c("PRCC", "SRC"), idx)) {
      save1(
        tryCatch(plot_indices(result, index = ix, output = arg_out), error = function(e) NULL),
        sprintf("%s%s_%s.png", prefix, tolower(ix), safe)
      )
    }
    if (all(c("S1", "ST") %in% idx)) {
      save1(
        tryCatch(plot_sobol(result, output = arg_out), error = function(e) NULL),
        sprintf("%ssobol_%s.png", prefix, safe)
      )
    }
    if ("mu_star" %in% idx) {
      save1(
        tryCatch(plot_morris(result, output = arg_out), error = function(e) NULL),
        sprintf("%smorris_%s.png", prefix, safe)
      )
    }
  }
  invisible(written)
}

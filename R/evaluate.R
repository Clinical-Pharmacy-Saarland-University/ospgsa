#' Wrap an R function as a GSA model evaluator
#'
#' @param fun A function taking a numeric matrix `X` (n x k, columns named by
#'   parameter) and returning either a numeric vector of length `n` (single
#'   output) or a numeric matrix `n x m` (multiple outputs).
#' @param output_names Optional names for the outputs.
#' @param vectorized If `TRUE` (default) `fun` is called once on the whole
#'   matrix; if `FALSE` it is applied row by row.
#' @return A function with the evaluator contract, usable by [gsa_sobol()],
#'   [gsa_morris()] and [gsa()].
#' @examples
#' f <- function(M) sin(M[, 1]) + M[, 2]^2
#' ev <- function_evaluator(f, output_names = "y")
#' ev(matrix(runif(20), ncol = 2, dimnames = list(NULL, c("x1", "x2"))))
#' @export
function_evaluator <- function(fun, output_names = NULL, vectorized = TRUE) {
  if (!is.function(fun)) {
    desc <- .obj_desc(fun)
    .stop(c(
      "{.arg fun} must be a function mapping an input matrix to outputs.",
      "x" = "You supplied {desc}."
    ))
  }
  function(X) {
    X <- as.matrix(X)
    Y <- if (vectorized) {
      fun(X)
    } else {
      tmp <- lapply(seq_len(nrow(X)), function(i) fun(X[i, , drop = FALSE]))
      do.call(rbind, lapply(tmp, function(v) matrix(v, nrow = 1L)))
    }
    if (is.null(dim(Y))) {
      Y <- matrix(Y, ncol = 1L)
    }
    if (!is.null(output_names)) {
      colnames(Y) <- output_names
    }
    if (is.null(colnames(Y))) {
      colnames(Y) <- if (ncol(Y) == 1L) "y" else paste0("y", seq_len(ncol(Y)))
    }
    Y
  }
}

#' Evaluate a design with an evaluator (optionally checkpointed and resumable)
#'
#' Runs a design matrix through a model evaluator (see [function_evaluator()],
#' [ospsuite_evaluator()]) and collects the outputs.
#'
#' For long or fragile runs -- thousands of PK-Sim solves that can take hours --
#' pass `checkpoint_dir` to evaluate the design in **resumable blocks**. Each
#' block of `block_size` rows is written to disk as soon as it finishes, keyed by
#' a content hash of the design; a re-run loads the completed blocks instead of
#' recomputing them. A crash, an interruption, or a deliberate stop therefore
#' costs at most one block, and re-running the same call picks up where it left
#' off. This is the engine the bundled `gsa.R` driver scripts use.
#'
#' @param X A [gsa_sample()] or a numeric input matrix (natural units).
#' @param evaluator A model evaluator (see [function_evaluator()],
#'   [ospsuite_evaluator()]).
#' @param checkpoint_dir Optional path. When `NULL` (default) the whole design is
#'   evaluated in a single call (no checkpointing). When a directory is given the
#'   evaluation is blocked and resumable; the directory is created if needed.
#' @param block_size Rows evaluated and checkpointed per block (only used when
#'   `checkpoint_dir` is set). Smaller blocks checkpoint more often -- less work
#'   lost on a crash -- at negligible overhead.
#' @param tag Short prefix for this design's checkpoint files. Use distinct tags
#'   to keep several designs' checkpoints (e.g. an independent and a correlated
#'   design) side by side in one directory.
#' @param crash_skip If `TRUE`, guard each block with a marker file so that a
#'   block which *hard-crashes* the R process (e.g. a solver segfault on a stiff
#'   parameter draw, which cannot be caught with [tryCatch()]) is filled with
#'   `NA` rows on the next resume instead of crashing again -- the estimators
#'   drop `NA` rows. Off by default because, while on, interrupting a block
#'   mid-evaluation also turns that block into `NA` on resume. Turn it on for
#'   unattended runs of models that crash on some draws, and re-run until
#'   completion (e.g. in a loop / watchdog). Requires `checkpoint_dir`.
#' @param quiet Suppress per-block progress messages.
#' @return A list with `X` (the input matrix), `Y` (the output matrix) and
#'   `n_failed` (rows with any non-finite output). When `checkpoint_dir` is set
#'   it additionally carries `n_blocks`, `checkpoint_dir`, `tag` and `key` (the
#'   design content hash used in the checkpoint file names).
#' @seealso [ospsuite_evaluator()], [gsa()]
#' @examples
#' f <- function_evaluator(function(M) M[, 1] + M[, 2]^2, "y")
#' X <- matrix(runif(20), ncol = 2, dimnames = list(NULL, c("a", "b")))
#' gsa_evaluate(X, f)$Y
#'
#' # Resumable, crash-safe evaluation of a long run. A second call with the same
#' # design and directory loads the finished blocks instead of recomputing them.
#' dir <- tempfile("ckpt")
#' r1 <- gsa_evaluate(X, f, checkpoint_dir = dir, block_size = 5, tag = "demo")
#' r2 <- gsa_evaluate(X, f, checkpoint_dir = dir, block_size = 5, tag = "demo")
#' identical(r1$Y, r2$Y)
#' @export
gsa_evaluate <- function(
  X,
  evaluator,
  checkpoint_dir = NULL,
  block_size = 1000L,
  tag = "block",
  crash_skip = FALSE,
  quiet = FALSE
) {
  des <- .coerce_design(X)
  Xm <- des$X
  .check_evaluator(evaluator)
  if (is.null(checkpoint_dir)) {
    Y <- .coerce_outputs(evaluator(Xm))
    if (nrow(Y) != nrow(Xm)) {
      .abort_eval_rows(nrow(Y), nrow(Xm))
    }
    return(list(X = Xm, Y = Y, n_failed = sum(!stats::complete.cases(Y))))
  }
  .assert_string(checkpoint_dir)
  block_size <- .assert_count(block_size, min = 1L)
  .assert_string(tag)
  .assert_flag(crash_skip)
  .assert_flag(quiet)
  .gsa_evaluate_checkpointed(Xm, evaluator, checkpoint_dir, block_size, tag, crash_skip, quiet)
}

#' Set non-positive outputs to NA for log-scale analysis
#'
#' AUC / Cmax and similar PK metrics are analysed on the log scale
#' (`log_output = TRUE`), which is undefined for non-positive values. A stiff or
#' failed model run can occasionally return a zero or negative output; this
#' helper sets those to `NA` (the estimators drop `NA` rows) so a few bad values
#' do not contaminate the log transform. Accepts either an output matrix or an
#' evaluation result (the list returned by [gsa_evaluate()]).
#'
#' @param Y A numeric output matrix, or an evaluation result list with a `Y`.
#' @param quiet Suppress the count message.
#' @return The same shape as `Y` with non-positive entries set to `NA` (and, for
#'   an evaluation result, `n_failed` refreshed).
#' @seealso [gsa_evaluate()]
#' @examples
#' Y <- matrix(c(1, -2, 3, 0), ncol = 1)
#' gsa_sanitize_positive(Y, quiet = TRUE)
#' @export
gsa_sanitize_positive <- function(Y, quiet = FALSE) {
  is_result <- is.list(Y) && !is.null(Y$Y)
  M <- if (is_result) Y$Y else Y
  M <- as.matrix(M)
  bad <- is.finite(M) & M <= 0
  nbad <- sum(bad)
  if (nbad && !quiet) {
    cli::cli_inform(c(
      "!" = "Set {nbad} non-positive output value{?s} to {.code NA} before log analysis."
    ))
  }
  M[bad] <- NA_real_
  if (is_result) {
    Y$Y <- M
    Y$n_failed <- sum(!stats::complete.cases(M))
    return(Y)
  }
  M
}

.gsa_design_key <- function(X, call = parent.frame()) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    .stop(
      c(
        "Package {.pkg digest} is required for checkpointed evaluation.",
        "i" = "Install it with {.code install.packages(\"digest\")}."
      ),
      call = call
    )
  }
  h <- digest::digest(list(dim(X), colnames(X), unname(X)), algo = "xxhash64")
  sprintf("%s_%d", h, nrow(X))
}

.gsa_block_cols <- function(known_cols, dir, call = parent.frame()) {
  if (!is.null(known_cols)) {
    return(known_cols)
  }
  for (f in list.files(dir, pattern = "\\.rds$", full.names = TRUE)) {
    Yb <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.matrix(Yb) && ncol(Yb) > 0L && !is.null(colnames(Yb))) {
      return(colnames(Yb))
    }
  }
  .stop(
    c(
      "Cannot fill a crashed block with {.code NA}: the output column layout is unknown.",
      "i" = "No completed checkpoint block exists yet to infer the columns from.",
      "i" = "Evaluate at least one block successfully first, or use {.code crash_skip = FALSE}."
    ),
    call = call
  )
}

.gsa_evaluate_checkpointed <- function(X, evaluator, dir, block_size, tag, crash_skip, quiet) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  n <- nrow(X)
  key <- .gsa_design_key(X)
  nblk <- as.integer(ceiling(n / block_size))
  say <- function(msg) if (!quiet) cli::cli_inform(msg, .envir = parent.frame())
  say(c(
    "i" = paste0(
      "Checkpointed evaluation of {n} row{?s} in {nblk} ",
      "block{?s} of <= {block_size} (tag {.val {tag}})."
    )
  ))

  blocks <- vector("list", nblk)
  known_cols <- NULL
  block_file <- function(b) {
    # block_size in the name: changing it starts a fresh set, not misaligned blocks.
    file.path(dir, sprintf("%s_%s_bs%d_block%04d.rds", tag, key, block_size, b))
  }

  for (b in seq_len(nblk)) {
    i0 <- (b - 1L) * block_size + 1L
    i1 <- min(b * block_size, n)
    nb <- i1 - i0 + 1L
    fp <- block_file(b)
    mk <- paste0(fp, ".crashed")

    if (file.exists(fp)) {
      Yb <- tryCatch(readRDS(fp), error = function(e) NULL)
      ok <- is.matrix(Yb) &&
        nrow(Yb) == nb &&
        (is.null(known_cols) || identical(colnames(Yb), known_cols))
      if (ok) {
        if (is.null(known_cols)) {
          known_cols <- colnames(Yb)
        }
        if (crash_skip && file.exists(mk)) {
          file.remove(mk)
        }
        blocks[[b]] <- Yb
        say(c(">" = "block {b}/{nblk}: loaded {nb} row{?s} from checkpoint."))
        next
      }
      say(c("!" = "block {b}/{nblk}: stale checkpoint -> recomputing."))
    }

    # A surviving crash marker means the previous attempt died mid-block.
    if (crash_skip && file.exists(mk)) {
      cols <- .gsa_block_cols(known_cols, dir)
      Yb <- matrix(NA_real_, nb, length(cols), dimnames = list(NULL, cols))
      saveRDS(Yb, fp)
      file.remove(mk)
      known_cols <- cols
      blocks[[b]] <- Yb
      say(c("x" = "block {b}/{nblk}: previous attempt crashed -> {nb} NA row{?s} skipped."))
      next
    }

    if (crash_skip) {
      file.create(mk)
    }
    Yb <- .coerce_outputs(evaluator(X[i0:i1, , drop = FALSE]))
    if (nrow(Yb) != nb) {
      .abort_eval_rows(nrow(Yb), nb)
    }
    saveRDS(Yb, fp)
    if (crash_skip && file.exists(mk)) {
      file.remove(mk)
    }
    if (is.null(known_cols)) {
      known_cols <- colnames(Yb)
    }
    blocks[[b]] <- Yb
    nf <- sum(!stats::complete.cases(Yb))
    say(c("v" = "block {b}/{nblk}: computed {nb} row{?s} ({nf} failed) -> saved."))
    invisible(gc(FALSE))
  }

  Y <- do.call(rbind, blocks)
  if (nrow(Y) != n) {
    .stop(c(
      "Checkpoint reassembly produced the wrong number of rows.",
      "x" = "Expected {n} but got {nrow(Y)}.",
      "i" = "Delete the checkpoint directory {.path {dir}} to recompute from scratch."
    ))
  }
  list(
    X = X,
    Y = Y,
    n_failed = sum(!stats::complete.cases(Y)),
    n_blocks = nblk,
    checkpoint_dir = dir,
    tag = tag,
    key = key
  )
}

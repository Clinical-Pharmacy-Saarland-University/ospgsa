# Evaluator contract, output hygiene, tables, plots.

test_that("function_evaluator vectorized and non-vectorized agree, and name outputs", {
  f <- function(M) cbind(M[, 1] + M[, 2], M[, 1] * M[, 2])
  X <- matrix(stats::runif(20), ncol = 2, dimnames = list(NULL, c("a", "b")))
  ev_v <- function_evaluator(f, output_names = c("s", "p"), vectorized = TRUE)
  ev_r <- function_evaluator(f, output_names = c("s", "p"), vectorized = FALSE)
  expect_equal(ev_v(X), ev_r(X))
  expect_equal(colnames(ev_v(X)), c("s", "p"))

  e1 <- function_evaluator(function(M) M[, 1])
  e2 <- function_evaluator(function(M) cbind(M[, 1], M[, 2]))
  X2 <- matrix(stats::runif(6), ncol = 2)
  expect_equal(colnames(e1(X2)), "y")
  expect_equal(colnames(e2(X2)), c("y1", "y2"))
  expect_error(function_evaluator(42), "must be a function")
})

test_that("gsa_evaluate errors when the evaluator returns the wrong row count", {
  X <- matrix(stats::runif(10), ncol = 2, dimnames = list(NULL, c("a", "b")))
  bad <- function(M) matrix(0, nrow(M) - 1L, 1L, dimnames = list(NULL, "y"))
  expect_error(gsa_evaluate(X, bad), "row")
})

test_that("gsa_sanitize_positive reports a count and no-ops on clean input", {
  expect_message(gsa_sanitize_positive(matrix(c(1, -1), ncol = 1)), "non-positive")
  M <- matrix(c(1, 2, 3), ncol = 1)
  expect_identical(gsa_sanitize_positive(M, quiet = TRUE), M)
})

test_that("gsa_table filters, rounds and orders", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  ev <- function_evaluator(function(M) 3 * M[, "a"] + M[, "b"], "y")
  res <- gsa(p, ev, method = "delta", n = 800, boot = 0, seed = 1)
  tb <- gsa_table(res, index = "delta", digits = 2)
  expect_true(all(tb$index == "delta"))
  expect_true(all(is.na(tb$time)))
  expect_equal(tb$estimate, round(tb$estimate, 2))
})

test_that("plot_indices auto-picks a non-R2 index and errors on a missing one", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  ev <- function_evaluator(function(M) 2 * M[, "a"], "y")
  res <- gsa(p, ev, method = "delta", n = 600, boot = 0, seed = 1)
  expect_s3_class(plot_indices(res), "ggplot")
  expect_error(plot_indices(res, index = "ST"), "Nothing to plot")
})

test_that("plot_time_heatmap and plot_convergence build", {
  d <- .indices_dt(
    output = "y",
    method = "delta",
    index = "delta",
    parameter = rep(c("a", "b"), each = 3),
    time = rep(c(1, 2, 3), 2),
    estimate = stats::runif(6)
  )
  expect_s3_class(plot_time_heatmap(new_ospgsa_result(d), index = "delta"), "ggplot")

  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = -pi, max = pi),
    gsa_parameter("x2", "uniform", min = -pi, max = pi)
  )
  ev <- function_evaluator(function(M) sin(M[, "x1"]) + 7 * sin(M[, "x2"])^2, "y")
  conv <- gsa_convergence(p, ev, method = "sobol", n_seq = c(256, 512), boot = 0)
  expect_s3_class(plot_convergence(conv), "ggplot")
})

test_that("crash_skip infers columns from an on-disk block when known_cols is NULL", {
  X <- matrix(stats::runif(40), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
  ev <- function_evaluator(function(M) cbind(a = M[, 1], b = M[, 2]))
  dir <- tempfile("ckpt")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  gsa_evaluate(
    X,
    ev,
    checkpoint_dir = dir,
    block_size = 5,
    tag = "t",
    crash_skip = TRUE,
    quiet = TRUE
  )
  files <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))
  file.remove(files[1:2]) # drop blocks 1-2, keep a later block on disk
  file.create(paste0(files[2L], ".crashed")) # block 2 "crashed" last run
  res <- gsa_evaluate(
    X,
    ev,
    checkpoint_dir = dir,
    block_size = 5,
    tag = "t",
    crash_skip = TRUE,
    quiet = TRUE
  )
  expect_equal(colnames(res$Y), c("a", "b")) # 2 cols inferred from disk, not 1
  expect_true(all(is.na(res$Y[6:10, ])))
})

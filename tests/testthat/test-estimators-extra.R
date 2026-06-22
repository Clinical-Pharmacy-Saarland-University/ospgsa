# Estimator behaviors beyond the analytic-benchmark tests: seeds, log guard,
# correlation warning, the regression index families, bootstrap CI sanity.

test_that("gsa_delta is reproducible under a seed", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1)
  )
  s <- gsa_sample(p, 2000, seed = 3)
  y <- 3 * s$X[, "x1"] + s$X[, "x2"]
  r1 <- gsa_delta(s$X, y, boot = 25, seed = 1)
  r2 <- gsa_delta(s$X, y, boot = 25, seed = 1)
  expect_equal(r1$indices$estimate, r2$indices$estimate)
})

test_that("log_output rejects non-positive outputs", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1)
  )
  s <- gsa_sample(p, 500, seed = 3)
  y <- s$X[, "x1"] - 0.5 # spans negative values
  expect_error(gsa_delta(s$X, y, boot = 0, log_output = TRUE), "positive")
})

test_that("gsa_sobol ignores a supplied correlation with a warning", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  ev <- function_evaluator(function(M) M[, "a"] + M[, "b"], "y")
  R <- gsa_correlation(p, c("a", "b", 0.5))
  expect_warning(gsa_sobol(p, ev, n = 256, boot = 0, correlation = R), "ignores")
})

test_that("gsa_regression returns SRC, SRRC, PCC and PRCC with correct signs", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  s <- gsa_sample(p, 2000, seed = 4)
  y <- 2 * s$X[, "a"] - s$X[, "b"] + stats::rnorm(2000, 0, 0.02)
  r <- gsa_regression(s$X, y, methods = c("SRC", "SRRC", "PCC", "PRCC"), parameters = p, boot = 0)
  d <- as.data.frame(r$indices)
  for (ix in c("SRC", "SRRC", "PCC", "PRCC")) {
    v <- setNames(d$estimate[d$index == ix], d$parameter[d$index == ix])
    expect_true(all(c("a", "b") %in% names(v)), info = ix)
    expect_gt(v[["a"]], 0)
    expect_lt(v[["b"]], 0)
  }
})

test_that("bootstrap delta CIs are finite and ordered", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1)
  )
  s <- gsa_sample(p, 1500, seed = 6)
  y <- 3 * s$X[, "x1"] + s$X[, "x2"]
  r <- gsa_delta(s$X, y, boot = 40, seed = 1)
  d <- r$indices[r$indices$index == "delta", ]
  expect_true(all(is.finite(d$conf_low)))
  expect_true(all(is.finite(d$conf_high)))
  expect_true(all(d$conf_low <= d$conf_high))
  expect_true(all(is.finite(d$std_error)))
})

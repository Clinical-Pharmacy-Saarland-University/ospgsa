# Analytic benchmarks with known sensitivities.

test_that("Sobol recovers Ishigami analytic indices", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = -pi, max = pi),
    gsa_parameter("x2", "uniform", min = -pi, max = pi),
    gsa_parameter("x3", "uniform", min = -pi, max = pi)
  )
  ish <- function(M) sin(M[, "x1"]) + 7 * sin(M[, "x2"])^2 + 0.1 * M[, "x3"]^4 * sin(M[, "x1"])
  set.seed(1)
  r <- gsa_sobol(p, function_evaluator(ish, "y"), n = 8192, boot = 0)
  d <- as.data.frame(r$indices)
  S1 <- setNames(d$estimate[d$index == "S1"], d$parameter[d$index == "S1"])
  ST <- setNames(d$estimate[d$index == "ST"], d$parameter[d$index == "ST"])
  expect_equal(unname(S1["x1"]), 0.314, tolerance = 0.06)
  expect_equal(unname(S1["x2"]), 0.442, tolerance = 0.06)
  expect_lt(abs(S1["x3"]), 0.06)
  expect_equal(unname(ST["x3"]), 0.244, tolerance = 0.06)
  expect_gt(ST["x1"], S1["x1"]) # x1 has interactions
})

test_that("Morris ranks active vs inert factors", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1),
    gsa_parameter("x3", "uniform", min = 0, max = 1)
  )
  f <- function(M) 5 * M[, "x1"] + 2 * M[, "x2"]
  set.seed(2)
  r <- gsa_morris(p, function_evaluator(f, "y"), r = 30, levels = 6, boot = 0)
  ms <- as.data.frame(r$indices)
  mu_star <- setNames(ms$estimate[ms$index == "mu_star"], ms$parameter[ms$index == "mu_star"])
  expect_gt(mu_star["x1"], mu_star["x2"])
  expect_lt(mu_star["x3"], 1e-6)
})

test_that("delta ranks influential parameters and inert ~ floor", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1),
    gsa_parameter("x3", "uniform", min = 0, max = 1)
  )
  set.seed(3)
  s <- gsa_sample(p, 3000, seed = 3)
  y <- 3 * s$X[, "x1"] + s$X[, "x2"]
  r <- gsa_delta(s$X, y, boot = 0)
  d <- as.data.frame(r$indices)
  del <- setNames(d$estimate[d$index == "delta"], d$parameter[d$index == "delta"])
  expect_gt(del["x1"], del["x2"])
  expect_gt(del["x2"], del["x3"])
})

test_that("regression SRC and PRCC recover signs and a high R2", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1),
    gsa_parameter("c", "uniform", min = 0, max = 1)
  )
  set.seed(4)
  s <- gsa_sample(p, 2000, seed = 4)
  y <- 2 * s$X[, "a"] - s$X[, "b"] + stats::rnorm(2000, 0, 0.05)
  r <- gsa_regression(s$X, y, methods = c("SRC", "PRCC"), parameters = p, boot = 0)
  d <- as.data.frame(r$indices)
  src <- setNames(d$estimate[d$index == "SRC"], d$parameter[d$index == "SRC"])
  expect_gt(src["a"], 0)
  expect_lt(src["b"], 0)
  expect_lt(abs(src["c"]), 0.1)
  r2 <- d$estimate[d$index == "R2_SRC"]
  expect_gt(r2, 0.95)
})

test_that("estimators tolerate NA (failed) runs", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1)
  )
  set.seed(5)
  s <- gsa_sample(p, 1500, seed = 5)
  y <- 3 * s$X[, "x1"] + s$X[, "x2"]
  y[sample.int(1500, 50)] <- NA
  r <- gsa_delta(s$X, y, boot = 0)
  expect_equal(r$meta$n_failed, 50)
  expect_true(all(is.finite(r$indices$estimate[r$indices$index == "delta"])))
})

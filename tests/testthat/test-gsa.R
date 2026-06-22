test_that("two-stage delta identifies an indirect-only (correlation-driven) parameter", {
  p <- gsa_parameters(
    gsa_parameter("k2", "uniform", min = 0, max = 1),
    gsa_parameter("IC50", "uniform", min = 0, max = 1)
  )
  R <- gsa_correlation(p, c("k2", "IC50", 0.9))
  ev <- function_evaluator(function(M) 2 * M[, "k2"], "CT") # only k2 acts
  res <- gsa(p, ev, method = "delta", n = 3000, correlation = R, boot = 40, seed = 1)
  cls <- delta_classification(res)
  cls <- as.data.frame(cls)
  expect_equal(cls$class[cls$parameter == "k2"], "both")
  expect_equal(cls$class[cls$parameter == "IC50"], "indirect-only")
  expect_equal(cls$class[cls$parameter == "dummy"], "non-influential")
})

test_that("gsa() combines several methods", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1)
  )
  ev <- function_evaluator(function(M) 3 * M[, "x1"] + M[, "x2"], "y")
  res <- gsa(p, ev, method = c("delta", "sobol", "morris"), n = 1024, r = 15, boot = 0, seed = 1)
  expect_true(all(c("delta", "sobol", "morris") %in% res$indices$method))
  expect_s3_class(res, "ospgsa_result")
})

test_that("convergence returns an increasing-N table", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = -pi, max = pi),
    gsa_parameter("x2", "uniform", min = -pi, max = pi),
    gsa_parameter("x3", "uniform", min = -pi, max = pi)
  )
  ev <- function_evaluator(function(M) sin(M[, "x1"]) + 7 * sin(M[, "x2"])^2, "y")
  conv <- gsa_convergence(p, ev, method = "sobol", n_seq = c(256, 1024), boot = 0)
  expect_s3_class(conv, "ospgsa_convergence")
  expect_setequal(unique(conv$n), c(256, 1024))
})

test_that("summary and coercion work", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = 0, max = 1),
    gsa_parameter("x2", "uniform", min = 0, max = 1)
  )
  ev <- function_evaluator(function(M) M[, "x1"] + 2 * M[, "x2"], "y")
  res <- gsa(p, ev, method = "delta", n = 1000, boot = 0, seed = 1)
  expect_s3_class(summary(res), "data.table")
  expect_s3_class(as.data.frame(res), "data.frame")
})

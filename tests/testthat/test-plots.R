test_that("plot functions return ggplot objects", {
  p <- gsa_parameters(
    gsa_parameter("x1", "uniform", min = -pi, max = pi),
    gsa_parameter("x2", "uniform", min = -pi, max = pi),
    gsa_parameter("x3", "uniform", min = -pi, max = pi)
  )
  ish <- function(M) sin(M[, "x1"]) + 7 * sin(M[, "x2"])^2 + 0.1 * M[, "x3"]^4 * sin(M[, "x1"])
  ev <- function_evaluator(ish, "y")
  res <- gsa(p, ev, method = c("sobol", "morris"), n = 1024, r = 12, boot = 10, seed = 1)
  expect_s3_class(plot_indices(res, index = "ST"), "ggplot")
  expect_s3_class(plot_sobol(res), "ggplot")
  expect_s3_class(plot_morris(res), "ggplot")
  expect_s3_class(plot(res, type = "sobol"), "ggplot")

  s <- gsa_sample(p, 500, seed = 2)
  expect_s3_class(plot_scatter(s$X, ish(s$X)), "ggplot")
})

test_that("two-stage delta plot builds", {
  p <- gsa_parameters(
    gsa_parameter("k2", "uniform", min = 0, max = 1),
    gsa_parameter("IC50", "uniform", min = 0, max = 1)
  )
  R <- gsa_correlation(p, c("k2", "IC50", 0.9))
  ev <- function_evaluator(function(M) 2 * M[, "k2"], "CT")
  res <- gsa(p, ev, method = "delta", n = 1500, correlation = R, boot = 10, seed = 1)
  expect_s3_class(plot_delta_two_stage(res), "ggplot")
})

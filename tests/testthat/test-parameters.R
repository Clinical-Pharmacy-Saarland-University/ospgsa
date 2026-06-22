test_that("distributions build and quantiles are correct", {
  expect_equal(gsa_parameter("a", "uniform", min = 2, max = 6)$qfun(0.5), 4)
  expect_equal(gsa_parameter("a", "loguniform", min = 1, max = 100)$qfun(0.5), 10)
  expect_equal(gsa_parameter("a", "normal", mean = 3, sd = 2)$qfun(0.5), 3)
  lp <- gsa_parameter("a", "lognormal", median = 5, gsd = 2)
  expect_equal(lp$qfun(0.5), 5, tolerance = 1e-8)
  # lognormal via mean + cv
  lp2 <- gsa_parameter("a", "lognormal", mean = 10, cv = 0.3)
  s <- lp2$qfun(stats::runif(20000))
  expect_equal(mean(s), 10, tolerance = 0.05)
  # truncnorm respects bounds
  tn <- gsa_parameter("a", "truncnorm", mean = 0, sd = 5, lower = -1, upper = 2)
  u <- tn$qfun(c(0.001, 0.999))
  expect_true(u[1] >= -1 && u[2] <= 2)
})

test_that("default scale is log for log-distributed parameters", {
  expect_equal(gsa_parameter("a", "lognormal", median = 1, gsd = 2)$scale, "log")
  expect_equal(gsa_parameter("a", "loguniform", min = 1, max = 10)$scale, "log")
  expect_equal(gsa_parameter("a", "uniform", min = 0, max = 1)$scale, "natural")
})

test_that("invalid distribution arguments error", {
  expect_error(gsa_parameter("a", "uniform", min = 5, max = 1))
  expect_error(gsa_parameter("a", "loguniform", min = -1, max = 10))
  expect_error(gsa_parameter("a", "lognormal", median = 5))
})

test_that("gsa_parameters enforces unique names", {
  expect_error(gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("a", "uniform", min = 0, max = 1)
  ))
})

test_that("correlation builder validates and repairs", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1),
    gsa_parameter("c", "uniform", min = 0, max = 1)
  )
  R <- gsa_correlation(p, c("a", "b", 0.5), c("a", "c", -0.2))
  expect_equal(R["a", "b"], 0.5)
  expect_equal(R["a", "c"], -0.2)
  expect_true(isSymmetric(R))
  expect_error(gsa_correlation(p, c("a", "b", 1.5)))
  expect_warning(gsa_correlation(p, c("a", "b", 0.99), c("a", "c", 0.99), c("b", "c", -0.99)))
})

test_that("dummy parameter is flagged", {
  d <- gsa_dummy()
  expect_true(d$is_dummy)
  p <- gsa_parameters(gsa_parameter("a", "uniform", min = 0, max = 1), gsa_dummy("noise0"))
  expect_equal(.dummy_name(p), "noise0")
})

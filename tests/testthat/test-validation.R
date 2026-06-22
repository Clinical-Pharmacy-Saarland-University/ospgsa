# Input validation, error paths, and less-travelled sampling branches.

test_that("gsa_parameter enforces the log-scale positive-support guard", {
  expect_error(gsa_parameter("a", "normal", mean = 0, sd = 1, scale = "log"), "positive")
  expect_error(gsa_parameter("a", "uniform", min = -1, max = 1, scale = "log"), "positive")
  # explicit natural override is allowed on a log-default distribution
  expect_equal(
    gsa_parameter("a", "lognormal", median = 2, gsd = 2, scale = "natural")$scale,
    "natural"
  )
})

test_that("gsa_parameter rejects out-of-domain distribution arguments", {
  expect_error(gsa_parameter("a", "lognormal", median = 5, gsd = 1), "gsd")
  expect_error(gsa_parameter("a", "lognormal", mean = 0, cv = 0.3), "mean")
  expect_error(gsa_parameter("a", "lognormal", mean = 10, cv = 0), "cv")
  expect_error(gsa_parameter("a", "normal", mean = 0, sd = 0), "sd")
  expect_error(gsa_parameter("a", "truncnorm", mean = 0, sd = 1, lower = 2, upper = 1), "lower")
})

test_that("gsa_parameters rejects empty / wrong-type input and unwraps a list", {
  expect_error(gsa_parameters(), "At least one")
  expect_error(gsa_parameters(1, "x"), "gsa_parameter")
  lst <- list(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  expect_length(gsa_parameters(lst), 2L)
})

test_that("gsa_dummy is uniform[0,1] with an NA path", {
  d <- gsa_dummy("noise")
  expect_equal(d$dist, "uniform")
  expect_equal(unname(d$dist_args$min), 0)
  expect_equal(unname(d$dist_args$max), 1)
  expect_true(is.na(d$path))
  expect_equal(d$qfun(c(0, 0.5, 1)), c(0, 0.5, 1))
})

test_that("gsa_correlation accepts and reorders a full matrix", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1),
    gsa_parameter("c", "uniform", min = 0, max = 1)
  )
  M <- diag(3)
  M[1, 2] <- M[2, 1] <- 0.4
  R <- gsa_correlation(p, M)
  expect_equal(dimnames(R), list(c("a", "b", "c"), c("a", "b", "c")))
  expect_equal(R["a", "b"], 0.4)
  Ms <- matrix(0, 3, 3, dimnames = list(c("c", "a", "b"), c("c", "a", "b")))
  diag(Ms) <- 1
  Ms["a", "b"] <- Ms["b", "a"] <- 0.4
  R2 <- gsa_correlation(p, Ms)
  expect_equal(rownames(R2), c("a", "b", "c"))
  expect_equal(R2["a", "b"], 0.4)
  expect_error(gsa_correlation(p, diag(2)), "wrong dimension")
})

test_that("gsa_correlation errors on non-PD with repair = FALSE and on bad specs", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1),
    gsa_parameter("c", "uniform", min = 0, max = 1)
  )
  expect_error(
    gsa_correlation(p, c("a", "b", 0.99), c("a", "c", 0.99), c("b", "c", -0.99), repair = FALSE),
    "positive semi-definite"
  )
  expect_error(gsa_correlation(p, c("a", "ZZ", 0.5)), "Unknown parameter")
  expect_error(gsa_correlation(p, c("a", "b")), "length-3")
})

test_that("gsa_sample validates n", {
  p <- gsa_parameters(gsa_parameter("a", "uniform", min = 0, max = 1))
  expect_error(gsa_sample(p, 1), "whole number")
  expect_error(gsa_sample(p, 2.5), "whole number")
  expect_error(gsa_sample(p, NA), "whole number")
})

test_that("gsa_sample is reproducible and method = 'random' works", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 10),
    gsa_parameter("b", "uniform", min = 0, max = 10)
  )
  expect_equal(gsa_sample(p, 1000, seed = 11)$X, gsa_sample(p, 1000, seed = 11)$X)
  s1 <- gsa_sample(p, 3000, method = "random", seed = 7)
  expect_equal(s1$X, gsa_sample(p, 3000, method = "random", seed = 7)$X)
  expect_equal(mean(s1$X[, "a"]), 5, tolerance = 0.4)
  R <- gsa_correlation(p, c("a", "b", 0.5))
  sc <- gsa_sample(p, 4000, method = "random", correlation = R, seed = 7)
  expect_equal(unname(cor(sc$X, method = "spearman")["a", "b"]), 0.5, tolerance = 0.08)
})

test_that("gaussian copula with corr_type = 'pearson' targets the latent normal directly", {
  p <- gsa_parameters(
    gsa_parameter("a", "normal", mean = 0, sd = 1),
    gsa_parameter("b", "normal", mean = 0, sd = 1)
  )
  R <- gsa_correlation(p, c("a", "b", 0.6))
  s <- gsa_sample(
    p,
    8000,
    correlation = R,
    correlation_method = "gaussian",
    corr_type = "pearson",
    seed = 1
  )
  expect_equal(unname(cor(s$X)["a", "b"]), 0.6, tolerance = 0.06)
})

test_that("Iman-Conover with k = 1 returns the marginal unchanged", {
  p <- gsa_parameters(gsa_parameter("a", "uniform", min = 0, max = 1))
  R <- gsa_correlation(p)
  s0 <- gsa_sample(p, 500, seed = 3)
  s1 <- gsa_sample(p, 500, correlation = R, correlation_method = "iman-conover", seed = 3)
  expect_equal(s1$X, s0$X)
})

test_that(".to_estimator_inputs log-transforms only log-scale columns", {
  p <- gsa_parameters(
    gsa_parameter("lin", "uniform", min = 1, max = 5), # natural scale
    gsa_parameter("lg", "lognormal", median = 2, gsd = 2) # log scale
  )
  X <- matrix(c(2, 4, exp(1), exp(2)), ncol = 2, dimnames = list(NULL, c("lin", "lg")))
  Z <- .to_estimator_inputs(X, p)
  expect_equal(Z[, "lin"], X[, "lin"])
  expect_equal(Z[, "lg"], log(X[, "lg"]))
})

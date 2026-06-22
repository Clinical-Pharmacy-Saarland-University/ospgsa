test_that("independent sampling reproduces marginals", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 10),
    gsa_parameter("b", "lognormal", median = 5, gsd = 2)
  )
  s <- gsa_sample(p, 5000, seed = 1)
  expect_equal(mean(s$X[, "a"]), 5, tolerance = 0.2)
  expect_equal(median(s$X[, "b"]), 5, tolerance = 0.3)
})

test_that("Iman-Conover induces the target rank correlation, marginals preserved", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "lognormal", median = 2, gsd = 1.5),
    gsa_parameter("c", "normal", mean = 0, sd = 1)
  )
  R <- gsa_correlation(p, c("a", "b", 0.6), c("a", "c", -0.3))
  s0 <- gsa_sample(p, 4000, seed = 1)
  s <- gsa_sample(p, 4000, correlation = R, correlation_method = "iman-conover", seed = 1)
  cc <- cor(s$X, method = "spearman")
  expect_equal(cc["a", "b"], 0.6, tolerance = 0.06)
  expect_equal(cc["a", "c"], -0.3, tolerance = 0.06)
  # marginals (sorted values) identical to independent draw with same seed
  expect_equal(sort(s$X[, "b"]), sort(s0$X[, "b"]))
})

test_that("Gaussian copula reaches target Spearman", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  R <- gsa_correlation(p, c("a", "b", 0.7))
  s <- gsa_sample(p, 5000, correlation = R, correlation_method = "gaussian", seed = 1)
  expect_equal(cor(s$X, method = "spearman")["a", "b"], 0.7, tolerance = 0.06)
})

test_that("LHS covers each margin's strata", {
  p <- gsa_parameters(gsa_parameter("a", "uniform", min = 0, max = 1))
  s <- gsa_sample(p, 100, method = "lhs", seed = 1)
  h <- hist(s$X[, "a"], breaks = seq(0, 1, by = 0.1), plot = FALSE)$counts
  expect_true(all(h >= 8 & h <= 12))
})

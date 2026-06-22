# Custom summary_fn output reducer and its debugger.

test_that(".align_summary maps a named return onto the fixed columns", {
  cn <- c("tmax", "auc")
  expect_equal(.align_summary(c(auc = 2, tmax = 1), cn), c(tmax = 1, auc = 2))
  expect_equal(unname(.align_summary(NULL, cn)), c(NA_real_, NA_real_))
  # missing name -> NA, extra name dropped
  r <- .align_summary(c(tmax = 5, junk = 9), cn)
  expect_equal(r[["tmax"]], 5)
  expect_true(is.na(r[["auc"]]))
  # unnamed return maps positionally when the length matches
  expect_equal(unname(.align_summary(c(1, 2), cn)), c(1, 2))
})

test_that("ospsuite_evaluator summary_fn computes a custom metric end to end", {
  skip_if_not_installed("ospsuite")
  sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
  skip_if(!nzchar(sim))
  p <- gsa_parameters(gsa_parameter(
    "Lipophilicity",
    "normal",
    mean = -0.1,
    sd = 0.3,
    path = "Aciclovir|Lipophilicity"
  ))
  tmax <- function(profiles, sr) {
    d <- profiles[[1]]
    c(tmax = d$time[which.max(d$value)])
  }
  ev <- ospsuite_evaluator(sim, p, summary_fn = tmax)
  Y <- ev(gsa_sample(p, 4, seed = 1)$X)
  expect_equal(colnames(Y), "tmax")
  expect_equal(nrow(Y), 4L)
  expect_true(all(is.finite(Y)))
})

test_that("ospsuite_evaluator rejects a non-function summary_fn", {
  skip_if_not_installed("ospsuite")
  sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
  skip_if(!nzchar(sim))
  p <- gsa_parameters(gsa_parameter(
    "Lipophilicity",
    "normal",
    mean = -0.1,
    sd = 0.3,
    path = "Aciclovir|Lipophilicity"
  ))
  expect_error(ospsuite_evaluator(sim, p, summary_fn = 42), "must be a function")
})

test_that("ospsuite_test_summary returns the profiles, value and results", {
  skip_if_not_installed("ospsuite")
  sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
  skip_if(!nzchar(sim))
  p <- gsa_parameters(gsa_parameter(
    "Lipophilicity",
    "normal",
    mean = -0.1,
    sd = 0.3,
    path = "Aciclovir|Lipophilicity"
  ))
  res <- ospsuite_test_summary(sim, p, function(profiles, sr) {
    c(tmax = profiles[[1]]$time[which.max(profiles[[1]]$value)])
  })
  expect_true(all(c("value", "profiles", "sr") %in% names(res)))
  expect_named(res$value, "tmax")
  expect_s3_class(res$profiles[[1]], "data.frame")
})

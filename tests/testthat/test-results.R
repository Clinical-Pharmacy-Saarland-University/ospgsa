# Result assembly: the dummy noise-floor rule, two-stage classification, combine.

test_that(".finalize_indices flags significance against the dummy CI band", {
  d <- .indices_dt(
    output = "y",
    method = "delta",
    index = "delta",
    stage = "full",
    parameter = c("p_big", "p_small", "dummy"),
    estimate = c(0.40, 0.06, 0.05),
    conf_low = c(0.35, 0.02, 0.02),
    conf_high = c(0.45, 0.09, 0.08)
  )
  out <- .finalize_indices(d, dummy = "dummy")
  sig <- setNames(out$significant, out$parameter)
  expect_true(sig[["p_big"]]) # CI strictly above the dummy band
  expect_false(sig[["p_small"]]) # CI overlaps the dummy band
  expect_false(sig[["dummy"]]) # the dummy is never significant
  expect_equal(out$rank[out$parameter == "p_big"], 1L) # ranked by estimate
})

test_that(".finalize_indices without a dummy uses the CI-excludes-zero rule", {
  d <- .indices_dt(
    output = "y",
    method = "delta",
    index = "delta",
    stage = "full",
    parameter = c("a", "b"),
    estimate = c(0.3, 0.0),
    conf_low = c(0.1, -0.2),
    conf_high = c(0.5, 0.2)
  )
  out <- .finalize_indices(d, dummy = NULL)
  expect_equal(setNames(out$significant, out$parameter), c(a = TRUE, b = FALSE))
})

test_that("delta_classification covers all four classes and errors without two stages", {
  mk <- function(stage, est, sig) {
    .indices_dt(
      output = "y",
      parameter = c("p1", "p2", "p3", "p4"),
      method = "delta",
      index = "delta",
      stage = stage,
      estimate = est,
      significant = sig
    )
  }
  ind <- mk("independent", c(.3, .3, .0, .0), c(TRUE, TRUE, FALSE, FALSE))
  full <- mk("full", c(.3, .0, .3, .0), c(TRUE, FALSE, TRUE, FALSE))
  res <- new_ospgsa_result(data.table::rbindlist(list(ind, full), fill = TRUE))
  cls <- as.data.frame(delta_classification(res))
  cl <- setNames(cls$class, cls$parameter)
  expect_equal(cl[["p1"]], "both")
  expect_equal(cl[["p2"]], "causal")
  expect_equal(cl[["p3"]], "indirect-only")
  expect_equal(cl[["p4"]], "non-influential")
  expect_error(delta_classification(new_ospgsa_result(full)), "two-stage")
})

test_that(".combine_results drops NULLs and errors when empty", {
  expect_error(.combine_results(list(NULL, NULL)), "no results")
})

test_that("summary / as.data.table / print on a result behave", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  ev <- function_evaluator(function(M) 3 * M[, "a"] + M[, "b"], "y")
  res <- gsa_sobol(p, ev, n = 256, boot = 0)
  s <- summary(res)
  expect_s3_class(s, "data.table")
  expect_true(all(c("n_param", "n_significant", "max_index") %in% names(s)))
  dt <- as.data.table(res)
  dt[, estimate := -999]
  expect_false(any(res$indices$estimate == -999)) # a copy, not an alias
  expect_output(print(res), "ospgsa_result")
})

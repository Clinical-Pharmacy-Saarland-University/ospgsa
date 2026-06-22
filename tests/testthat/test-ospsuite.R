# PK-Sim integration smoke test (skipped when ospsuite is unavailable).

test_that("ospsuite_evaluator runs PK-Sim and gsa() works end-to-end", {
  skip_if_not_installed("ospsuite")
  sim_path <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
  skip_if(sim_path == "", "Aciclovir.pkml example not found.")

  p <- gsa_parameters(
    gsa_parameter(
      "Lipophilicity",
      "normal",
      mean = -0.097,
      sd = 0.3,
      path = "Aciclovir|Lipophilicity"
    ),
    gsa_parameter(
      "FractionUnbound",
      "truncnorm",
      mean = 0.85,
      sd = 0.1,
      lower = 0.3,
      upper = 0.999,
      path = "Aciclovir|Fraction unbound (plasma)"
    )
  )
  ev <- ospsuite_evaluator(
    sim_path,
    p,
    pk_parameters = c("AUC_tEnd", "C_max"),
    n_cores = 1,
    silent = TRUE
  )

  Y <- ev(gsa_sample(p, 6, seed = 1)$X)
  expect_equal(nrow(Y), 6L)
  expect_equal(ncol(Y), 2L) # 1 output path x 2 PK parameters (AUC, Cmax)
  expect_false(anyNA(Y))

  res <- gsa(p, ev, method = "delta", n = 24, boot = 5, seed = 2, log_output = TRUE)
  expect_s3_class(res, "ospgsa_result")
  expect_true("dummy" %in% res$indices$parameter)
  expect_true(all(c("Lipophilicity", "FractionUnbound") %in% res$indices$parameter))
})

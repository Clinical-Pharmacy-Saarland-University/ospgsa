# Parameter-setup, design, and output helpers.

test_that(".gsa_anchor_parameter matches direct gsa_parameter construction", {
  a <- .gsa_anchor_parameter("x", "lognormal", nominal = 5, width = 1.5)
  b <- gsa_parameter("x", "lognormal", median = 5, gsd = 1.5)
  expect_equal(a$dist, "lognormal")
  expect_equal(a$dist_args, b$dist_args)

  a <- .gsa_anchor_parameter(
    "x",
    "truncnorm",
    nominal = 0.85,
    width = 0.1,
    lower = 0.3,
    upper = 0.999
  )
  b <- gsa_parameter(
    "x",
    "truncnorm",
    mean = 0.85,
    sd = max(0.1 * 0.85, 1e-6),
    lower = 0.3,
    upper = 0.999
  )
  expect_equal(a$dist_args, b$dist_args)

  a <- .gsa_anchor_parameter("x", "normal", nominal = -0.1, width = 0.3)
  b <- gsa_parameter("x", "normal", mean = -0.1, sd = 0.3)
  expect_equal(a$dist_args, b$dist_args)

  expect_error(.gsa_anchor_parameter("x", "uniform", 0.5, 0.1), "Unsupported")
})

test_that(".resolve_unique resolves one path and errors on 0 / >1", {
  paths <- c("A|kcat", "A|Km", "B|kcat")
  expect_equal(.resolve_unique("x", "^A\\|kcat$", paths), "A|kcat")
  expect_error(.resolve_unique("x", "NOPE", paths), "matched no")
  expect_error(.resolve_unique("x", "kcat", paths), "matched 2")
})

test_that("gsa_correlation accepts a list of triplets", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1),
    gsa_parameter("c", "uniform", min = 0, max = 1)
  )
  R1 <- gsa_correlation(p, c("a", "b", 0.5), c("a", "c", -0.2))
  R2 <- gsa_correlation(p, list(c("a", "b", 0.5), c("a", "c", -0.2)))
  expect_equal(R1, R2)
})

test_that("gsa_two_stage_designs builds an independent + correlated pair", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  R <- gsa_correlation(p, c("a", "b", 0.6))
  d <- gsa_two_stage_designs(p, 500, correlation = R, seed = 1)
  expect_named(d, c("full", "ind"))
  expect_equal(dim(d$full$X), c(500L, 2L))
  expect_false(is.null(d$full$correlation)) # Stage 2 is correlated
  expect_true(is.null(d$ind$correlation)) # Stage 1 is independent
  # distinct seeds: full uses `seed`, ind uses `seed + 1`
  expect_equal(d$full$X, gsa_sample(p, 500, correlation = R, seed = 1)$X)
  expect_equal(d$ind$X, gsa_sample(p, 500, seed = 2)$X)
  # reproducible
  expect_equal(gsa_two_stage_designs(p, 500, correlation = R, seed = 1)$full$X, d$full$X)
})

test_that("gsa_sanitize_positive sets non-positive outputs to NA", {
  Y <- matrix(c(1, -2, 3, 0, 5, 6), ncol = 2)
  s <- gsa_sanitize_positive(Y, quiet = TRUE)
  expect_true(is.na(s[2, 1])) # -2
  expect_true(is.na(s[1, 2])) # 0
  expect_equal(sum(is.na(s)), 2L)

  res <- list(X = matrix(1, 3, 1), Y = matrix(c(1, -1, 2), ncol = 1), n_failed = 0L)
  s2 <- gsa_sanitize_positive(res, quiet = TRUE)
  expect_true(is.na(s2$Y[2, 1]))
  expect_equal(s2$n_failed, 1L)
})

test_that("gsa_save_plots writes PNGs for the available index families", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  f <- function_evaluator(function(M) 2 * M[, "a"] + M[, "b"], "y")

  # Sobol result -> S1/ST plot
  rs <- gsa_sobol(p, f, n = 256, boot = 0)
  dir1 <- tempfile("plots")
  on.exit(unlink(dir1, recursive = TRUE), add = TRUE)
  files1 <- gsa_save_plots(rs, dir = dir1, quiet = TRUE)
  expect_gte(length(files1), 1L)
  expect_true(all(file.exists(file.path(dir1, files1))))
  expect_true(any(grepl("sobol", files1)))

  # Two-stage delta result -> delta_two_stage plot
  R <- gsa_correlation(p, c("a", "b", 0.8))
  d <- gsa_two_stage_designs(p, 1500, correlation = R, seed = 1)
  rd <- gsa_delta_two_stage(d$full$X, f(d$full$X), d$ind$X, f(d$ind$X), boot = 0)
  dir2 <- tempfile("plots")
  on.exit(unlink(dir2, recursive = TRUE), add = TRUE)
  files2 <- gsa_save_plots(rd, dir = dir2, quiet = TRUE)
  expect_true(any(grepl("delta_two_stage", files2)))

  expect_error(gsa_save_plots(list(), dir = tempfile()), "ospgsa_result")
})

test_that("ospsuite_parameters / ospsuite_resolve_paths anchor at model values", {
  skip_if_not_installed("ospsuite")
  sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
  skip_if(!nzchar(sim))
  specs <- data.frame(
    name = c("Lipophilicity", "FractionUnbound"),
    pattern = c("Aciclovir\\|Lipophilicity$", "Aciclovir\\|Fraction unbound \\(plasma\\)$"),
    dist = c("normal", "truncnorm"),
    width = c(0.3, 0.15),
    lower = c(NA, 1e-4),
    upper = c(NA, 0.999),
    stringsAsFactors = FALSE
  )
  params <- ospsuite_parameters(sim, specs, dummy = "dummy")
  expect_s3_class(params, "gsa_parameters")
  expect_length(params, 3L) # 2 specs + dummy
  resolved <- attr(params, "resolved")
  expect_true(all(c("Lipophilicity", "FractionUnbound") %in% resolved$name))
  # the normal marginal must be anchored at the model's nominal value
  lip <- params[["Lipophilicity"]]
  expect_equal(unname(lip$dist_args$mean), resolved$nominal[resolved$name == "Lipophilicity"])

  paths <- ospsuite_resolve_paths(sim, c(lip = "Aciclovir\\|Lipophilicity$"))
  expect_match(unname(paths[["lip"]]), "Lipophilicity")
})

test_that(".path_near_matches surfaces likely intended paths", {
  paths <- c(
    "Organism|Liver|Volume",
    "Aciclovir|Lipophilicity",
    "Aciclovir|Fraction unbound (plasma)"
  )
  # a "|" left unescaped matches nothing as a regex, but the tokens still point home
  near <- .path_near_matches("Aciclovir|Lipophilicity", paths)
  expect_true("Aciclovir|Lipophilicity" %in% near)
  expect_length(.path_near_matches("zzz nothing here", paths), 0L)
})

test_that(".resolve_unique gives an actionable zero-match error", {
  paths <- c("A|kcat", "A|Km")
  err <- expect_error(.resolve_unique("x", "Aciclovir|Lipophilicity", paths))
  expect_match(conditionMessage(err), "ospsuite_parameter_paths")
})

test_that("gsa_delta_two_stage records run metadata", {
  p <- gsa_parameters(
    gsa_parameter("a", "uniform", min = 0, max = 1),
    gsa_parameter("b", "uniform", min = 0, max = 1)
  )
  R <- gsa_correlation(p, c("a", "b", 0.7))
  d <- gsa_two_stage_designs(p, 1000, correlation = R, seed = 1)
  f <- function_evaluator(function(M) 2 * M[, "a"], "y")
  res <- gsa_delta_two_stage(d$full$X, f(d$full$X), d$ind$X, f(d$ind$X), boot = 0, dummy = NULL)
  expect_equal(res$meta$n_runs, 1000L)
  expect_equal(res$meta$n_failed, 0L)
  expect_true(isTRUE(res$meta$correlated))
  expect_identical(res$meta$methods, "delta_two_stage")
})

test_that("ospsuite_parameter_paths lists and searches model paths", {
  skip_if_not_installed("ospsuite")
  sim <- system.file("extdata", "Aciclovir.pkml", package = "ospsuite")
  skip_if(!nzchar(sim))
  all_paths <- ospsuite_parameter_paths(sim)
  expect_type(all_paths, "character")
  expect_gt(length(all_paths), 1L)
  hit <- ospsuite_parameter_paths(sim, "lipophilicity") # ignore_case default
  expect_true(any(grepl("Lipophilicity", hit)))
  tbl <- ospsuite_parameter_paths(sim, "Lipophilicity", fixed = TRUE, value = TRUE)
  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("path", "value") %in% names(tbl)))
})

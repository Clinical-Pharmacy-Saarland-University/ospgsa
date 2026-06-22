# Resumable, crash-safe checkpointed evaluation (gsa_evaluate).

design <- function(n = 23L, seed = 1L) {
  set.seed(seed)
  matrix(stats::runif(n * 2), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
}

model <- function_evaluator(function(M) {
  cbind(a = M[, "x1"] + M[, "x2"], b = M[, "x1"] * M[, "x2"])
})

test_that("checkpointed result matches one-shot evaluation", {
  X <- design()
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  one <- gsa_evaluate(X, model)
  ck <- gsa_evaluate(X, model, checkpoint_dir = dir, block_size = 5, tag = "t", quiet = TRUE)
  expect_equal(ck$Y, one$Y)
  expect_equal(ck$n_blocks, as.integer(ceiling(nrow(X) / 5)))
  expect_equal(ck$n_failed, 0L)
  expect_length(list.files(dir, pattern = "\\.rds$"), ck$n_blocks)
  expect_identical(ck$tag, "t")
})

test_that("a second run loads every block and recomputes nothing", {
  X <- design()
  env <- new.env()
  env$calls <- 0L
  ev <- function(M) {
    env$calls <- env$calls + 1L
    cbind(a = M[, 1], b = M[, 2])
  }
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  r1 <- gsa_evaluate(X, ev, checkpoint_dir = dir, block_size = 5, tag = "t", quiet = TRUE)
  expect_gt(env$calls, 0L)
  env$calls <- 0L
  r2 <- gsa_evaluate(X, ev, checkpoint_dir = dir, block_size = 5, tag = "t", quiet = TRUE)
  expect_equal(env$calls, 0L)
  expect_equal(r2$Y, r1$Y)
})

test_that("partial resume recomputes only the missing blocks", {
  X <- design()
  env <- new.env()
  env$nblocks <- 0L
  ev <- function(M) {
    env$nblocks <- env$nblocks + 1L
    cbind(a = M[, 1] + M[, 2])
  }
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  full <- gsa_evaluate(X, ev, checkpoint_dir = dir, block_size = 5, tag = "t", quiet = TRUE)
  files <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))
  file.remove(files[c(2L, 4L)]) # drop two completed blocks
  env$nblocks <- 0L
  again <- gsa_evaluate(X, ev, checkpoint_dir = dir, block_size = 5, tag = "t", quiet = TRUE)
  expect_equal(again$Y, full$Y)
  expect_equal(env$nblocks, 2L) # exactly the two missing blocks
})

test_that("a stale checkpoint block (wrong row count) is recomputed", {
  X <- design()
  ev <- function_evaluator(function(M) cbind(a = M[, 1]))
  dir <- tempfile("ckpt")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  full <- gsa_evaluate(X, ev, checkpoint_dir = dir, block_size = 5, tag = "t", quiet = TRUE)
  f1 <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))[1L]
  saveRDS(matrix(0, 2L, 1L, dimnames = list(NULL, "a")), f1) # wrong nrow
  again <- gsa_evaluate(X, ev, checkpoint_dir = dir, block_size = 5, tag = "t", quiet = TRUE)
  expect_equal(again$Y, full$Y)
})

test_that("distinct tags keep separate checkpoints in one directory", {
  Xa <- design(seed = 1L)
  Xb <- design(seed = 2L)
  ev <- function_evaluator(function(M) cbind(a = M[, 1] + M[, 2]))
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  ra <- gsa_evaluate(Xa, ev, checkpoint_dir = dir, block_size = 10, tag = "corr", quiet = TRUE)
  rb <- gsa_evaluate(Xb, ev, checkpoint_dir = dir, block_size = 10, tag = "indep", quiet = TRUE)
  expect_true(any(grepl("^corr_", list.files(dir))))
  expect_true(any(grepl("^indep_", list.files(dir))))
  expect_false(isTRUE(all.equal(ra$Y, rb$Y)))
  reload <- gsa_evaluate(Xa, ev, checkpoint_dir = dir, block_size = 10, tag = "corr", quiet = TRUE)
  expect_equal(reload$Y, ra$Y)
})

test_that("changing the design starts a fresh checkpoint set (new key)", {
  X1 <- design(seed = 1L)
  ev <- function_evaluator(function(M) cbind(a = M[, 1]))
  dir <- tempfile("ckpt")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  r1 <- gsa_evaluate(X1, ev, checkpoint_dir = dir, block_size = 10, tag = "t", quiet = TRUE)
  X2 <- X1
  X2[1L, 1L] <- X2[1L, 1L] + 1 # perturb a single value
  r2 <- gsa_evaluate(X2, ev, checkpoint_dir = dir, block_size = 10, tag = "t", quiet = TRUE)
  expect_false(identical(r1$key, r2$key))
})

test_that("crash_skip NA-fills a block whose crash marker survived", {
  X <- design(n = 20L) # 4 blocks of 5
  ev <- function_evaluator(function(M) cbind(a = M[, 1] + M[, 2]))
  dir <- tempfile("ckpt")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  full <- gsa_evaluate(
    X,
    ev,
    checkpoint_dir = dir,
    block_size = 5,
    tag = "t",
    crash_skip = TRUE,
    quiet = TRUE
  )
  files <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))
  file.remove(files[2:4]) # keep block 1 (for column inference)
  file.create(paste0(files[2L], ".crashed")) # block 2 "crashed" last time
  res <- gsa_evaluate(
    X,
    ev,
    checkpoint_dir = dir,
    block_size = 5,
    tag = "t",
    crash_skip = TRUE,
    quiet = TRUE
  )
  expect_equal(nrow(res$Y), nrow(X))
  expect_true(all(is.na(res$Y[6:10, ]))) # block 2 -> NA
  expect_false(anyNA(res$Y[c(1:5, 11:20), ])) # the rest are intact
  expect_equal(res$n_failed, 5L)
  expect_equal(res$Y[c(1:5, 11:20), ], full$Y[c(1:5, 11:20), ])
})

test_that("crash_skip errors when columns cannot yet be inferred", {
  X <- design(n = 20L)
  ev <- function_evaluator(function(M) cbind(a = M[, 1]))
  dir <- tempfile("ckpt")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  key <- .gsa_design_key(X)
  marker <- file.path(dir, sprintf("solo_%s_bs5_block0001.rds.crashed", key))
  file.create(marker) # crash on block 1, nothing completed
  expect_error(
    gsa_evaluate(
      X,
      ev,
      checkpoint_dir = dir,
      block_size = 5,
      tag = "solo",
      crash_skip = TRUE,
      quiet = TRUE
    ),
    "column layout is unknown"
  )
})

test_that("checkpoint args are validated", {
  X <- design()
  expect_error(gsa_evaluate(X, model, checkpoint_dir = tempfile(), block_size = 0), "block_size")
  expect_error(gsa_evaluate(X, "not a function"), "evaluator")
})

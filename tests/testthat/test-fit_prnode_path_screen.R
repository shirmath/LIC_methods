## The strong rule's re-admission gate (estimator plan § Implementation status
## Phase 4b, test R7b; the finding is Z6 of the tuning plan). R7 pins that
## screening is lossless on a 30-point fixture; on the tune_prnode() example's
## 20-point path the re-admission of a wrongly screened group never ran,
## because it waited for a duality gap computed over all groups, which that
## very group kept open. Oracles: the un-screened solver and a cold solve.

.screen_example <- function(){
  key <- "screen_example"
  if(!exists(key, envir = .fx)){
    public <- .fixture_dataset_small()$public
    lambda_max <- lambda_max_prnode(public, mu = c(S = 1, M = 1, D = 1), bool_pooled = TRUE)
    assign(key, list(public = public, lambda = lambda_max*10^seq(0, -2, length.out = 20)), envir = .fx)
  }
  get(key, envir = .fx)
}

test_that("R7b: screening converges everywhere on the example's 20-point path and matches the un-screened solver", {
  .skip_unless_built("fit_prnode_path")
  ex <- .screen_example()
  expect_no_warning(with_screen <- fit_prnode_path(ex$public, lambda = ex$lambda, screen = TRUE))
  expect_no_warning(without <- fit_prnode_path(ex$public, lambda = ex$lambda, screen = FALSE))
  expect_true(all(without$converged))
  expect_true(all(with_screen$converged))
  expect_equal(dim(with_screen$converged), c(4, 20))
  expect_equal(with_screen$objective, without$objective, tolerance = 1e-8)
  for(what in c("s", "d_M", "d_D")){
    expect_equal(with_screen[[what]], without[[what]], tolerance = 1e-6, info = what)
  }
  for(cell in c("S", "M", "D")){
    expect_identical(with_screen$networks[[cell]], without$networks[[cell]], info = cell)
  }
})

test_that("R7b: the warm-started point that stalled ends at the cold solve's objective", {
  .skip_unless_built("fit_prnode_path")
  ex <- .screen_example()
  with_screen <- suppressWarnings(fit_prnode_path(ex$public, lambda = ex$lambda, screen = TRUE))
  cold <- fit_prnode_path(ex$public, lambda = ex$lambda[14], screen = TRUE)
  expect_true(all(cold$converged))
  expect_equal(with_screen$objective[, 14], cold$objective[, 1], tolerance = 1e-8)
  expect_true(all(with_screen$iterations[, 14] < 1e4))
})

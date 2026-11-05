# tune_prnode_spar() (0.9.0): the smoothing level and the penalty chosen
# jointly by the known-variance BIC, and the raw-observation branch of
# fit_prnode() that calls it. The oracle is tune_prnode() on one hand-built
# design per level.

.spar_args <- function(raw){
  list(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time, Y_D = raw$Y_D,
       times_e = raw$times_e, num_lambda = 20, num_decades = 2)
}

## S1 [oracle]: the chosen level is the argmin over the grid of each level's
## minimum criterion, computed here from a design and a tuning per level; the
## returned design and tuning are that level's.
test_that("S1: tune_prnode_spar() returns the level with the smallest minimum criterion and that level's design and tuning", {
  raw <- .fixture_dataset_small()$raw
  grid <- c(0, 0.25, 0.5)
  out <- .with_warnings(do.call(tune_prnode_spar, c(.spar_args(raw), list(spar_grid = grid))))$value
  by_hand <- lapply(grid, function(sp){
    pub <- make_prnode_design(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, q = 1, spar = sp,
                              times_e = raw$times_e)
    list(public = pub, tuning = .with_warnings(tune_prnode(pub, num_lambda = 20, num_decades = 2))$value)
  })
  crit <- vapply(by_hand, function(z) min(z$tuning$criterion), numeric(1))
  expect_equal(unname(out$spar_criterion), crit, tolerance = 1e-10)
  expect_identical(names(out$spar_criterion), as.character(grid))
  best <- which.min(crit)
  expect_equal(out$spar, grid[best])
  expect_equal(out$public$spar, grid[best])
  expect_equal(out$public$x, by_hand[[best]]$public$x)
  expect_equal(out$tuning$pick_idx, by_hand[[best]]$tuning$pick_idx)
  expect_equal(out$tuning$spar, grid[best])
  expect_s3_class(out$tuning, "prnode_tuning")
  # a one-level grid is that level
  one <- .with_warnings(do.call(tune_prnode_spar, c(.spar_args(raw), list(spar_grid = 0.25))))$value
  expect_equal(one$spar, 0.25)
})

## S2: the common weight does not depend on the level (diff2 noise), which is
## what makes the criteria comparable across levels.
test_that("S2: with the default noise estimate the weights are the same at every level", {
  raw <- .fixture_dataset_small()$raw
  w <- vapply(c(0, 0.5, 0.75), function(sp){
    make_prnode_design(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, q = 1, spar = sp,
                       times_e = raw$times_e)$omega[["M"]]
  }, numeric(1))
  expect_equal(w[2], w[1])
  expect_equal(w[3], w[1])
})

## S3: the guards.
test_that("S3: tune_prnode_spar() refuses an empty, duplicated or out-of-range grid", {
  raw <- .fixture_dataset_small()$raw
  args <- .spar_args(raw)
  expect_error(do.call(tune_prnode_spar, c(args, list(spar_grid = numeric(0)))))
  expect_error(do.call(tune_prnode_spar, c(args, list(spar_grid = c(0.5, 0.5)))))
  expect_error(do.call(tune_prnode_spar, c(args, list(spar_grid = c(0.5, 1.5)))))
  expect_error(do.call(tune_prnode_spar, c(args, list(spar_grid = c(0.5, NA)))))
})

## S4: fit_prnode() on raw observations with nothing fixed searches the grid;
## its fit equals the fit at the search's design and tuning. A numeric `spar`
## or a given `lambda` skips the search and builds the design at that level
## (0.5 when only `lambda` is given).
test_that("S4: fit_prnode(raw) searches the spar grid; a numeric spar or a given lambda does not", {
  raw <- .fixture_dataset_small()$raw
  searched <- .with_warnings(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                                        Y_D = raw$Y_D, times_e = raw$times_e, spar_grid = c(0.25, 0.5)))$value
  out <- .with_warnings(tune_prnode_spar(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                         spar_grid = c(0.25, 0.5), times_e = raw$times_e))$value
  expect_equal(searched$tuning$spar, out$spar)
  expect_equal(searched$lambda, out$tuning$lambda)
  expect_equal(searched$d_ode_est, .with_warnings(fit_prnode(out$public, lambda = out$tuning))$value$d_ode_est)

  fixed <- .with_warnings(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                                     Y_D = raw$Y_D, times_e = raw$times_e, spar = 0.25))$value
  expect_equal(fixed$tuning$spar, 0.25)
  numeric_lambda <- .with_warnings(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M,
                                              d_obs_time = raw$d_obs_time, Y_D = raw$Y_D, times_e = raw$times_e,
                                              lambda = fixed$lambda))$value
  expect_null(numeric_lambda$tuning)
  expect_equal(numeric_lambda$path$design$spar, 0.5)
})

## S5: ties in the criterion go to the level listed first, and a daughter
## under 3 frames is refused by the design before any tuning.
test_that("S5: spar ties go to the first listed level; a too-short daughter is refused", {
  raw <- .fixture_dataset_small()$raw
  a <- .with_warnings(tune_prnode_spar(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                       spar_grid = c(0.5, 0.5 + 1e-12), times_e = raw$times_e))$value
  expect_equal(a$spar, 0.5)
  expect_error(tune_prnode_spar(raw$m_obs_time, raw$Y_M, raw$d_obs_time[1:2], raw$Y_D[1:2, , drop = FALSE],
                                times_e = raw$times_e))
})

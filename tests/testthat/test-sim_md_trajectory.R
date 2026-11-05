# Plan section 3: sim_md_trajectory() rewritten around a sim_md_network() truth
# (R/sim_md_trajectory.R). The tests of the retired triplet generator that
# carried this name before Phase 4 (M8) left the package on 2026-10-22.

.md_grid <- .fx_times_coarse
.md_obs <- .fx_obs_time

# p = 3 has six off-diagonal positions: three shared edges leave room for one
# mother-only and one daughter-only edge
.md_truth <- function(seed = 1, num_mother_only = 1, num_daughter_only = 1, frac_unequal = 0.5, p = 3){
  sim_md_network(p = p, q = 1, num_shared = p, num_mother_only = num_mother_only,
                 num_daughter_only = num_daughter_only, frac_unequal = frac_unequal, seed_number = seed)
}

.md_sim <- function(truth, birth_time = .md_obs[37], sigma = rep(0.1, truth$p), seed = 2, ...){
  sim_md_trajectory(truth = truth, obs_time = .md_obs, birth_time = birth_time, times_e = .md_grid,
                    sigma = sigma, seed_number = seed, ...)
}

## M1/M9: the daughter starts from the mother's noiseless state (never from
## the noisy observation), on the grid from her birth on.
test_that("M1: the daughter starts from the mother's noiseless state at birth and is observed from birth on", {
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("sim_md_network")
  truth <- .md_truth()
  birth_time <- .md_obs[37]
  sim <- .md_sim(truth, birth_time = birth_time, sigma = rep(0.5, 3))
  birth_idx <- which(abs(.md_grid - birth_time) < 1e-10)
  expect_identical(sim$daughter$X[1, ], sim$mother$X[birth_idx, ])
  expect_identical(sim$daughter$Y_noiseless[1, ], sim$mother$X[birth_idx, ])
  expect_false(isTRUE(all.equal(sim$daughter$Y_noiseless[1, ], sim$mother$Y[which(.md_obs == birth_time), ])))
  expect_equal(nrow(sim$daughter$X), sum(.md_grid >= birth_time - 1e-12))
  expect_equal(sim$daughter$obs_time[1], birth_time)
  expect_equal(sim$daughter$obs_time, .md_obs[.md_obs >= birth_time])
  expect_equal(sim$daughter$birth_time, birth_time)
  expect_equal(sim$daughter$birth_idx, which(.md_obs == birth_time))
  expect_equal(nrow(sim$daughter$Y), 12)
  expect_equal(nrow(sim$mother$Y), 48)
})

## M2: with d_delta = 0 the daughter continues the mother's trajectory; the
## solver is exact so the tolerance is tight (revision 1's d_same test at 1e-6
## was Euler-limited).
test_that("M2: identical networks continue the mother to 1e-10", {
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("sim_md_network")
  truth <- .md_truth(num_mother_only = 0, num_daughter_only = 0, frac_unequal = 0)
  expect_true(all(truth$d_delta == 0))
  birth_time <- .md_obs[25]
  sim <- .md_sim(truth, birth_time = birth_time, sigma = rep(0, 3))
  birth_idx <- which(abs(.md_grid - birth_time) < 1e-10)
  expect_equal(sim$daughter$X, sim$mother$X[birth_idx:nrow(sim$mother$X), ], tolerance = 1e-10)
})

## M3: different networks diverge continuously from birth.
test_that("M3: with a daughter-only edge the departure from the mother grows from zero at birth", {
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("sim_md_network")
  truth <- .md_truth(num_mother_only = 0, num_daughter_only = 1, frac_unequal = 0)
  birth_time <- .md_obs[25]
  sim <- .md_sim(truth, birth_time = birth_time, sigma = rep(0, 3))
  birth_idx <- which(abs(.md_grid - birth_time) < 1e-10)
  gap <- abs(sim$daughter$X - sim$mother$X[birth_idx:nrow(sim$mother$X), ])
  n <- nrow(gap)
  expect_equal(gap[1, ], rep(0, 3))
  expect_lt(max(gap[1:10, ]), max(gap[(n - 9):n, ]))
  expect_gt(max(gap), 1e-3)
})

## M4: the noise model, checked against its own definition.
test_that("M4: realized noise sd, column independence, one sigma for both cells, and the AR(1) option", {
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("sim_md_network")
  truth <- .md_truth()
  sigma <- c(0.2, 0.5, 1)
  sim <- .md_sim(truth, sigma = sigma, seed = 4)
  realized <- apply(sim$mother$Y - sim$mother$Y_noiseless, 2, stats::sd)
  expect_equal(realized, sigma, tolerance = 0.15)
  # one sigma vector: the daughter gets the same absolute sd (one microscope)
  expect_equal(sim$daughter$sigma, sigma)
  expect_equal(sim$mother$sigma, sigma)

  dense_obs <- seq(0, 1, by = 5e-4)
  sim_dense <- sim_md_trajectory(truth = truth, obs_time = dense_obs, birth_time = dense_obs[1001],
                                 times_e = dense_obs, sigma = sigma, seed_number = 5)
  err <- sim_dense$mother$Y - sim_dense$mother$Y_noiseless
  expect_equal(apply(err, 2, stats::sd), sigma, tolerance = 0.03)
  expect_lt(max(abs(stats::cor(err)[upper.tri(diag(3))])), 0.2)
  err_d <- sim_dense$daughter$Y - sim_dense$daughter$Y_noiseless
  expect_equal(apply(err_d, 2, stats::sd), sigma, tolerance = 0.05)
  # a list(M, D) sigma gives each cell its own
  sim_two <- sim_md_trajectory(truth = truth, obs_time = dense_obs, birth_time = dense_obs[1001],
                               times_e = dense_obs, sigma = list(M = sigma, D = 2*sigma), seed_number = 5)
  expect_equal(apply(sim_two$daughter$Y - sim_two$daughter$Y_noiseless, 2, stats::sd), 2*sigma, tolerance = 0.05)

  sim_ar <- sim_md_trajectory(truth = truth, obs_time = dense_obs, birth_time = dense_obs[1001],
                              times_e = dense_obs, sigma = sigma, noise = "ar1", rho = 0.5, seed_number = 6)
  err_ar <- sim_ar$mother$Y - sim_ar$mother$Y_noiseless
  lag1 <- apply(err_ar, 2, function(e) stats::cor(e[-1], e[-length(e)]))
  expect_equal(lag1, rep(0.5, 3), tolerance = 0.2)
  expect_equal(apply(err_ar, 2, stats::sd), sigma, tolerance = 0.06)
})

## M5: the RNG contract; sigma = 0 still consumes the stream so that a dataset
## built at two noise levels shares its truth draws (revision 1's S4).
test_that("M5: same seed identical; sigma = 0 advances the stream as sigma > 0 does", {
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("sim_md_network")
  truth <- .md_truth()
  s1 <- .md_sim(truth, seed = 9)
  s2 <- .md_sim(truth, seed = 9)
  expect_identical(s1, s2)
  .md_sim(truth, sigma = rep(0, 3), seed = 9)
  after_zero <- .Random.seed
  .md_sim(truth, sigma = rep(1, 3), seed = 9)
  after_one <- .Random.seed
  expect_identical(after_zero, after_one)
  s0 <- .md_sim(truth, sigma = rep(0, 3), seed = 9)
  expect_identical(s0$mother$Y, s0$mother$Y_noiseless)
  expect_identical(s0$mother$X, s1$mother$X)
})

## M6: the endpoints. A daughter born at the last frame has one row (allowed
## here; refused by profile_cell(), C8); one born at the first frame is a
## second mother with her own network.
test_that("M6: birth at the last frame gives one row; birth at the first frame gives a second mother", {
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("sim_md_network")
  .skip_unless_built("sim_ode_trajectory")
  truth <- .md_truth()
  last <- .md_sim(truth, birth_time = .md_obs[48], sigma = rep(0, 3))
  expect_equal(nrow(last$daughter$X), 1)
  expect_equal(nrow(last$daughter$Y), 1)
  expect_identical(c(last$daughter$Y_noiseless), last$mother$X[nrow(last$mother$X), ])
  first <- .md_sim(truth, birth_time = .md_obs[1], sigma = rep(0, 3))
  expect_equal(nrow(first$daughter$Y), 48)
  expect_equal(nrow(first$daughter$X), length(.md_grid))
  second_mother <- sim_ode_trajectory(truth$d_ode_pars, truth$d_intercept, first$mother$X[1, ], .md_grid, q = 1)
  expect_equal(first$daughter$X, second_mother, tolerance = 1e-10, ignore_attr = TRUE)
})

## M7: exact-index lookup (section 6.2 of the derivation notes step 1); map_to() is no longer
## called here (U1).
test_that("M7: an observation time off the grid errors naming it; one on the grid to 1e-12 passes", {
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("sim_md_network")
  truth <- .md_truth()
  bad_obs <- .md_obs
  bad_obs[10] <- bad_obs[10] + 3e-6
  expect_error(sim_md_trajectory(truth = truth, obs_time = bad_obs, birth_time = bad_obs[37], times_e = .md_grid,
                                 sigma = rep(0.1, 3), seed_number = 1),
               as.character(signif(bad_obs[10], 3)))
  ok_obs <- .md_obs + 1e-13
  expect_silent(sim_md_trajectory(truth = truth, obs_time = ok_obs, birth_time = ok_obs[37], times_e = .md_grid,
                                  sigma = rep(0.1, 3), seed_number = 1))
  expect_error(sim_md_trajectory(truth = truth, obs_time = .md_obs, birth_time = 0.5, times_e = .md_grid,
                                 sigma = rep(0.1, 3), seed_number = 1),
               "birth_time")
})

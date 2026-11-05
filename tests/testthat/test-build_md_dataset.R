# Plan sections 4 and 5: the regression-level simulator (Q1-Q4) and the
# dataset pipeline with the unified list(meta, raw, truth, public) format
# (B3-B23), in R/sim_md_dataset.R. B17 and B18 (the inheritance level
# search, a bisection on the retired ratio measure) left with it in package 0.8.0. B1, B2 and C10 compared the retired
# four-block pipeline with the package; they left with it on 2026-10-22.

.bank_obs_time <- seq(0, 1, length.out = 48)
.bank_times_e <- seq(0, 1, by = 1e-5)

# a bank-shaped observed pair for the smoother test (B3): two proteins, 48
# mother frames, the daughter born at frame 37
.bank_shaped_pair <- function(seed = 1){
  set.seed(seed)
  t_M <- .bank_obs_time
  t_D <- t_M[37:48]
  f <- function(tt) cbind(sin(2*pi*tt) + 0.5, cos(3*pi*tt) - 1)
  list(m_obs_time = t_M,
       d_obs_time = t_D,
       birth_time = t_M[37],
       Y_M = f(t_M) + matrix(stats::rnorm(96, sd = 0.1), 48),
       Y_D = f(t_D) + 0.3 + matrix(stats::rnorm(24, sd = 0.1), 12))
}

.source_grade_smoother <- function(env){
  path <- test_path("..", "..", "GRADE", "Sources", "Smoother.R")
  skip_if_not(file.exists(path), "GRADE/Sources/Smoother.R is not reachable from the test directory")
  sys.source(path, envir = env)
  env
}

## Q1: the integrated form holds on the regression-level design: with the true
## trajectories, Z gamma + Psi theta reproduces the noiseless observations up
## to trapezoid error (the derivation notes eq. 3). gamma_c = (x_c(t_c1), b_c).
test_that("Q1: Z gamma + Psi theta equals Y_noiseless to 1e-6 relative in both cells", {
  ds <- .fixture_dataset("regression", num_d = 36)
  for(cell in c("M", "D")){
    hc <- .hand_cell(ds, cell)
    theta <- if(cell == "M") ds$truth$m_ode_pars else ds$truth$d_ode_pars
    intercept <- if(cell == "M") ds$truth$m_intercept else ds$truth$d_intercept
    Y0 <- if(cell == "M") ds$truth$Y_M_noiseless else ds$truth$Y_D_noiseless
    for(j in seq_len(ncol(Y0))){
      gamma <- c(Y0[1, j], intercept[j])
      fitted <- c(hc$Z %*% gamma + hc$Psi %*% theta[j, ])
      expect_equal(fitted, Y0[, j], tolerance = 1e-6, info = paste(cell, j))
    }
    expect_equal(ds$truth[[paste0("gamma_", cell)]][1, ], Y0[1, ], tolerance = 1e-10, info = cell, ignore_attr = TRUE)
    expect_equal(ds$truth[[paste0("gamma_", cell)]][2, ], intercept, tolerance = 1e-10, info = cell, ignore_attr = TRUE)
  }
})

## Q2: profiled OLS on the regression-level design of the noiseless response
## recovers the generator's coefficients: the fixture check that revision 1's
## C5 could only do for the daughter (no rank deficiency now, N4).
test_that("Q2: profiled OLS of the noiseless y recovers theta to 1e-5 and gamma to 1e-6 in both cells", {
  ds <- .fixture_dataset("regression", num_d = 36)
  for(cell in c("M", "D")){
    hc <- .hand_cell(ds, cell)
    theta <- if(cell == "M") ds$truth$m_ode_pars else ds$truth$d_ode_pars
    intercept <- if(cell == "M") ds$truth$m_intercept else ds$truth$d_intercept
    Y0 <- if(cell == "M") ds$truth$Y_M_noiseless else ds$truth$Y_D_noiseless
    y_tilde <- qr.resid(hc$qr_Z, Y0)
    Psi_tilde <- qr.resid(hc$qr_Z, hc$Psi)
    theta_hat <- t(qr.coef(qr(Psi_tilde), y_tilde))
    expect_equal(theta_hat, theta, tolerance = 1e-5, info = cell, ignore_attr = TRUE)
    gamma_hat <- qr.coef(hc$qr_Z, Y0 - hc$Psi %*% t(theta_hat))
    expect_equal(gamma_hat[1, ], Y0[1, ], tolerance = 1e-6, info = cell, ignore_attr = TRUE)
    expect_equal(gamma_hat[2, ], intercept, tolerance = 1e-6, info = cell, ignore_attr = TRUE)
  }
})

## Q3: the two levels differ only in `public`; `meta$level` is the label of
## the level itself and is the one meta element allowed to differ.
test_that("Q3: meta, raw and truth are identical between the ODE and regression levels at one seed", {
  ds_reg <- .fixture_dataset("regression", num_d = 12)
  ds_ode <- .fixture_dataset("ode", num_d = 12)
  expect_identical(ds_reg$raw, ds_ode$raw)
  expect_identical(ds_reg$truth, ds_ode$truth)
  meta_reg <- ds_reg$meta
  meta_ode <- ds_ode$meta
  expect_identical(meta_reg$level, "regression")
  expect_identical(meta_ode$level, "ode")
  meta_reg$level <- meta_ode$level <- NULL
  expect_identical(meta_reg, meta_ode)
})

## Q4/Z16: smoothing is the gap between the levels. The same noisy y fitted by
## profiled OLS on the smoothed design has a larger coefficient error than on
## the true design.
test_that("Q4: public$x differs between levels and the ODE-level OLS error is larger", {
  ds_reg <- .fixture_dataset("regression", num_d = 36)
  ds_ode <- .fixture_dataset("ode", num_d = 36)
  expect_false(isTRUE(all.equal(ds_reg$public$x, ds_ode$public$x)))
  err <- sapply(list(reg = ds_reg, ode = ds_ode), function(ds){
    hc <- .hand_cell(ds, "M")
    y_tilde <- qr.resid(hc$qr_Z, ds$raw$Y_M)
    Psi_tilde <- qr.resid(hc$qr_Z, hc$Psi)
    theta_hat <- t(qr.coef(qr(Psi_tilde), y_tilde))
    sqrt(sum((theta_hat - ds$truth$m_ode_pars)^2))
  })
  expect_gt(err[["ode"]], err[["reg"]])
})

## B3: the package smoother equals GRADE's smoothX(type_smooth = "smoothspline")
## identically, and spar = NULL is the GCV spline that the deleted
## smooth_trajectory() computed (smooth.spline(all.knots = TRUE), no spar).
test_that("B3: smooth_md_trajectory() equals GRADE's smoothX() and spar = NULL is the GCV spline", {
  env <- .source_grade_smoother(new.env())
  pair <- .bank_shaped_pair()
  grade <- env$smoothX(list(cbind(pair$m_obs_time, pair$Y_M)), times_e = .bank_times_e,
                       type_smooth = "smoothspline", type_data = "Replicates")
  ours <- smooth_md_trajectory(pair$m_obs_time, pair$Y_M, .bank_times_e, spar = 0.5)
  expect_identical(ours, grade$Xhat[[1]])
  gcv <- sapply(1:2, function(j){
    stats::predict(stats::smooth.spline(pair$m_obs_time, pair$Y_M[, j], all.knots = TRUE), .bank_times_e)$y
  })
  expect_identical(smooth_md_trajectory(pair$m_obs_time, pair$Y_M, .bank_times_e, spar = NULL), gcv)
  expect_false(isTRUE(all.equal(ours, gcv)))
})

## B4: the daughter grid is what every pipeline script built by hand.
test_that("B4: make_daughter_grid() equals seq(times_e[map_to(bt, times_e)], 1, by = 1e-5) and refuses a decreasing grid", {
  for(bt in c(.bank_times_e[37], 0.5 + 3e-6, 0, 1, .bank_obs_time[37])){
    expect_equal(make_daughter_grid(bt, .bank_times_e),
                 seq(.bank_times_e[map_to(bt, .bank_times_e)], 1, by = 1e-5),
                 tolerance = 0, info = bt)
  }
  expect_length(make_daughter_grid(1, .bank_times_e), 1)
  expect_error(make_daughter_grid(0.5, rev(.bank_times_e)))
})

## B5: the unified format, names exact and alphabetized; public is a pure
## function of raw and the settings.
test_that("B5: meta, raw, truth, public carry exactly the documented names and public rebuilds from raw", {
  ds <- .fixture_dataset("ode", num_d = 12)
  expect_identical(sort(names(ds)), c("meta", "public", "raw", "truth"))
  expect_identical(sort(names(ds$meta)),
                   sort(c("dataset_id", "seed", "num_d", "n_M", "n_D", "q", "level", "noise_level", "noise_level_ratio",
                          "noise_ratio", "noise_stat_m", "noise_stat_d", "n_attempts", "deviation_scale",
                          "birth_time", "length_ratio")))
  expect_identical(sort(names(ds$raw)),
                   sort(c("m_obs_time", "Y_M", "d_obs_time", "Y_D", "birth_time", "times_e", "protein_names")))
  expect_identical(sort(names(ds$truth)),
                   sort(c("s_pars", "m_delta", "d_delta", "m_ode_pars", "d_ode_pars", "m_intercept", "d_intercept",
                          "s_network", "m_network", "d_network", "edge_type", "X_M", "X_D", "Y_M_noiseless",
                          "Y_D_noiseless", "sigma_M", "sigma_D", "gamma_M", "gamma_D", "inheritance",
                          "inheritance_std", "deviation_scale")))
  expect_identical(sort(names(ds$public)),
                   sort(c("x", "y", "index", "M_groups", "D_groups", "omega", "profile", "response_scale",
                          "m_obs_time", "d_obs_time", "times_e", "q", "n_M", "n_D", "p", "protein_names", "spar",
                          "mean_state")))
  expect_identical(colnames(ds$raw$Y_M), ds$raw$protein_names)
  expect_identical(colnames(ds$raw$Y_D), ds$raw$protein_names)
  expect_equal(nrow(ds$raw$Y_D), 12)
  expect_equal(ds$raw$d_obs_time[1], ds$raw$birth_time)
  expect_equal(ds$meta$n_D, 12)
  expect_equal(ds$meta$n_M, 48)
  expect_equal(ds$meta$length_ratio, 12/48)
  rebuilt <- make_prnode_design(m_obs_time = ds$raw$m_obs_time, Y_M = ds$raw$Y_M, d_obs_time = ds$raw$d_obs_time,
                                Y_D = ds$raw$Y_D, times_e = ds$raw$times_e, q = ds$public$q)
  expect_equal(rebuilt, ds$public, tolerance = 1e-12)
})

## B6/B7: the rejection loop's stream contract and acceptance rule.
test_that("B6: draw_md_truth() depends on seed_number alone, reports the attempt, and errors after max_attempts", {
  .skip_unless_built("draw_md_truth", "truth_args")
  .skip_unless_new_pipeline()
  args <- list(num_d = 12, obs_time = .fx_obs_time, truth_args = .fx_truth_args, times_e = .fx_times_coarse)
  set.seed(1)
  stats::runif(3)
  d1 <- do.call(draw_md_truth, c(args, seed_number = 77))
  set.seed(2)
  d2 <- do.call(draw_md_truth, c(args, seed_number = 77))
  expect_identical(d1, d2)
  expect_true(is.numeric(d1$num_attempts) && d1$num_attempts >= 1)
  expect_error(do.call(draw_md_truth, c(args, seed_number = 77, max_attempts = 3, max_growth = 0)),
               "max_attempts")
})

test_that("B7: the accepted draw satisfies the growth bound and max_growth = Inf accepts the first draw", {
  .skip_unless_built("draw_md_truth", "truth_args")
  .skip_unless_new_pipeline()
  args <- list(num_d = 12, obs_time = .fx_obs_time, truth_args = .fx_truth_args, times_e = .fx_times_coarse)
  d <- do.call(draw_md_truth, c(args, seed_number = 78, max_growth = 1.2))
  Y_M <- d$sim$mother$Y_noiseless
  Y_D <- d$sim$daughter$Y_noiseless
  expect_lte(max(abs(Y_D)), 1.2*max(abs(Y_M)))
  d_inf <- do.call(draw_md_truth, c(args, seed_number = 78, max_growth = Inf))
  expect_equal(d_inf$num_attempts, 1)
})

## B8: add_md_noise(), checked against its own definition. The default
## convention changes in Phase 4 to one absolute sd for both cells (one
## microscope); the per-cell convention of the old bank stays reachable
## through bool_per_cell_sd = TRUE.
test_that("B8: realized noise sd matches the ratio, 0 returns the input, a scalar recycles, cells are independent", {
  set.seed(3)
  m_mat <- matrix(stats::rnorm(48*3), 48) + outer(seq_len(48)/48, c(3, 5, 7))
  d_mat <- 2*m_mat[37:48, ] + 1
  ratio <- c(0.1, 0.5, 1.4)
  out <- add_md_noise(m_mat, d_mat, ratio, seed_number = 4)
  m_sd <- apply(m_mat, 2, stats::sd)
  expect_equal(apply(out$Y_M - m_mat, 2, stats::sd), ratio*m_sd, tolerance = 0.15)
  big <- matrix(stats::rnorm(2000*3), 2000) + outer(seq_len(2000)/2000, c(3, 5, 7))
  out_big <- add_md_noise(big, big[1001:2000, ], ratio, seed_number = 5)
  expect_equal(apply(out_big$Y_M - big, 2, stats::sd), ratio*apply(big, 2, stats::sd), tolerance = 0.03)
  zero <- add_md_noise(m_mat, d_mat, 0, seed_number = 6)
  expect_identical(zero$Y_M, m_mat)
  expect_identical(zero$Y_D, d_mat)
  scalar <- add_md_noise(m_mat, d_mat, 0.5, seed_number = 7)
  vector <- add_md_noise(m_mat, d_mat, c(0.5, 0.5, 0.5), seed_number = 7)
  expect_identical(scalar, vector)
  expect_lt(abs(stats::cor(c(out_big$Y_M[1001:2000, ] - big[1001:2000, ]), c(out_big$Y_D - big[1001:2000, ]))), 0.05)
})

test_that("B8: the default gives the daughter the mother's absolute sd; bool_per_cell_sd = TRUE gives each cell its own", {
  set.seed(8)
  m_big <- matrix(stats::rnorm(2000*2), 2000) + outer(seq_len(2000)/2000, c(3, 5))
  d_big <- 4*m_big[1001:2000, ]
  ratio <- c(0.5, 1)
  by_default <- add_md_noise(m_big, d_big, ratio, seed_number = 9)
  expect_equal(apply(by_default$Y_D - d_big, 2, stats::sd), ratio*apply(m_big, 2, stats::sd), tolerance = 0.05)
  .skip_unless_built("add_md_noise", "bool_per_cell_sd")
  per_cell <- add_md_noise(m_big, d_big, ratio, bool_per_cell_sd = TRUE, seed_number = 9)
  expect_equal(apply(per_cell$Y_D - d_big, 2, stats::sd), ratio*apply(d_big, 2, stats::sd), tolerance = 0.05)
  one_sd <- add_md_noise(m_big, d_big, ratio, bool_per_cell_sd = FALSE, seed_number = 9)
  expect_identical(one_sd, by_default)
})

## B9: the per-protein ratio draw.
test_that("B9: draw_noise_ratio() returns the level with jitter 1 and log-uniform ratios otherwise", {
  expect_identical(draw_noise_ratio(6, 0.5, noise_jitter = 1, seed_number = 1), rep(0.5, 6))
  r <- draw_noise_ratio(2000, 0.5, noise_jitter = 1.5, seed_number = 2)
  expect_true(all(r >= 0.5/1.5) && all(r <= 0.5*1.5))
  ks <- stats::ks.test(log(r/0.5), "punif", -log(1.5), log(1.5))
  expect_gt(ks$p.value, 0.01)
})

## B10/B11: the noise statistic; NaN on a constant series is the documented
## behaviour (B11), not a guard.
test_that("B10: noise_stat() is 1 on white noise, small on a smooth series and NaN on a constant", {
  set.seed(10)
  expect_equal(noise_stat(stats::rnorm(1e5)), 1, tolerance = 0.02)
  expect_lt(noise_stat(sin(2*pi*seq(0, 1, length.out = 48))), 0.05)
  expect_true(is.nan(noise_stat(rep(1, 10))))
  expect_error(noise_stat(c(1, 2)))
})

## B12/B13(c): the truth's inheritance is the package's measure on the two
## true systems (package 0.8.0): flow_map_similarity() of m_ode_pars and
## d_ode_pars at the horizon of one sampling interval, the convention of
## inheritance_flow(). It is a function of the coefficient matrices alone, so
## it must not move with the noise, the daughter's length or the level, and it
## is 1 exactly when the daughter's system is the mother's.
test_that("B12: truth$inheritance is the flow-map similarity of the true systems at h = one frame, in [0, 1], 1 iff the systems coincide", {
  h <- stats::median(diff(.fx_obs_time))
  for(num_d in c(12, 36)){
    ds <- .fixture_dataset("regression", num_d = num_d)
    expect_identical(ds$truth$inheritance,
                     flow_map_similarity(ds$truth$m_ode_pars, ds$truth$d_ode_pars, h = h), info = num_d)
    expect_gte(ds$truth$inheritance, 0)
    expect_lt(ds$truth$inheritance, 1)
    expect_equal(ds$truth$deviation_scale, 1)
  }
  # no private edges and no unequal shared edge: both deltas are zero and the daughter is a clone
  same <- build_md_dataset(seed = 5, num_d = 12, noise_level_ratio = 0.5,
                           truth_args = list(p = 4, q = 1, num_shared = 4, num_mother_only = 0,
                                             num_daughter_only = 0, frac_unequal = 0),
                           obs_time = .fx_obs_time, times_e = .fx_times_coarse, level = "regression")
  expect_true(all(same$truth$d_delta == 0))
  expect_true(all(same$truth$m_delta == 0))
  expect_identical(same$truth$d_ode_pars, same$truth$m_ode_pars)
  expect_identical(same$truth$inheritance, 1)
})

test_that("B12b: truth$inheritance_std is the similarity of S^-1 A S with S the pooled sd of the two cells' noiseless states, and 1 for a clone", {
  h <- stats::median(diff(.fx_obs_time))
  ds <- .fixture_dataset("regression", num_d = 12)
  truth <- ds$truth
  s_vec <- sapply(seq_len(ncol(truth$Y_M_noiseless)),
                  function(j) stats::sd(c(truth$Y_M_noiseless[, j], truth$Y_D_noiseless[, j])))
  to_std <- function(A) diag(1/s_vec) %*% A %*% diag(s_vec)
  expect_equal(truth$inheritance_std,
               flow_map_similarity(to_std(truth$m_ode_pars), to_std(truth$d_ode_pars), h = h),
               tolerance = 1e-12)
  expect_gte(truth$inheritance_std, 0)
  expect_lte(truth$inheritance_std, 1)
  same <- build_md_dataset(seed = 5, num_d = 12, noise_level_ratio = 0.5,
                           truth_args = list(p = 4, q = 1, num_shared = 4, num_mother_only = 0,
                                             num_daughter_only = 0, frac_unequal = 0),
                           obs_time = .fx_obs_time, times_e = .fx_times_coarse, level = "regression")
  expect_equal(same$truth$inheritance_std, 1)
})

test_that("B12: truth$inheritance depends on the truth alone at q = 1: not on the noise, the level or num_d; at q = 3 each field is linearized at its own cell's mean noiseless state", {
  base <- .fixture_dataset("regression", num_d = 12)
  expect_identical(.fixture_dataset("regression", num_d = 36)$truth$inheritance, base$truth$inheritance)
  expect_identical(.fixture_dataset("ode", num_d = 12)$truth$inheritance, base$truth$inheritance)
  expect_identical(.fixture_dataset("regression", num_d = 12, noise_level_ratio = 1)$truth$inheritance,
                   base$truth$inheritance)
  # a given truth under two noise seeds and two noise levels
  truth <- sim_md_network(p = 6, q = 1, num_shared = 6, num_mother_only = 2, num_daughter_only = 2,
                          frac_unequal = 0.5, seed_number = 7)
  grid <- expand.grid(seed = c(305, 306), noise = c(0.1, 1))
  vals <- sapply(seq_len(nrow(grid)), function(r){
    build_md_dataset(seed = grid$seed[r], num_d = 12, noise_level_ratio = grid$noise[r], truth = truth,
                     obs_time = .fx_obs_time, times_e = .fx_times_coarse, level = "regression")$truth$inheritance
  })
  expect_true(all(vals == vals[1]))
  expect_identical(vals[1], flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars,
                                                h = stats::median(diff(.fx_obs_time))))
  # a cubic truth: each cell's Jacobian at the mean of its own noiseless
  # frames (package 0.13.0; until 0.12.0 both at the mother's), in raw units
  # for `inheritance` and in the pooled-sd units for `inheritance_std`
  ds3 <- build_md_dataset(seed = 3, num_d = 18, noise_level_ratio = 0.1,
                          truth_args = list(p = 4, q = 3, num_shared = 4, num_mother_only = 1,
                                            num_daughter_only = 1),
                          obs_time = seq(0, 1, length.out = 36), steps_per_frame = 100, q = 3,
                          level = "regression")
  h3 <- stats::median(diff(seq(0, 1, length.out = 36)))
  own3 <- list(M = colMeans(ds3$truth$Y_M_noiseless), D = colMeans(ds3$truth$Y_D_noiseless))
  expect_identical(ds3$truth$inheritance,
                   flow_map_similarity(ds3$truth$m_ode_pars, ds3$truth$d_ode_pars, h = h3, q = 3, x_ref = own3))
  expect_false(isTRUE(all.equal(ds3$truth$inheritance,
                                flow_map_similarity(ds3$truth$m_ode_pars, ds3$truth$d_ode_pars, h = h3, q = 3,
                                                    x_ref = list(M = own3$M, D = own3$M)))))
  expect_gte(ds3$truth$inheritance, 0)
  expect_lte(ds3$truth$inheritance, 1)
  s3 <- sapply(1:4, function(j) stats::sd(c(ds3$truth$Y_M_noiseless[, j], ds3$truth$Y_D_noiseless[, j])))
  to_std3 <- function(A) diag(1/s3) %*% A %*% diag(s3)
  expect_equal(ds3$truth$inheritance_std,
               flow_map_similarity(to_std3(LICmethods:::.jacobian_at(ds3$truth$m_ode_pars, 3, own3$M)),
                                   to_std3(LICmethods:::.jacobian_at(ds3$truth$d_ode_pars, 3, own3$D)), h = h3),
               tolerance = 1e-12)
})

## B14: input guards of the builder and the design, each a stopifnot/stop
## with a message.
test_that("B14: num_d outside [3, n_M] and malformed observed pairs are refused", {
  .skip_unless_new_pipeline()
  args <- list(seed = 1, noise_level_ratio = 0.1, truth_args = .fx_truth_args, obs_time = .fx_obs_time,
               times_e = .fx_times_coarse, level = "regression")
  expect_error(do.call(build_md_dataset, c(args, num_d = 2)), "num_d")
  expect_error(do.call(build_md_dataset, c(args, num_d = 49)), "num_d")
  ds <- .fixture_dataset("ode", num_d = 12)
  raw <- ds$raw
  call_design <- function(...){
    a <- utils::modifyList(list(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                                Y_D = raw$Y_D, times_e = raw$times_e, q = 1), list(...))
    do.call(make_prnode_design, a)
  }
  expect_error(call_design(Y_M = raw$Y_M[-1, ]), "m_obs_time")
  expect_error(call_design(m_obs_time = rev(raw$m_obs_time)))
  expect_error(call_design(Y_D = raw$Y_D[, 1:5]), "ncol")
  expect_error(call_design(d_obs_time = raw$d_obs_time + 3e-6))
  # every daughter frame on the grid (half a frame interval is 1000 grid
  # steps) but none of them a mother frame: refused naming d_obs_time
  half_frame <- (raw$m_obs_time[2] - raw$m_obs_time[1])/2
  shifted <- raw$d_obs_time[-length(raw$d_obs_time)] + half_frame
  expect_true(all(sapply(shifted, function(tt) any(abs(raw$times_e - tt) < 1e-10))))
  expect_false(any(sapply(shifted, function(tt) any(abs(raw$m_obs_time - tt) < 1e-10))))
  expect_error(call_design(d_obs_time = shifted, Y_D = raw$Y_D[-1, ]), "d_obs_time")
})

## B16: the banks are read back from disk; a round trip changes nothing.
test_that("B16: a dataset survives saveRDS/readRDS and fits identically", {
  .skip_unless_built("fit_prnode", "public")
  ds <- .fixture_dataset("ode", num_d = 12)
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(ds, path)
  back <- readRDS(path)
  expect_identical(back, ds)
  lambda <- .lambda_path(ds$public, num_lambda = 3)[2]
  fit_mem <- fit_prnode(ds$public, lambda = lambda)
  fit_disk <- fit_prnode(back$public, lambda = lambda)
  expect_identical(fit_disk$m_ode_est, fit_mem$m_ode_est)
  expect_identical(fit_disk$d_ode_est, fit_mem$d_ode_est)
})

## B19: the input guards of the daughter's deviation. The level search of
## packages 0.4.0 to 0.7.0 (`inheritance = target`, a bisection on the retired
## ratio measure) is gone in 0.8.0: the flow-map similarity is bounded away from 1
## by the mother's own deviation and is not monotone in the scale, so a dataset
## reads its similarity off instead of targeting it. The old argument must be
## refused, not silently ignored.
test_that("B19: deviation_scale takes one finite non-negative number and the level-search arguments are refused", {
  .skip_unless_new_pipeline()
  base <- list(seed = 304, num_d = 12, noise_level_ratio = 0.1, obs_time = .fx_obs_time,
               times_e = .fx_times_coarse, level = "regression", truth_args = list(p = 6))
  ds <- .fixture_dataset("regression", num_d = 12)
  expect_identical(ds$truth$deviation_scale, 1)
  expect_identical(ds$meta$deviation_scale, 1)
  for(arg in c("inheritance", "inheritance_tol", "inheritance_max_iter")){
    extra <- list(0.5)
    names(extra) <- arg
    expect_error(do.call(build_md_dataset, c(base, extra)), "unused argument", info = arg)
  }
  for(bad in list(-1, NA_real_, Inf, c(1, 2))){
    expect_error(do.call(build_md_dataset, c(base, list(deviation_scale = bad))),
                 "deviation_scale", info = paste(bad, collapse = ","))
  }
  expect_error(draw_md_truth(num_d = 12, obs_time = .fx_obs_time, truth_args = list(p = 6), deviation_scale = -1,
                             times_e = .fx_times_coarse), "deviation_scale")
})

## B22/B23: setting #3, a fixed multiplier on the drawn
## deviation, is the tuning bank's inheritance level: one helper behind
## draw_md_truth() and the given-truth path, so the two agree, and a larger
## deviation reads a lower similarity (for this seed; the similarity is not
## monotone in the scale in general, B23).
test_that("B22: deviation_scale multiplies d_delta, keeps the adjacency and lowers the similarity, in the draw and on a given truth", {
  .skip_unless_new_pipeline()
  base <- list(seed = 307, num_d = 24, noise_level_ratio = 0.3, truth_args = list(p = 6), obs_time = .fx_obs_time,
               times_e = .fx_times_coarse, level = "regression")
  as_drawn <- do.call(build_md_dataset, base)
  doubled <- do.call(build_md_dataset, c(base, list(deviation_scale = 2)))
  expect_equal(doubled$truth$d_delta, 2*as_drawn$truth$d_delta)
  expect_equal(doubled$truth$d_ode_pars, as_drawn$truth$s_pars + 2*as_drawn$truth$d_delta)
  expect_identical(doubled$truth$m_ode_pars, as_drawn$truth$m_ode_pars)
  expect_identical(doubled$truth$d_network, as_drawn$truth$d_network)
  expect_identical(doubled$truth$edge_type, as_drawn$truth$edge_type)
  expect_identical(doubled$truth$deviation_scale, 2)
  expect_identical(doubled$meta$deviation_scale, 2)
  expect_lt(doubled$truth$inheritance, as_drawn$truth$inheritance)
  # the mother is the same cell; the daughter's trajectory is not
  expect_identical(doubled$raw$Y_M, as_drawn$raw$Y_M)
  expect_false(isTRUE(all.equal(doubled$raw$Y_D, as_drawn$raw$Y_D)))
  # the draw itself: draw_md_truth() scales before the growth check and the stream is untouched by the scaling
  d1 <- draw_md_truth(num_d = 24, obs_time = .fx_obs_time, truth_args = list(p = 6), times_e = .fx_times_coarse,
                      seed_number = 307)
  d2 <- draw_md_truth(num_d = 24, obs_time = .fx_obs_time, truth_args = list(p = 6), deviation_scale = 2,
                      times_e = .fx_times_coarse, seed_number = 307)
  expect_equal(d2$sim$truth$d_delta, 2*d1$sim$truth$d_delta)
  expect_identical(d2$sim$truth$m_ode_pars, d1$sim$truth$m_ode_pars)
  expect_identical(d2$num_attempts, d1$num_attempts)
  # a given truth is scaled the same way
  given <- do.call(build_md_dataset, c(base[names(base) != "truth_args"],
                                       list(truth = d1$sim$truth, deviation_scale = 2)))
  expect_equal(given$truth$d_delta, 2*d1$sim$truth$d_delta)
  expect_identical(given$truth$deviation_scale, 2)
})

test_that("B23: deviation_scale = 0 gives the daughter the scaffold with relabelled edge types, a clone only of a mother without deviation", {
  .skip_unless_new_pipeline()
  base <- list(seed = 308, num_d = 12, noise_level_ratio = 0.5, truth_args = list(p = 6), obs_time = .fx_obs_time,
               times_e = .fx_times_coarse, level = "regression")
  none <- do.call(build_md_dataset, c(base, list(deviation_scale = 0)))
  expect_true(all(none$truth$d_delta == 0))
  expect_identical(none$truth$d_ode_pars, none$truth$s_pars)
  expect_false(any(none$truth$edge_type == "daughter_only"))
  expect_identical(none$truth$d_network, none$truth$s_network)
  expect_identical(none$truth$deviation_scale, 0)
  expect_identical(none$meta$deviation_scale, 0)
  # the mother keeps her own deviation, so the scaffold daughter is not her clone
  expect_identical(none$truth$inheritance,
                   flow_map_similarity(none$truth$m_ode_pars, none$truth$s_pars,
                                       h = stats::median(diff(.fx_obs_time))))
  expect_true(any(none$truth$m_delta != 0))
  expect_lt(none$truth$inheritance, 1)
  # the similarity is not monotone in the scale: it falls toward the chance level and turns back
  # (the reason the level search was retired), so no test asserts monotonicity beyond B22's pair.
  # One truth, drawn once and rescaled through `truth = `: draw_md_truth() would redraw the whole
  # pair if a scaled daughter tripped max_growth, and then the scales would not share a truth
  truth_once <- draw_md_truth(num_d = 12, obs_time = .fx_obs_time, truth_args = list(p = 6),
                              times_e = .fx_times_coarse, seed_number = 308)$sim$truth
  vals <- sapply(c(0, 1, 4, 16), function(scale){
    do.call(build_md_dataset, c(base[names(base) != "truth_args"],
                                list(truth = truth_once, deviation_scale = scale)))$truth$inheritance
  })
  expect_identical(vals[1], none$truth$inheritance)
  expect_true(all(vals >= 0 & vals <= 1))
  expect_lt(vals[3], vals[1])
})

## B20 (plan P1(a)): the documented defaults run together. With times_e =
## NULL the grid is derived from the mother's frames at steps_per_frame
## steps per frame interval, so every frame is an exact grid point (C2, M7),
## and public is the design built on that grid. make_prnode_design() derives
## the same grid on real-shaped frames, and unequally spaced frames still
## error, since the derived grid is equally spaced.
.max_off_grid <- function(frames, grid) max(sapply(frames, function(t) min(abs(grid - t))))
test_that("B20: times_e = NULL derives a frame-aligned grid from the mother's frames, and public is built on it", {
  .skip_unless_new_pipeline()
  base <- list(seed = 31, num_d = 12, noise_level_ratio = 0.1, truth_args = .fx_truth_args, obs_time = .fx_obs_time,
               level = "regression")
  ds <- do.call(build_md_dataset, base)
  n_frames <- length(.fx_obs_time)
  expect_length(ds$raw$times_e, (n_frames - 1)*2000 + 1)
  expect_equal(range(ds$raw$times_e), range(.fx_obs_time), tolerance = 1e-12)
  expect_lt(.max_off_grid(ds$raw$m_obs_time, ds$raw$times_e), 1e-12)
  expect_lt(.max_off_grid(ds$raw$d_obs_time, ds$raw$times_e), 1e-12)
  # the dataset is at the regression level, so the rebuild hands the design
  # the true trajectories; without x_list it would smooth the observations
  # instead, which Q4 asserts gives a different design
  rebuilt <- make_prnode_design(m_obs_time = ds$raw$m_obs_time, Y_M = ds$raw$Y_M, d_obs_time = ds$raw$d_obs_time,
                                Y_D = ds$raw$Y_D, times_e = ds$raw$times_e, q = ds$public$q,
                                x_list = list(M = ds$truth$X_M, D = ds$truth$X_D))
  expect_equal(ds$public, rebuilt, tolerance = 1e-12)
  coarse <- do.call(build_md_dataset, c(base, list(steps_per_frame = 200)))
  expect_length(coarse$raw$times_e, (n_frames - 1)*200 + 1)
  # the design's own derivation on real-shaped frames: 48 frames every 12 minutes
  frames_min <- seq(0, by = 12, length.out = n_frames)
  d_idx <- match(ds$raw$d_obs_time, ds$raw$m_obs_time)
  expect_false(anyNA(d_idx))
  real_shaped <- make_prnode_design(m_obs_time = frames_min, Y_M = ds$raw$Y_M, d_obs_time = frames_min[d_idx],
                                    Y_D = ds$raw$Y_D, q = 1)
  expect_length(real_shaped$times_e, (n_frames - 1)*2000 + 1)
  expect_equal(real_shaped$times_e[1 + 2000*(seq_len(n_frames) - 1)], frames_min, tolerance = 1e-12)
  # one stretched last interval puts the second frame (t = 12) off the derived grid
  uneven <- frames_min
  uneven[n_frames] <- uneven[n_frames] + 5
  expect_error(make_prnode_design(m_obs_time = uneven, Y_M = ds$raw$Y_M, d_obs_time = frames_min[d_idx],
                                  Y_D = ds$raw$Y_D, q = 1), "12")
})

## B21: a truth drawn by hand with
## sim_md_network() can be handed to build_md_dataset(), so one truth runs
## through a whole document (the vignette hands its hand-drawn truth
## to the dataset this way). The given networks are carried through unchanged,
## `seed` still governs the noise alone, `truth_args` is refused beside it,
## the growth bound errors instead of redrawing (there is nothing to redraw),
## and `deviation_scale` still rescales the given daughter deviation.
test_that("B21: build_md_dataset(truth = ) carries the given networks through and refuses truth_args", {
  .skip_unless_new_pipeline()
  truth <- sim_md_network(p = 6, q = 1, num_shared = 6, num_mother_only = 2, num_daughter_only = 2,
                          frac_unequal = 0.5, seed_number = 7)
  base <- list(seed = 305, num_d = 12, noise_level_ratio = 0.1, obs_time = .fx_obs_time,
               times_e = .fx_times_coarse, level = "regression")
  ds <- do.call(build_md_dataset, c(base, list(truth = truth)))
  for(nm in c("s_pars", "m_delta", "d_delta", "m_ode_pars", "d_ode_pars", "m_intercept", "d_intercept",
              "s_network", "m_network", "d_network", "edge_type")){
    expect_identical(ds$truth[[nm]], truth[[nm]], info = nm)
  }
  expect_true(is.na(ds$meta$n_attempts))
  expect_identical(sort(names(ds$meta)), sort(names(.fixture_dataset("regression", num_d = 12)$meta)))
  # the seed governs the noise: the same seed repeats the data, another seed changes only the noise
  ds_again <- do.call(build_md_dataset, c(base, list(truth = truth)))
  expect_identical(ds_again$raw, ds$raw)
  ds_other <- do.call(build_md_dataset, c(base[names(base) != "seed"], list(seed = 306, truth = truth)))
  expect_false(identical(ds_other$raw$Y_M, ds$raw$Y_M))
  # the mother's initial state is drawn from the stream too (sim_md_trajectory() draws rnorm(p) when
  # init_vals is NULL), so another seed moves the noiseless trajectories as well as the noise
  expect_false(identical(ds_other$truth$Y_M_noiseless[1, ], ds$truth$Y_M_noiseless[1, ]))
  expect_identical(ds_other$truth$m_ode_pars, truth$m_ode_pars)
  # refused inputs
  expect_error(do.call(build_md_dataset, c(base, list(truth = truth, truth_args = list(p = 6)))), "truth_args")
  expect_error(do.call(build_md_dataset, c(base, list(truth = truth$m_ode_pars))), "truth")
  expect_error(do.call(build_md_dataset, c(base, list(truth = truth, max_growth = 0))), "max_growth")
  # deviation_scale on a given truth: the mother is untouched, the daughter's deviation is rescaled
  ds_half <- do.call(build_md_dataset, c(base, list(truth = truth, deviation_scale = 0.5)))
  expect_identical(ds_half$truth$m_ode_pars, truth$m_ode_pars)
  expect_equal(ds_half$truth$d_delta, 0.5*truth$d_delta)
  expect_identical(ds_half$truth$inheritance,
                   flow_map_similarity(truth$m_ode_pars, truth$s_pars + 0.5*truth$d_delta,
                                       h = stats::median(diff(.fx_obs_time))))
})

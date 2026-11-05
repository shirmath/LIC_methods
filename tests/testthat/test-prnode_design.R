# Plan sections 6 and 7: make_integrated_basis(), profile_cell() (C1-C10) and
# the two-block design make_prnode_design() (D1-D8), in
# R/prnode_design.R. C10 reads the frozen four-block get_cov_mat() from
# helper-four-block.R.

# a synthetic observed pair with its true trajectories on the aligned grid, so
# that the regression route (x_list) needs no smoothing and n_D = 3 is allowed
.design_pair <- function(p = 6, num_d = 12, seed = 1, q = 1, grid = .fx_times_coarse, sigma = 0.1){
  .skip_unless_built("sim_md_network")
  .skip_unless_built("sim_ode_trajectory")
  # p(p - 1) off-diagonal positions: p = 2 has room for two edges only
  counts <- if(p == 1) c(0, 0, 0) else if(p == 2) c(1, 1, 0) else c(p, 1, 1)
  truth <- sim_md_network(p = p, q = q, num_shared = counts[1], num_mother_only = counts[2],
                          num_daughter_only = counts[3], seed_number = seed)
  obs_time <- .fx_obs_time
  birth_time <- obs_time[48 - num_d + 1]
  set.seed(seed)
  x0 <- stats::rnorm(p)
  X_M <- sim_ode_trajectory(truth$m_ode_pars, truth$m_intercept, x0, grid, q = q)
  d_grid <- grid[grid >= birth_time - 1e-12]
  birth_idx <- which(abs(grid - birth_time) < 1e-10)
  X_D <- sim_ode_trajectory(truth$d_ode_pars, truth$d_intercept, X_M[birth_idx, ], d_grid, q = q)
  m_idx <- sapply(obs_time, function(tt) which(abs(grid - tt) < 1e-10))
  d_obs <- obs_time[obs_time >= birth_time]
  d_idx <- sapply(d_obs, function(tt) which(abs(d_grid - tt) < 1e-10))
  Y_M <- X_M[m_idx, , drop = FALSE] + matrix(stats::rnorm(48*p, sd = sigma), 48)
  Y_D <- X_D[d_idx, , drop = FALSE] + matrix(stats::rnorm(num_d*p, sd = sigma), num_d)
  list(truth = truth, m_obs_time = obs_time, d_obs_time = d_obs, Y_M = Y_M, Y_D = Y_D, times_e = grid,
       d_times_e = d_grid, X_M = X_M, X_D = X_D, birth_time = birth_time, p = p, q = q)
}

.design_from_pair <- function(pair, q = pair$q, bool_regression = TRUE, ...){
  .skip_unless_built("make_prnode_design")
  make_prnode_design(m_obs_time = pair$m_obs_time, Y_M = pair$Y_M, d_obs_time = pair$d_obs_time, Y_D = pair$Y_D,
                     times_e = pair$times_e, q = q,
                     x_list = if(bool_regression) list(M = pair$X_M, D = pair$X_D) else NULL, ...)
}

## C1: closed-form integrals of the cumulative trapezoid rule; the first row
## is exactly zero by construction (section 6.2 of the derivation notes step 1).
test_that("C1: the integrated basis matches closed forms and starts at exactly zero", {
  .skip_unless_built("make_integrated_basis")
  grid <- seq(0, 1, by = 1e-5)
  obs <- grid[c(1, 10001, 25001, 50001, 100001)]
  Psi_lin <- make_integrated_basis(matrix(grid, ncol = 1), grid, obs, q = 2)
  expect_equal(dim(Psi_lin), c(5, 2))
  expect_identical(Psi_lin[1, ], c(0, 0))
  expect_equal(Psi_lin[, 1], obs^2/2, tolerance = 1e-9)
  expect_equal(Psi_lin[, 2], obs^3/3, tolerance = 1e-9)
  Psi_sin <- make_integrated_basis(matrix(sin(2*pi*grid), ncol = 1), grid, obs, q = 1)
  expect_equal(Psi_sin[, 1], (1 - cos(2*pi*obs))/(2*pi), tolerance = 1e-8)
  expect_identical(Psi_sin[1, 1], 0)
})

## C2: exact-index lookup replaces the nearest-neighbour map_to() of the old
## design (revision 1's C6/U1); floating seq() arithmetic must still pass.
test_that("C2: a frame off the grid by 1e-6 errors naming it; one on the grid to 1e-12 passes", {
  .skip_unless_built("make_integrated_basis")
  grid <- seq(0, 1, by = 1e-4)
  x_mat <- matrix(grid, ncol = 1)
  obs_ok <- c(0, 0.25, 0.5 + 1e-13, 1)
  expect_silent(make_integrated_basis(x_mat, grid, obs_ok, q = 1))
  obs_bad <- c(0, 0.25, 0.5 + 1e-6, 1)
  expect_error(make_integrated_basis(x_mat, grid, obs_bad, q = 1), "0.5")
})

## C3: the daughter's integral starts at her own first frame: with a daughter
## trajectory equal to the mother's on the overlap, Psi_D = Psi_M - Psi_M(t_b).
test_that("C3: Psi_D(t_i) equals Psi_M(t_i) - Psi_M(t_b) to 1e-12 on a common trajectory", {
  .skip_unless_built("make_integrated_basis")
  grid <- .fx_times_coarse
  obs <- .fx_obs_time
  x_M <- cbind(sin(2*pi*grid), exp(-grid))
  birth_time <- obs[30]
  d_grid <- grid[grid >= birth_time - 1e-12]
  x_D <- x_M[grid >= birth_time - 1e-12, ]
  d_obs <- obs[obs >= birth_time]
  Psi_M <- make_integrated_basis(x_M, grid, obs, q = 2)
  Psi_D <- make_integrated_basis(x_D, d_grid, d_obs, q = 2)
  Psi_M_at_birth <- Psi_M[30, ]
  expect_equal(Psi_D, sweep(Psi_M[30:48, ], 2, Psi_M_at_birth), tolerance = 1e-12)
  expect_identical(Psi_D[1, ], rep(0, 4))
})

## C4 (T2): both cells' bases share Delta and psi. The signature derives the
## daughter's grid from the mother's, so differing steps can only enter
## through x_list; a daughter matrix not on that grid is refused.
test_that("C4: an x_list trajectory that is not on the cell's grid is refused", {
  pair <- .design_pair(p = 2, num_d = 12)
  expect_silent(.design_from_pair(pair))
  expect_error(make_prnode_design(m_obs_time = pair$m_obs_time, Y_M = pair$Y_M, d_obs_time = pair$d_obs_time,
                                  Y_D = pair$Y_D, times_e = pair$times_e, q = 1,
                                  x_list = list(M = pair$X_M, D = pair$X_D[-1, ])),
               "x_list")
  expect_error(make_prnode_design(m_obs_time = pair$m_obs_time, Y_M = pair$Y_M, d_obs_time = pair$d_obs_time,
                                  Y_D = pair$Y_D, times_e = pair$times_e, q = 1,
                                  x_list = list(M = pair$X_M[seq(1, nrow(pair$X_M), by = 2), ], D = pair$X_D)),
               "x_list")
})

## C5: fixture check at the basis level (Q2 restated through profile_cell()).
test_that("C5: profiled OLS on Psi_tilde from the true trajectory recovers theta to 1e-5 in both cells", {
  .skip_unless_built("profile_cell")
  pair <- .design_pair(p = 6, num_d = 36, grid = .fx_times_e, sigma = 0)
  for(cell in c("M", "D")){
    X <- if(cell == "M") pair$X_M else pair$X_D
    grid <- if(cell == "M") pair$times_e else pair$d_times_e
    obs <- if(cell == "M") pair$m_obs_time else pair$d_obs_time
    Y <- if(cell == "M") pair$Y_M else pair$Y_D
    theta <- if(cell == "M") pair$truth$m_ode_pars else pair$truth$d_ode_pars
    Psi <- make_integrated_basis(X, grid, obs, q = 1)
    prof <- profile_cell(Y, Psi, obs)
    theta_hat <- t(qr.coef(qr(prof$Psi_tilde), prof$y_tilde))
    expect_equal(theta_hat, theta, tolerance = 1e-5, info = cell, ignore_attr = TRUE)
  }
})

## C6: what profiling returns is orthogonal to Z, and y - y_tilde lies in
## span(Z).
test_that("C6: y_tilde and Psi_tilde are orthogonal to Z and the removed part lies in span(Z)", {
  .skip_unless_built("profile_cell")
  set.seed(6)
  obs <- sort(stats::runif(20))
  Psi <- matrix(stats::rnorm(20*6), 20)
  y <- matrix(stats::rnorm(20*3), 20)
  prof <- profile_cell(y, Psi, obs)
  expect_equal(prof$Z, cbind(1, obs - obs[1]), ignore_attr = TRUE)
  expect_equal(prof$n, 20)
  expect_lt(max(abs(crossprod(prof$Z, prof$y_tilde))), 1e-10)
  expect_lt(max(abs(crossprod(prof$Z, prof$Psi_tilde))), 1e-10)
  removed <- y - prof$y_tilde
  expect_equal(qr.fitted(prof$qr_Z, removed), removed, tolerance = 1e-10)
  expect_equal(qr.resid(prof$qr_Z, y), prof$y_tilde, tolerance = 1e-12)
})

## C7 (T1, Lemma 1): of the derivation notes profiling absorbs a time shift, a constant
## offset of a Psi column and any multiple of t - t_1.
test_that("C7: shifting time, adding a constant or c (t - t_1) to a Psi column leaves the profiled quantities unchanged", {
  .skip_unless_built("profile_cell")
  set.seed(7)
  obs <- sort(stats::runif(25))
  Psi <- matrix(stats::rnorm(25*4), 25)
  y <- matrix(stats::rnorm(25*2), 25)
  base <- profile_cell(y, Psi, obs)
  shifted <- profile_cell(y, Psi, obs + 3.7)
  expect_equal(shifted$Psi_tilde, base$Psi_tilde, tolerance = 1e-12)
  expect_equal(shifted$y_tilde, base$y_tilde, tolerance = 1e-12)
  Psi_off <- Psi
  Psi_off[, 2] <- Psi_off[, 2] + 5
  Psi_off[, 3] <- Psi_off[, 3] - 2.5*(obs - obs[1])
  altered <- profile_cell(y + 1.5, Psi_off, obs)
  expect_equal(altered$Psi_tilde, base$Psi_tilde, tolerance = 1e-12)
  expect_equal(altered$y_tilde, base$y_tilde, tolerance = 1e-12)
})

## C8 (resolves Z6): at least three frames and two distinct times.
test_that("C8: fewer than 3 frames, coincident times and a row mismatch each error by message", {
  .skip_unless_built("profile_cell")
  expect_error(profile_cell(matrix(1:2, 2), matrix(1:2, 2), c(0, 1)), "3")
  expect_error(profile_cell(matrix(1:3, 3), matrix(1:3, 3), c(0.5, 0.5, 0.5)), "rank")
  expect_error(profile_cell(matrix(1:4, 4), matrix(1:3, 3), c(0, 0.3, 0.6, 1)))
})

## C9: shapes and regulator-major column order.
test_that("C9: Psi is n x pq with regulators as the outer index", {
  .skip_unless_built("make_integrated_basis")
  grid <- seq(0, 1, by = 1e-3)
  grid_sweep <- expand.grid(q = 1:3, p = c(1, 2, 6), n = c(3, 12, 48))
  for(r in seq_len(nrow(grid_sweep))){
    g <- grid_sweep[r, ]
    set.seed(r)
    x_mat <- sapply(seq_len(g$p), function(k) sin((k + 1)*pi*grid) + k)
    x_mat <- matrix(x_mat, ncol = g$p)
    obs <- grid[round(seq(1, length(grid), length.out = g$n))]
    Psi <- make_integrated_basis(x_mat, grid, obs, q = g$q)
    expect_equal(dim(Psi), c(g$n, g$p*g$q), info = paste(g, collapse = " "))
    if(g$q == 3){
      Psi1 <- make_integrated_basis(x_mat, grid, obs, q = 1)
      expect_identical(Psi[, seq(1, g$p*3, by = 3), drop = FALSE], Psi1, info = paste(g, collapse = " "))
    }
  }
})

## C10 (the closed form of the retired four-block get_cov_mat()) left with
## that helper on 2026-10-22.

## D1 (T4, Lemma 3(a)): of the derivation notes the design identity, exactly.
## D1 is an identity of the design, not of the weights (D2/D8 test those), so
## it passes omega = c(1, 1): with omega = NULL the n_D = 3 case would error in
## estimate_noise_sd() (H8), which D6 tests.
test_that("D1: the scaffold block is exactly the sum of the two deviation blocks, with the stated dimensions", {
  grid_sweep <- expand.grid(q = c(1, 3), p = c(1, 6), num_d = c(3, 12, 48))
  for(r in seq_len(nrow(grid_sweep))){
    g <- grid_sweep[r, ]
    label <- paste("q", g$q, "p", g$p, "n_D", g$num_d)
    pair <- .design_pair(p = g$p, num_d = g$num_d, seed = r)
    design <- .design_from_pair(pair, q = g$q, omega = c(1, 1))
    pq <- g$p*g$q
    N <- 48 + g$num_d
    expect_equal(dim(design$x), c(N, 3*pq), info = label)
    expect_equal(dim(design$y), c(N, g$p), info = label)
    S_block <- design$x[, seq_len(pq), drop = FALSE]
    M_block <- design$x[, pq + seq_len(pq), drop = FALSE]
    D_block <- design$x[, 2*pq + seq_len(pq), drop = FALSE]
    expect_identical(S_block, M_block + D_block, info = label)
    expect_true(all(M_block[48 + seq_len(g$num_d), ] == 0), info = label)
    expect_true(all(D_block[seq_len(48), ] == 0), info = label)
    expect_equal(design$index, rep(seq_len(3*g$p), each = g$q), info = label)
    expect_equal(design$M_groups, g$p + seq_len(g$p), info = label)
    expect_equal(design$D_groups, 2*g$p + seq_len(g$p), info = label)
    expect_equal(c(design$n_M, design$n_D, design$p, design$q), c(48, g$num_d, g$p, g$q), info = label)
  }
})

## D2: the row weights are sqrt(omega_c).
test_that("D2: omega = c(4, 1) doubles the mother rows of x and y and leaves the daughter rows", {
  pair <- .design_pair(p = 3, num_d = 12)
  unit <- .design_from_pair(pair, omega = c(M = 1, D = 1))
  weighted <- .design_from_pair(pair, omega = c(M = 4, D = 1))
  expect_equal(weighted$x[1:48, ], 2*unit$x[1:48, ], tolerance = 1e-12)
  expect_equal(weighted$y[1:48, ], 2*unit$y[1:48, ], tolerance = 1e-12)
  expect_equal(weighted$x[49:60, ], unit$x[49:60, ], tolerance = 1e-12)
  expect_equal(weighted$y[49:60, ], unit$y[49:60, ], tolerance = 1e-12)
  expect_equal(unname(weighted$omega), c(4, 1))
})

## D3: the negation of revision 1's contract: no unpenalized columns, no
## pooled rows.
test_that("D3: index has no NA and the design has n_M + n_D rows", {
  pair <- .design_pair(p = 6, num_d = 12)
  design <- .design_from_pair(pair)
  expect_false(anyNA(design$index))
  expect_equal(nrow(design$x), 60)
  expect_equal(nrow(design$y), 60)
})

## D4: protein relabelling moves each regulator block with its protein, at
## both levels.
test_that("D4: permuting the proteins permutes the columns of y and the regulator blocks of x consistently", {
  pair <- .design_pair(p = 4, num_d = 12)
  perm <- c(3, 1, 4, 2)
  permute <- function(pr){
    pr$Y_M <- pr$Y_M[, perm]
    pr$Y_D <- pr$Y_D[, perm]
    pr$X_M <- pr$X_M[, perm]
    pr$X_D <- pr$X_D[, perm]
    pr
  }
  for(bool_regression in c(TRUE, FALSE)){
    base <- .design_from_pair(pair, bool_regression = bool_regression, omega = c(1, 1))
    permuted <- .design_from_pair(permute(pair), bool_regression = bool_regression, omega = c(1, 1))
    expect_equal(permuted$y, base$y[, perm], tolerance = 1e-12, info = bool_regression)
    col_perm <- c(perm, 4 + perm, 8 + perm)
    expect_equal(permuted$x, base$x[, col_perm], tolerance = 1e-12, info = bool_regression)
    expect_equal(permuted$response_scale, base$response_scale[perm], tolerance = 1e-12, info = bool_regression,
                 ignore_attr = TRUE)
    # each cell's stored mean state moves with its proteins
    expect_equal(permuted$mean_state$M, base$mean_state$M[perm], tolerance = 1e-12, info = bool_regression,
                 ignore_attr = TRUE)
    expect_equal(permuted$mean_state$D, base$mean_state$D[perm], tolerance = 1e-12, info = bool_regression,
                 ignore_attr = TRUE)
  }
})

## D5: the x_list route is the smoothed route when handed the smoothed
## matrices, and differs from it on noisy data.
test_that("D5: x_list = the smoothed matrices reproduces the smoothed route identically; the true trajectories differ", {
  pair <- .design_pair(p = 3, num_d = 12)
  smoothed <- .design_from_pair(pair, bool_regression = FALSE)
  d_grid <- make_daughter_grid(pair$birth_time, pair$times_e)
  x_list <- list(M = smooth_md_trajectory(pair$m_obs_time, pair$Y_M, pair$times_e, spar = 0.5),
                 D = smooth_md_trajectory(pair$d_obs_time, pair$Y_D, d_grid, spar = 0.5))
  via_x_list <- make_prnode_design(m_obs_time = pair$m_obs_time, Y_M = pair$Y_M, d_obs_time = pair$d_obs_time,
                                   Y_D = pair$Y_D, times_e = pair$times_e, q = 1, x_list = x_list)
  expect_identical(via_x_list$x, smoothed$x)
  # everything but the recorded smoothing level, which is NA when x_list supplies the trajectories
  keep <- setdiff(names(smoothed), "spar")
  expect_identical(via_x_list[keep], smoothed[keep])
  expect_true(is.na(via_x_list$spar))
  expect_equal(smoothed$spar, 0.5)
  truth_route <- .design_from_pair(pair, bool_regression = TRUE)
  expect_false(isTRUE(all.equal(truth_route$x, smoothed$x)))
})

## D6: guards, each by message, naming the cell where the plan asks.
test_that("D6: mismatched proteins, off-grid frames, too few frames (3 without omega, 2 at all), a daughter before her mother each error", {
  pair <- .design_pair(p = 3, num_d = 12)
  call_design <- function(...){
    a <- utils::modifyList(list(m_obs_time = pair$m_obs_time, Y_M = pair$Y_M, d_obs_time = pair$d_obs_time,
                                Y_D = pair$Y_D, times_e = pair$times_e, q = 1,
                                x_list = list(M = pair$X_M, D = pair$X_D)), list(...))
    do.call(make_prnode_design, a)
  }
  expect_error(call_design(Y_D = pair$Y_D[, 1:2], x_list = list(M = pair$X_M, D = pair$X_D[, 1:2])), "ncol")
  expect_error(call_design(d_obs_time = pair$d_obs_time + 2e-6))
  short <- .design_pair(p = 3, num_d = 2)
  expect_error(make_prnode_design(m_obs_time = short$m_obs_time, Y_M = short$Y_M, d_obs_time = short$d_obs_time,
                                  Y_D = short$Y_D, times_e = short$times_e, q = 1,
                                  x_list = list(M = short$X_M, D = short$X_D)),
               "daughter")
  # three frames are enough for the design (C8) but not for the noise weights
  # (H8): refused with omega = NULL, naming the cell and the minimum of 4
  three <- .design_pair(p = 3, num_d = 3)
  three_args <- list(m_obs_time = three$m_obs_time, Y_M = three$Y_M, d_obs_time = three$d_obs_time,
                     Y_D = three$Y_D, times_e = three$times_e, q = 1, x_list = list(M = three$X_M, D = three$X_D))
  expect_error(do.call(make_prnode_design, three_args), "daughter")
  expect_error(do.call(make_prnode_design, three_args), "4")
  expect_no_error(do.call(make_prnode_design, c(three_args, list(omega = c(1, 1)))))
  expect_error(call_design(d_obs_time = pair$m_obs_time[1:12], Y_D = pair$Y_M[1:12, ],
                           x_list = list(M = pair$X_M, D = pair$X_M)),
               "birth|first frame|precede")
})

## D7(c)/D8: the response scaling and the pooled weights, in the order the
## documentation must state: scale each protein by its pooled sd, estimate the
## noise sd on the unscaled series and divide by that sd, pool 1 / mean_j(.)^2
## per cell.
test_that("D8: response_scale, omega and y follow the documented arithmetic", {
  .skip_unless_built("estimate_noise_sd")
  pair <- .design_pair(p = 4, num_d = 12, sigma = 0.3)
  scaled <- .design_from_pair(pair, bool_scale_response = TRUE)
  s_vec <- sapply(1:4, function(j) stats::sd(c(pair$Y_M[, j], pair$Y_D[, j])))
  expect_equal(unname(scaled$response_scale), s_vec, tolerance = 1e-12)
  sigma_M <- estimate_noise_sd(pair$m_obs_time, pair$Y_M)
  sigma_D <- estimate_noise_sd(pair$d_obs_time, pair$Y_D)
  omega_M <- 1/mean((sigma_M/s_vec)^2)
  omega_D <- 1/mean((sigma_D/s_vec)^2)
  # one common weight: the frame-weighted pool of the per-cell estimates
  n_M <- nrow(pair$Y_M)
  n_D <- nrow(pair$Y_D)
  w <- (n_M + n_D)/(n_M/omega_M + n_D/omega_D)
  expect_equal(unname(scaled$omega), c(w, w), tolerance = 1e-12)
  # the scaling of y, read with a fixed omega so that only the scale moves
  fixed_scaled <- .design_from_pair(pair, bool_scale_response = TRUE, omega = c(1, 1))
  fixed_raw <- .design_from_pair(pair, bool_scale_response = FALSE, omega = c(1, 1))
  for(j in 1:4){
    expect_equal(fixed_scaled$y[, j], fixed_raw$y[, j]/s_vec[j], tolerance = 1e-12, info = j)
  }
  # the regulators are scaled too (D7(c)): regulator k's column in each of
  # the three blocks is divided by s_k (q = 1), so x is not invariant
  for(block_start in c(0, 4, 8)) for(k in 1:4){
    expect_equal(fixed_scaled$x[, block_start + k], fixed_raw$x[, block_start + k]/s_vec[k], tolerance = 1e-12,
                 info = paste(block_start, k))
  }
  raw <- .design_from_pair(pair, bool_scale_response = FALSE)
  expect_equal(unname(raw$response_scale), rep(1, 4))
  # the mean states are the column means of the observations in their
  # unscaled units, whatever the response scaling (0.13.0: what
  # inheritance_flow() linearizes each cell's field at when q > 1)
  mean_state <- list(M = colMeans(pair$Y_M), D = colMeans(pair$Y_D))
  expect_equal(scaled$mean_state, mean_state, tolerance = 1e-12)
  expect_equal(raw$mean_state, mean_state, tolerance = 1e-12)
  w_raw <- (n_M + n_D)/(n_M*mean(sigma_M^2) + n_D*mean(sigma_D^2))
  expect_equal(unname(raw$omega), c(w_raw, w_raw), tolerance = 1e-12)
})

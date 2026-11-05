# Plan section 13: invariances across the pipeline (I1-I8), through
# make_prnode_design(), fit_prnode_path() and fit_prnode().

.inv_raw <- function(num_d = 12){
  .fixture_dataset("regression", num_d = num_d)
}

.inv_fit_raw <- function(raw, x_list, lambda, mu = c(S = 1, M = 1, D = 1), ...){
  .skip_unless_built("fit_prnode", "public")
  public <- make_prnode_design(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                               Y_D = raw$Y_D, times_e = raw$times_e, q = 1, x_list = x_list, ...)
  fit_prnode(public, lambda = lambda, mu = mu)
}

## I1 (D4 + F2 end to end): protein relabelling.
test_that("I1: permuting the proteins permutes rows and columns of every network and coefficient matrix", {
  ds <- .inv_raw()
  perm <- c(4, 1, 6, 2, 5, 3)
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  base <- .inv_fit_raw(ds$raw, list(M = ds$truth$X_M, D = ds$truth$X_D), lambda)
  raw_p <- ds$raw
  raw_p$Y_M <- raw_p$Y_M[, perm]
  raw_p$Y_D <- raw_p$Y_D[, perm]
  moved <- .inv_fit_raw(raw_p, list(M = ds$truth$X_M[, perm], D = ds$truth$X_D[, perm]), lambda)
  for(what in c("m_network_est", "d_network_est", "s_network_est", "m_ode_est", "d_ode_est", "s_ode_est")){
    expect_equal(unname(moved[[what]]), unname(base[[what]][perm, perm]), tolerance = 1e-10, info = what)
  }
  expect_equal(unname(moved$intercepts), unname(base$intercepts[, perm]), tolerance = 1e-10)
  expect_equal(unname(moved$omega), unname(base$omega), tolerance = 1e-12)
  expect_equal(unname(moved$response_scale), unname(base$response_scale[perm]), tolerance = 1e-12)
})

## I2: time scale. Multiplying every time by c multiplies Psi by c, so the
## standardized problem is identical: supports equal, coefficients and
## intercept slopes divided by c, initial conditions unchanged.
test_that("I2: scaling all times by c divides the ODE coefficients and slopes by c and keeps the supports", {
  ds <- .inv_raw()
  c_time <- 12
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  base <- .inv_fit_raw(ds$raw, list(M = ds$truth$X_M, D = ds$truth$X_D), lambda)
  raw_c <- ds$raw
  raw_c$m_obs_time <- c_time*raw_c$m_obs_time
  raw_c$d_obs_time <- c_time*raw_c$d_obs_time
  raw_c$times_e <- c_time*raw_c$times_e
  raw_c$birth_time <- c_time*raw_c$birth_time
  scaled <- .inv_fit_raw(raw_c, list(M = ds$truth$X_M, D = ds$truth$X_D), lambda)
  expect_identical(.as_adj(scaled$m_network_est), .as_adj(base$m_network_est))
  expect_identical(.as_adj(scaled$d_network_est), .as_adj(base$d_network_est))
  expect_equal(scaled$m_ode_est, base$m_ode_est/c_time, tolerance = 1e-8)
  expect_equal(scaled$d_ode_est, base$d_ode_est/c_time, tolerance = 1e-8)
  expect_equal(scaled$intercepts, base$intercepts/c_time, tolerance = 1e-8)
  expect_equal(scaled$m_init_conditions, base$m_init_conditions, tolerance = 1e-8)
})

## I3 (T14): response scaling. A global unit change with omega re-estimated
## leaves the standardized problem unchanged under either setting of
## bool_scale_response; a per-protein rescaling is absorbed exactly by the
## response scaling (supports unchanged, coefficients a_j / a_k); with a fixed
## user omega and no response scaling, a global scale moves the supports
## unless lambda scales with it.
test_that("I3: y -> a y leaves the supports unchanged with omega re-estimated, and per-protein scaling is absorbed by response scaling", {
  ds <- .inv_raw(num_d = 36)
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  x_list <- list(M = ds$truth$X_M, D = ds$truth$X_D)
  a_global <- 7
  scale_raw <- function(raw, a_vec){
    raw$Y_M <- sweep(raw$Y_M, 2, a_vec, "*")
    raw$Y_D <- sweep(raw$Y_D, 2, a_vec, "*")
    raw
  }
  scale_x <- function(a_vec) list(M = sweep(x_list$M, 2, a_vec, "*"), D = sweep(x_list$D, 2, a_vec, "*"))
  for(bool_scale in c(TRUE, FALSE)){
    base <- .inv_fit_raw(ds$raw, x_list, lambda, bool_scale_response = bool_scale)
    scaled <- .inv_fit_raw(scale_raw(ds$raw, rep(a_global, 6)), scale_x(rep(a_global, 6)), lambda,
                           bool_scale_response = bool_scale)
    expect_identical(.as_adj(scaled$d_network_est), .as_adj(base$d_network_est), info = bool_scale)
    expect_identical(.as_adj(scaled$m_network_est), .as_adj(base$m_network_est), info = bool_scale)
    expect_equal(scaled$m_ode_est, base$m_ode_est, tolerance = 1e-8, info = bool_scale)
    expect_equal(scaled$intercepts, a_global*base$intercepts, tolerance = 1e-8, info = bool_scale)
  }
  a_vec <- c(0.5, 1, 2, 4, 0.25, 3)
  base <- .inv_fit_raw(ds$raw, x_list, lambda, bool_scale_response = TRUE)
  per_protein <- .inv_fit_raw(scale_raw(ds$raw, a_vec), scale_x(a_vec), lambda, bool_scale_response = TRUE)
  expect_identical(.as_adj(per_protein$d_network_est), .as_adj(base$d_network_est))
  expect_identical(.as_adj(per_protein$m_network_est), .as_adj(base$m_network_est))
  expect_equal(per_protein$m_ode_est, base$m_ode_est*outer(a_vec, 1/a_vec), tolerance = 1e-8)
  fixed_base <- .inv_fit_raw(ds$raw, x_list, lambda, omega = c(1, 1), bool_scale_response = FALSE)
  fixed_scaled <- .inv_fit_raw(scale_raw(ds$raw, rep(a_global, 6)), scale_x(rep(a_global, 6)), lambda,
                               omega = c(1, 1), bool_scale_response = FALSE)
  expect_false(identical(.as_adj(fixed_scaled$d_network_est), .as_adj(fixed_base$d_network_est)) &&
                 identical(.as_adj(fixed_scaled$m_network_est), .as_adj(fixed_base$m_network_est)))
  fixed_rescaled <- .inv_fit_raw(scale_raw(ds$raw, rep(a_global, 6)), scale_x(rep(a_global, 6)), a_global*lambda,
                                 omega = c(1, 1), bool_scale_response = FALSE)
  expect_identical(.as_adj(fixed_rescaled$d_network_est), .as_adj(fixed_base$d_network_est))
  expect_equal(fixed_rescaled$m_ode_est, fixed_base$m_ode_est, tolerance = 1e-8)
})

## I4 (T3, Lemma 2(a)): reparametrizing a regulator's basis in both cells.
test_that("I4: Psi_k -> Psi_k A in both cells maps theta_k to A^-1 theta_k with supports and objective unchanged", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_design_q3(12)
  lambda <- .lambda_path(design, num_lambda = 4)[3]
  base <- fit_prnode_path(design, lambda = lambda)
  set.seed(4)
  A <- matrix(stats::rnorm(9), 3) + 3*diag(3)
  transformed <- design
  p <- design$p
  for(block in 0:2) for(k in seq_len(p)){
    cols <- block*p*3 + (k - 1)*3 + 1:3
    transformed$x[, cols] <- design$x[, cols] %*% A
  }
  moved <- fit_prnode_path(transformed, lambda = lambda)
  for(j in seq_len(p)){
    expect_equal(base$objective[j, 1], moved$objective[j, 1], tolerance = 1e-8, info = j)
    for(k in seq_len(p)){
      cols <- (k - 1)*3 + 1:3
      expect_equal(.path_mat(moved, "theta_M")[j, cols], c(solve(A, .path_mat(base, "theta_M")[j, cols])),
                   tolerance = 1e-6, info = paste(j, k))
    }
  }
  for(cell in c("S", "M", "D")) expect_identical(.path_net(moved, cell), .path_net(base, cell), info = cell)
})

## I5 (T2, Lemma 1's last sentence): the cells must share units.
test_that("I5: rescaling one cell's Psi_tilde by 2 changes the fit", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda <- .lambda_path(design, num_lambda = 4)[3]
  base <- fit_prnode_path(design, lambda = lambda)
  altered <- design
  d_rows <- design$n_M + seq_len(design$n_D)
  pq <- design$p*design$q
  altered$x[d_rows, c(seq_len(pq), 2*pq + seq_len(pq))] <- 2*design$x[d_rows, c(seq_len(pq), 2*pq + seq_len(pq))]
  moved <- fit_prnode_path(altered, lambda = lambda)
  same_support <- identical(.path_net(moved, "D"), .path_net(base, "D")) && identical(.path_net(moved, "M"), .path_net(base, "M"))
  same_coef <- isTRUE(all.equal(.path_mat(moved, "theta_D"), .path_mat(base, "theta_D"), tolerance = 1e-6))
  expect_false(same_support && same_coef)
})

## I6 (T12, Goal 3): swap symmetry at equal deviation multipliers, done at
## the solver level because the design guard requires the daughter to start
## no earlier than the mother.
test_that("I6: swapping the cells' rows and deviation blocks swaps theta_M and theta_D and keeps s, at n_M = n_D and n_M != n_D", {
  .skip_unless_built("fit_prnode_path")
  for(num_d in c(48, 12)){
    design <- .fixture_dataset("regression", num_d = num_d)$public
    lambda <- .lambda_path(design, mu = c(S = 1, M = 0.9, D = 0.9), num_lambda = 4)[3]
    base <- fit_prnode_path(design, lambda = lambda, mu = c(S = 1, M = 0.9, D = 0.9))
    pq <- design$p*design$q
    m_rows <- seq_len(design$n_M)
    d_rows <- design$n_M + seq_len(design$n_D)
    swapped <- design
    swapped$x <- design$x[c(d_rows, m_rows), c(seq_len(pq), 2*pq + seq_len(pq), pq + seq_len(pq))]
    swapped$y <- design$y[c(d_rows, m_rows), , drop = FALSE]
    swapped$n_M <- design$n_D
    swapped$n_D <- design$n_M
    swapped$omega <- design$omega[2:1]
    swapped$m_obs_time <- design$d_obs_time
    swapped$d_obs_time <- design$m_obs_time
    swapped$profile <- list(M = design$profile$D, D = design$profile$M)
    moved <- fit_prnode_path(swapped, lambda = lambda, mu = c(S = 1, M = 0.9, D = 0.9))
    expect_equal(.path_mat(moved, "theta_M"), .path_mat(base, "theta_D"), tolerance = 1e-8, info = num_d)
    expect_equal(.path_mat(moved, "theta_D"), .path_mat(base, "theta_M"), tolerance = 1e-8, info = num_d)
    expect_equal(.path_mat(moved, "s"), .path_mat(base, "s"), tolerance = 1e-8, info = num_d)
  }
})

## I7: the loss is a sum over rows.
test_that("I7: permuting the rows within a cell leaves the coefficients identical", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda <- .lambda_path(design, num_lambda = 4)[3]
  base <- fit_prnode_path(design, lambda = lambda)
  set.seed(7)
  perm <- c(sample(design$n_M), design$n_M + sample(design$n_D))
  shuffled <- design
  shuffled$x <- design$x[perm, ]
  shuffled$y <- design$y[perm, , drop = FALSE]
  moved <- fit_prnode_path(shuffled, lambda = lambda)
  for(what in c("s", "d_M", "d_D")) expect_equal(moved[[what]], base[[what]], tolerance = 1e-10, info = what)
})

## I8 (revision 1's I4; R10 through fit_prnode()): identical cells.
test_that("I8: with n_M = n_D, the same data and mu > 1/sqrt(2), theta_M = theta_D and d = 0", {
  .skip_unless_built("fit_prnode", "public")
  ds <- .fixture_dataset("regression", num_d = 12)
  raw <- ds$raw
  fit <- fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$m_obs_time, Y_D = raw$Y_M,
                   times_e = raw$times_e, q = 1, lambda = .lambda_path(ds$public, num_lambda = 4)[3],
                   mu = c(S = 1, M = 0.8, D = 0.8))
  expect_equal(fit$m_ode_est, fit$d_ode_est, tolerance = 1e-10)
  expect_identical(.as_adj(fit$m_network_est), .as_adj(fit$d_network_est))
  expect_true(all(.path_mat(fit$path, "d_M") == 0))
  expect_true(all(.path_mat(fit$path, "d_D") == 0))
  expect_true(any(fit$m_ode_est != 0))
})

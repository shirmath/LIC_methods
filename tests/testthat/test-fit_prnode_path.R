# Plan sections 8 and 9: estimate_noise_sd(), lambda_max_prnode() (H1-H7) and
# the reference solver fit_prnode_path() (R1-R17), in R/fit_prnode_path.R.
# Every test skips until the function it exercises is built.

.mu_settings <- list(c(S = 1, M = 1, D = 1), c(S = 1, M = 0.5, D = 0.5), c(S = 1, M = 2, D = 1), c(S = 1, M = 1, D = 2))

# the hand KKT form of lambda_max per protein: max_g ||Q_g' y_j|| / (sqrt(q) mu_g)
.hand_lambda_max <- function(design, mu = c(S = 1, M = 1, D = 1)){
  groups <- .design_groups(design)
  w_vec <- .group_weights(design, mu)
  sapply(seq_len(ncol(design$y)), function(j){
    max(sapply(seq_along(groups), function(g) sqrt(sum(crossprod(groups[[g]]$Q, design$y[, j])^2))/w_vec[g]))
  })
}

# every coefficient block of protein j at path point l is exactly zero
.all_zero <- function(path, j, l = 1){
  all(.path_beta(path, j, l) == 0)
}

## H1 (T9): lambda_max is the boundary, for several multipliers and both q.
test_that("H1: every group is zero just above lambda_max and some group is nonzero just below it", {
  .skip_unless_built("lambda_max_prnode")
  .skip_unless_built("fit_prnode_path")
  designs <- list(q1 = .fixture_dataset("regression", num_d = 36)$public, q3 = .fixture_design_q3(36))
  for(q_name in names(designs)) for(mu in .mu_settings){
    design <- designs[[q_name]]
    label <- paste(q_name, paste(mu, collapse = ","))
    lambda_vec <- lambda_max_prnode(design, mu = mu)
    expect_length(lambda_vec, design$p)
    expect_equal(lambda_max_prnode(design, mu = mu, bool_pooled = TRUE), max(lambda_vec), info = label)
    above <- fit_prnode_path(design, lambda = (1 + 1e-6)*max(lambda_vec), mu = mu)
    for(j in seq_len(design$p)) expect_true(.all_zero(above, j), info = paste(label, "protein", j))
    below_path <- sort(unique((1 - 1e-3)*lambda_vec), decreasing = TRUE)
    below <- fit_prnode_path(design, lambda = below_path, mu = mu)
    for(j in seq_len(design$p)){
      l <- which(abs(below_path - (1 - 1e-3)*lambda_vec[j]) < 1e-12*lambda_vec[j])[1]
      expect_false(.all_zero(below, j, l), info = paste(label, "protein", j))
    }
  }
})

## H2: lambda_max honours the multipliers and equals the KKT formula at b = 0;
## the multiplier settings are chosen so that the maximizing group is a
## scaffold group (doubling the deviation multipliers changes nothing) or a
## deviation group (doubling halves it).
test_that("H2: lambda_max equals the hand KKT formula and responds to the deviation multipliers as the argmax group dictates", {
  .skip_unless_built("lambda_max_prnode")
  design <- .fixture_dataset("regression", num_d = 12)$public
  for(mu in .mu_settings){
    expect_equal(lambda_max_prnode(design, mu = mu), .hand_lambda_max(design, mu), tolerance = 1e-12,
                 info = paste(mu, collapse = ","), ignore_attr = TRUE)
  }
  scaffold_wins <- lambda_max_prnode(design, mu = c(S = 1, M = 100, D = 100))
  expect_equal(lambda_max_prnode(design, mu = c(S = 1, M = 200, D = 200)), scaffold_wins, tolerance = 1e-12)
  deviation_wins <- lambda_max_prnode(design, mu = c(S = 100, M = 1, D = 1))
  expect_equal(lambda_max_prnode(design, mu = c(S = 100, M = 2, D = 2)), deviation_wins/2, tolerance = 1e-12)
})

## H3/H5: the noise-sd estimator against the generator's sigma; diff2 is
## biased upward on a fast sinusoid, and the size of that bias is what the
## H5 revisit reads.
test_that("H3: both methods are within 20% of sigma on a smooth signal at 48 frames; diff2's bias on 6 cycles is pinned", {
  .skip_unless_built("estimate_noise_sd")
  obs_time <- .fx_obs_time
  signal <- sin(2*pi*obs_time) + 0.3*obs_time
  for(ratio in c(0.1, 0.5, 1.4)){
    sigma <- ratio*stats::sd(signal)
    set.seed(round(100*ratio))
    Y <- signal + matrix(stats::rnorm(48*10, sd = sigma), 48)
    for(method in c("diff2", "spline")){
      est <- estimate_noise_sd(obs_time, Y, method = method)
      expect_length(est, 10)
      expect_equal(mean(est), sigma, tolerance = 0.2, info = paste(method, ratio))
      expect_true(all(abs(est/sigma - 1) < 0.35), info = paste(method, ratio))
    }
  }
  fast <- sin(2*pi*6*obs_time)
  sigma <- 0.5*stats::sd(fast)
  set.seed(6)
  Y_fast <- fast + matrix(stats::rnorm(48*50, sd = sigma), 48)
  est_fast <- estimate_noise_sd(obs_time, Y_fast, method = "diff2")
  closed_form <- sqrt(stats::sd(diff(fast, differences = 2))^2/6 + sigma^2)
  expect_gt(mean(est_fast), sigma)
  expect_equal(mean(est_fast), closed_form, tolerance = 0.1)
  expect_equal(estimate_noise_sd(obs_time, matrix(fast, ncol = 1), method = "diff2"),
               stats::sd(diff(fast, differences = 2))/sqrt(6), tolerance = 1e-12, ignore_attr = TRUE)
  expect_snapshot(cat("diff2 bias (estimate / sigma) on a 6-cycle sinusoid at noise ratio 0.5:", round(mean(est_fast)/sigma, 2)))
})

## H4 (Lemma 2(b)'s of the derivation notes premise): with the pooled omega_c the rows have
## unit noise variance on average over proteins; per protein the ratio is the
## identity omega_c (sigma_cj / s_j)^2, whose mean over proteins is 1 by
## construction and which is 1 for every protein only when the scaled noise
## sds are equal.
##
## The variance claim is about the weight arithmetic, so it is read on a
## design whose omega comes from the realized noise sd (Y minus the noiseless
## Y): on this fixture (noise ratio 0.1, 2 to 5 cycles per movie) the default
## `diff2` estimate is curvature-dominated, 1.0 to 1.6 times sigma for the
## mother and 2 to 11 times for the daughter, which is the H5 revisit item
## (H3 pins the bias). The identity clauses are read on the dataset's own
## omega, whatever its estimate.
test_that("H4: the pooled weights make the profiled residual variance one on average over proteins", {
  .skip_unless_built("estimate_noise_sd")
  ds <- .fixture_dataset("regression", num_d = 36)
  public <- ds$public
  s_vec <- public$response_scale
  omega_realized <- sapply(list(ds$raw$Y_M - ds$truth$Y_M_noiseless, ds$raw$Y_D - ds$truth$Y_D_noiseless),
                           function(err_mat) public$p/sum((apply(err_mat, 2, stats::sd)/s_vec)^2))
  design <- make_prnode_design(m_obs_time = ds$raw$m_obs_time, Y_M = ds$raw$Y_M, d_obs_time = ds$raw$d_obs_time,
                               Y_D = ds$raw$Y_D, times_e = ds$raw$times_e, q = public$q, omega = omega_realized,
                               x_list = list(M = ds$truth$X_M, D = ds$truth$X_D))
  for(cell in c("M", "D")){
    rows <- if(cell == "M") seq_len(48) else 48 + seq_len(36)
    theta <- if(cell == "M") ds$truth$m_ode_pars else ds$truth$d_ode_pars
    Y <- if(cell == "M") ds$raw$Y_M else ds$raw$Y_D
    obs_time <- if(cell == "M") ds$raw$m_obs_time else ds$raw$d_obs_time
    pq <- design$p*design$q
    cell_block <- design$x[rows, (if(cell == "M") pq else 2*pq) + seq_len(pq), drop = FALSE]
    # the true coefficients in the design's units: theta * s_k / s_j
    s_vec <- design$response_scale
    theta_scaled <- sweep(sweep(theta, 2, s_vec, "*"), 1, s_vec, "/")
    resid_var <- sapply(seq_len(design$p), function(j){
      r <- design$y[rows, j] - c(cell_block %*% theta_scaled[j, ])
      sum(r^2)/(length(rows) - 2)
    })
    expect_equal(mean(resid_var), 1, tolerance = 0.15, info = cell)
    # the per-cell pool of the estimates; the design then shares one weight
    # between the cells (D8 pins the pooling), so the identity is read on the
    # per-cell value
    sigma_hat <- estimate_noise_sd(obs_time, Y)
    omega_c <- 1/mean((sigma_hat/s_vec)^2)
    identity_vec <- omega_c*(sigma_hat/s_vec)^2
    expect_equal(mean(identity_vec), 1, tolerance = 1e-12, info = cell)
    expect_false(isTRUE(all.equal(identity_vec, rep(1, design$p))), info = cell)
  }
})

## H7: estimate_noise_sd() is scale-equivariant and shift-invariant for both
## methods, which is what lets H6 estimate on the unscaled series and divide.
test_that("H7: estimate_noise_sd() scales with the series and ignores shifts of the series and of time", {
  .skip_unless_built("estimate_noise_sd")
  set.seed(77)
  obs_time <- .fx_obs_time
  Y <- sin(2*pi*obs_time) %o% c(1, 2, 3) + matrix(stats::rnorm(48*3, sd = 0.2), 48)
  s_vec <- c(0.5, 2, 10)
  for(method in c("diff2", "spline")){
    base <- estimate_noise_sd(obs_time, Y, method = method)
    expect_equal(estimate_noise_sd(obs_time, Y %*% diag(1/s_vec), method = method), base/s_vec,
                 tolerance = 1e-12, info = method, ignore_attr = TRUE)
    expect_equal(estimate_noise_sd(obs_time, sweep(Y, 2, c(5, -3, 100), "+"), method = method), base,
                 tolerance = 1e-10, info = method, ignore_attr = TRUE)
    expect_equal(estimate_noise_sd(obs_time + 12.5, Y, method = method), base, tolerance = 1e-10, info = method,
                 ignore_attr = TRUE)
  }
})

## H8 (review pass 5, audit item 3): below four frames neither method has a
## noise estimate (diff2 needs two second differences; smooth.spline needs
## four points), so the function errors naming the minimum rather than
## returning NA into the weights.
test_that("H8: estimate_noise_sd() errors below four frames, naming the minimum, for both methods", {
  .skip_unless_built("estimate_noise_sd")
  Y3 <- matrix(stats::rnorm(6), 3, 2)
  for(method in c("diff2", "spline")){
    expect_error(estimate_noise_sd(c(0, 0.5, 1), Y3, method = method), "4", info = method)
  }
  Y4 <- matrix(stats::rnorm(8), 4, 2)
  for(method in c("diff2", "spline")){
    est <- estimate_noise_sd(c(0, 0.3, 0.6, 1), Y4, method = method)
    expect_length(est, 2)
    expect_false(anyNA(est), info = method)
  }
})

## R1: at the solution every block is the minimizer of F over that block with
## the others fixed (the block update is the exact block minimizer); an
## optimizer started off the solution cannot do better.
test_that("R1: no block can be improved by a generic optimizer at the solution, q = 1 and q = 3", {
  .skip_unless_built("fit_prnode_path")
  designs <- list(q1 = .fixture_dataset("regression", num_d = 12)$public, q3 = .fixture_design_q3(12))
  for(q_name in names(designs)){
    design <- designs[[q_name]]
    lambda <- .lambda_path(design, num_lambda = 3)[3]
    fit <- fit_prnode_path(design, lambda = lambda, mu = c(S = 1, M = 1, D = 1))
    groups <- .design_groups(design)
    for(j in c(1, 4)){
      beta <- .path_beta(fit, j)
      f_at <- .objective(design, j, beta, lambda)
      for(g in c(1, 2, 7, 13)){
        cols <- groups[[g]]$cols
        f_block <- function(b){
          beta_try <- beta
          beta_try[cols] <- b
          .objective(design, j, beta_try, lambda)
        }
        set.seed(g)
        start <- beta[cols] + stats::rnorm(length(cols), sd = 0.1*(1 + abs(beta[cols])))
        if(length(cols) == 1){
          opt <- stats::optimize(f_block, interval = c(beta[cols] - 3, beta[cols] + 3), tol = 1e-12)
          expect_gte(opt$objective + 1e-8*abs(f_at), f_at, label = paste(q_name, j, g))
          expect_equal(opt$minimum, beta[cols], tolerance = 1e-5, info = paste(q_name, j, g))
        } else {
          opt <- stats::optim(start, f_block, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
          expect_gte(opt$value + 1e-8*abs(f_at), f_at, label = paste(q_name, j, g))
        }
      }
    }
  }
})

## R2 (T5, Lemma 3(a)) was dropped in review pass 5 (2026-10-21): the
## gradient identity grad_s = grad_dM + grad_dD follows from the design
## identity x_S = x_M + x_D for any beta, so a test of it never calls the
## solver; D1 pins the design identity and R3 the solver's KKT conditions.

## R3 (T10): the optimality certificate at every point of a path, and
## back-transform consistency of the reported objective.
test_that("R3: gap <= tol F, kkt_max <= 1e-6 lambda, and objective equals F recomputed on the original design", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 36)$public
  mu <- c(S = 1, M = 0.9, D = 0.9)
  lambda_path <- .lambda_path(design, mu = mu, num_lambda = 20)
  fit <- fit_prnode_path(design, lambda = lambda_path, mu = mu, tol = 1e-10)
  for(l in seq_along(lambda_path)) for(j in seq_len(design$p)){
    label <- paste("lambda", l, "protein", j)
    beta <- .path_beta(fit, j, l)
    F_hand <- .objective(design, j, beta, lambda_path[l], mu)
    expect_equal(fit$objective[j, l], F_hand, tolerance = 1e-10, info = label)
    expect_lte(fit$gap[j, l], 1e-10*F_hand + 1e-14, label = label)
    expect_lte(fit$kkt_max[j, l], 1e-6*lambda_path[l], label = label)
    expect_lte(.kkt_violation(design, j, beta, lambda_path[l], mu), 1e-6*lambda_path[l], label = label)
    expect_true(fit$converged[j, l], info = label)
  }
})

## R4 (T10): independent solvers. p = 1: the exact enumeration; p = 6, q = 1:
## glmnet on the orthonormalized columns (a lasso, since q = 1), and FISTA for
## q = 1 and q = 3. CVXR is not used (review pass 5): FISTA covers q = 3.
test_that("R4: fit_prnode_path() matches the scalar closed form, glmnet and FISTA to 1e-8 relative in the objective", {
  .skip_unless_built("fit_prnode_path")
  sc <- .scalar_design()
  mu <- c(S = 1, M = 0.9, D = 0.9)
  lambda_vec <- .lambda_path(sc$design, mu = mu, num_lambda = 5)[2:5]
  fit <- fit_prnode_path(sc$design, lambda = lambda_vec, mu = mu)
  for(l in seq_along(lambda_vec)){
    exact <- .scalar_exact(sc$design, lambda_vec[l], mu)
    beta <- .path_beta(fit, 1, l)
    expect_equal(.objective(sc$design, 1, beta, lambda_vec[l], mu), exact$value, tolerance = 1e-8, info = l)
    expect_equal(beta, exact$beta, tolerance = 1e-6, info = l)
  }

  design <- .fixture_dataset("regression", num_d = 36)$public
  lambda <- .lambda_path(design, num_lambda = 4)[3]
  fit6 <- fit_prnode_path(design, lambda = lambda, mu = c(S = 1, M = 1, D = 1))
  groups <- .design_groups(design)
  Q <- do.call(cbind, lapply(groups, function(gr) gr$Q))
  N <- nrow(design$x)
  for(j in c(1, 3, 6)){
    beta <- .path_beta(fit6, j)
    F_pkg <- .objective(design, j, beta, lambda)
    if(requireNamespace("glmnet", quietly = TRUE)){
      # glmnet 5 moved thresh and maxit into `control` and warns on the old
      # names; older versions have no `control` argument
      gl_args <- list(x = Q, y = design$y[, j], lambda = c(10*lambda/N, lambda/N), standardize = FALSE,
                      intercept = FALSE)
      if("control" %in% names(formals(glmnet::glmnet))){
        gl_args$control <- list(thresh = 1e-14, maxit = 1e6)
      } else {
        gl_args$thresh <- 1e-14
        gl_args$maxit <- 1e6
      }
      gl <- do.call(glmnet::glmnet, gl_args)
      b_gl <- as.numeric(stats::coef(gl, s = lambda/N, exact = TRUE, x = Q, y = design$y[, j]))[-1]
      beta_gl <- numeric(ncol(design$x))
      for(g in seq_along(groups)) beta_gl[groups[[g]]$cols] <- b_gl[g]/groups[[g]]$R[1, 1]
      expect_equal(.objective(design, j, beta_gl, lambda), F_pkg, tolerance = 1e-8, info = paste("glmnet", j))
      expect_identical(abs(beta_gl) > 1e-8, abs(beta) > 1e-8, info = paste("glmnet", j))
    }
    beta_fista <- .fista_group_lasso(design, j, lambda)
    expect_equal(.objective(design, j, beta_fista, lambda), F_pkg, tolerance = 1e-8, info = paste("fista", j))
  }

  design3 <- .fixture_design_q3(36)
  lambda3 <- .lambda_path(design3, num_lambda = 4)[3]
  fit3 <- fit_prnode_path(design3, lambda = lambda3, mu = c(S = 1, M = 1, D = 1))
  for(j in c(2, 5)){
    beta <- .path_beta(fit3, j)
    beta_fista <- .fista_group_lasso(design3, j, lambda3)
    expect_equal(.objective(design3, j, beta_fista, lambda3), .objective(design3, j, beta, lambda3),
                 tolerance = 1e-8, info = paste("fista q3", j))
  }
})

## R5: the path is a characterization, not a law: R3 holds at every point
## (checked there), and the number of non-nested support steps is pinned as a
## snapshot on one seed.
test_that("R5: supports along the path are nearly nested [snapshot of the non-nested count]", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 36)$public
  lambda_path <- .lambda_path(design, num_lambda = 30)
  fit <- fit_prnode_path(design, lambda = lambda_path, mu = c(S = 1, M = 1, D = 1))
  non_nested <- 0
  for(l in 2:length(lambda_path)){
    for(cell in c("M", "D")){
      prev <- .path_net(fit, cell, l - 1)
      cur <- .path_net(fit, cell, l)
      if(any(prev & !cur)) non_nested <- non_nested + 1
    }
  }
  expect_lte(non_nested, 10)
  expect_snapshot(cat("non-nested support steps over a 30-point path (seed 101, n_D = 36):", non_nested))
})

## R6 (T11): warm and cold starts agree (uniqueness of the fitted values).
test_that("R6: the last point of a 10-point path equals a cold start in fitted values to 1e-8 and in supports exactly", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("ode", num_d = 12)$public
  lambda_path <- .lambda_path(design, num_lambda = 10)
  warm <- fit_prnode_path(design, lambda = lambda_path, mu = c(S = 1, M = 1, D = 1))
  cold <- fit_prnode_path(design, lambda = lambda_path[10], mu = c(S = 1, M = 1, D = 1))
  for(j in seq_len(design$p)){
    fitted_warm <- c(design$x %*% .path_beta(warm, j, 10))
    fitted_cold <- c(design$x %*% .path_beta(cold, j, 1))
    expect_equal(fitted_warm, fitted_cold, tolerance = 1e-8, info = j)
  }
  for(cell in c("S", "M", "D")) expect_identical(.path_net(warm, cell, 10), .path_net(cold, cell, 1), info = cell)
})

## R7: screening is lossless because of the KKT re-admission pass.
test_that("R7: screen = TRUE and FALSE give identical supports and coefficients to 1e-10 on a 30-point path", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda_path <- .lambda_path(design, num_lambda = 30)
  with_screen <- fit_prnode_path(design, lambda = lambda_path, screen = TRUE)
  without <- fit_prnode_path(design, lambda = lambda_path, screen = FALSE)
  for(what in c("s", "d_M", "d_D")) expect_equal(with_screen[[what]], without[[what]], tolerance = 1e-10, info = what)
  for(cell in c("S", "M", "D")) expect_identical(with_screen$networks[[cell]], without$networks[[cell]], info = cell)
})

## R8: a failed step must not look converged (pitfalls row 45).
test_that("R8: max_iter = 1 reports converged = FALSE with a warning naming the protein and lambda, and a finite iterate", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda <- .lambda_path(design, num_lambda = 3)[3]
  expect_warning(fit <- fit_prnode_path(design, lambda = lambda, max_iter = 1), "lambda")
  expect_false(all(fit$converged))
  bad <- which(!fit$converged[, 1])
  expect_true(length(bad) >= 1)
  for(j in bad){
    beta <- .path_beta(fit, j)
    expect_true(all(is.finite(beta)), info = j)
    expect_equal(fit$objective[j, 1], .objective(design, j, beta, lambda), tolerance = 1e-10, info = j)
  }
  msg <- tryCatch(fit_prnode_path(design, lambda = lambda, max_iter = 1), warning = function(w) conditionMessage(w))
  expect_match(msg, "protein")
})

## R9 (T6, Proposition 2): with mu <= 1/sqrt(2) the cells do not borrow and
## each equals its own single-cell standardized group lasso at weight
## lambda sqrt(q) mu; the scaffold is exactly zero.
test_that("R9: at mu in {0.5, 1/sqrt(2)} theta_c equals the single-cell fit of cell c and s = 0", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 36)$public
  pq <- design$p*design$q
  for(mu_val in c(0.5, 1/sqrt(2))){
    mu <- c(S = 1, M = mu_val, D = mu_val)
    lambda <- .lambda_path(design, mu = mu, num_lambda = 4)[3]
    fit <- fit_prnode_path(design, lambda = lambda, mu = mu)
    expect_true(all(.path_mat(fit, "s") == 0), info = mu_val)
    for(cell in c("M", "D")){
      rows <- if(cell == "M") seq_len(design$n_M) else design$n_M + seq_len(design$n_D)
      block <- (if(cell == "M") pq else 2*pq) + seq_len(pq)
      x_cell <- design$x[rows, block, drop = FALSE]
      index_cell <- design$index[block]
      theta_pkg <- .path_mat(fit, paste0("theta_", cell))
      for(j in seq_len(design$p)){
        theta_single <- .single_cell_group_lasso(x_cell, design$y[rows, j], index_cell, weight = lambda*sqrt(design$q)*mu_val)
        expect_equal(theta_pkg[j, ], theta_single, tolerance = 1e-6, info = paste(mu_val, cell, j))
      }
    }
  }
})

## R10 (T8, Proposition 5 of the derivation notes `prop:split`): duplicate cells. With y_D = y_M and Psi_D = Psi_M
## the two networks coincide, and with mu_S = 1 the scaffold carries them
## iff mu > 1/sqrt(2); the switch sits between 0.70 and 0.72.
test_that("R10: duplicated cells give theta_M = theta_D, with d = 0 at mu = 0.72 and s = 0 at mu = 0.70", {
  .skip_unless_built("fit_prnode_path")
  pair <- .fixture_dataset("regression", num_d = 12)
  raw <- pair$raw
  X <- pair$truth$X_M
  design <- make_prnode_design(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$m_obs_time, Y_D = raw$Y_M,
                               times_e = raw$times_e, q = 1, omega = c(1, 1), x_list = list(M = X, D = X))
  lambda <- .lambda_path(design, mu = c(S = 1, M = 0.72, D = 0.72), num_lambda = 4)[3]
  for(mu_val in c(0.70, 0.72)){
    fit <- fit_prnode_path(design, lambda = lambda, mu = c(S = 1, M = mu_val, D = mu_val))
    expect_equal(.path_mat(fit, "theta_M"), .path_mat(fit, "theta_D"), tolerance = 1e-10, info = mu_val)
    expect_true(any(.path_mat(fit, "theta_M") != 0), info = mu_val)
    if(mu_val > 1/sqrt(2)){
      expect_true(all(.path_mat(fit, "d_M") == 0) && all(.path_mat(fit, "d_D") == 0), info = mu_val)
    } else {
      expect_true(all(.path_mat(fit, "s") == 0), info = mu_val)
      expect_equal(.path_mat(fit, "d_M"), .path_mat(fit, "d_D"), tolerance = 1e-10, info = mu_val)
    }
  }
})

## R11: read-out identities, including gamma from the unprofiled Z and Psi
## rebuilt by hand (C6's last clause).
test_that("R11: theta = s + d, gamma = qr.coef(qr_Z, y - Psi theta), and the networks follow the norm rule", {
  .skip_unless_built("fit_prnode_path")
  ds <- .fixture_dataset("regression", num_d = 36)
  design <- ds$public
  lambda <- .lambda_path(design, num_lambda = 4)[3]
  fit <- fit_prnode_path(design, lambda = lambda, edge_tol = 1e-8)
  expect_equal(.path_mat(fit, "theta_M"), .path_mat(fit, "s") + .path_mat(fit, "d_M"), tolerance = 1e-14)
  expect_equal(.path_mat(fit, "theta_D"), .path_mat(fit, "s") + .path_mat(fit, "d_D"), tolerance = 1e-14)
  pq <- design$p*design$q
  for(cell in c("M", "D")){
    # theta is in the design's units, so Psi must carry the regulator scaling
    hc <- .hand_cell(ds, cell, bool_scale_regulators = TRUE)
    theta <- .path_mat(fit, paste0("theta_", cell))
    gamma_hand <- qr.coef(hc$qr_Z, hc$y - hc$Psi %*% t(theta))
    expect_equal(matrix(fit[[paste0("gamma_", cell)]][, , 1], 2), gamma_hand, tolerance = 1e-8, ignore_attr = TRUE,
                 info = cell)
    rows <- if(cell == "M") seq_len(design$n_M) else design$n_M + seq_len(design$n_D)
    block <- (if(cell == "M") pq else 2*pq) + seq_len(pq)
    net_hand <- matrix(FALSE, design$p, design$p)
    for(j in seq_len(design$p)) for(k in seq_len(design$p)){
      cols <- block[(k - 1)*design$q + seq_len(design$q)]
      contribution <- design$x[rows, cols, drop = FALSE] %*% theta[j, (k - 1)*design$q + seq_len(design$q)]
      net_hand[j, k] <- sqrt(sum(contribution^2)) > 1e-8*sqrt(sum(design$y[rows, j]^2))
    }
    expect_identical(.path_net(fit, cell), net_hand, info = cell)
  }
  s_mat <- .path_mat(fit, "s")
  s_net <- matrix(FALSE, design$p, design$p)
  for(j in seq_len(design$p)) for(k in seq_len(design$p)){
    cols <- (k - 1)*design$q + seq_len(design$q)
    s_net[j, k] <- sqrt(sum((design$x[, cols, drop = FALSE] %*% s_mat[j, cols])^2)) > 1e-8*sqrt(sum(design$y[, j]^2))
  }
  expect_identical(.path_net(fit, "S"), s_net)
})

## R12: the reference solver matches groups by label value, which the fork
## does not (Z10).
test_that("R12: relabelling index by 10 * label gives the same fit", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda <- .lambda_path(design, num_lambda = 4)[3]
  relabelled <- design
  relabelled$index <- 10*design$index
  relabelled$M_groups <- 10*design$M_groups
  relabelled$D_groups <- 10*design$D_groups
  base <- fit_prnode_path(design, lambda = lambda)
  moved <- fit_prnode_path(relabelled, lambda = lambda)
  for(what in c("s", "d_M", "d_D")) expect_equal(moved[[what]], base[[what]], tolerance = 1e-12, info = what)
})

## R13 (T7; Proposition 3 of the derivation notes `prop:threshold`, whose threshold
## inequality is eq. (10) `eq:threshold`, and Corollary 1 `cor:leak`, the
## leak probability): the inheritance threshold, draw by draw. The
## corollary's probability is evaluated on this fixture's own weights
## (plan P5(a)): the 0.190 of Table 4 of the derivation notes `tab:borrowing` (row
## mu = 0.926) belongs to its unprofiled hand-built columns, which the
## profiled design cannot reproduce exactly. Numbers are those of the PDF
## compiled 2026-10-20; the labels are stable.
test_that("R13: the daughter takes the mother-only edge exactly when eq. (10), the threshold inequality, is violated; the rate matches Corollary 1 (leak probability)", {
  .skip_unless_built("fit_prnode_path")
  .skip_unless_built("lambda_max_prnode")
  skip_on_cran()
  started <- Sys.time()
  sigma <- 0.5
  mu <- c(S = 1, M = 0.926, D = 0.926)
  base <- .scalar_design(sigma = sigma, seed = 1)
  design <- base$design
  hc <- list(M = qr(cbind(1, base$obs_time - base$obs_time[1])), D = qr(cbind(1, base$d_obs - base$d_obs[1])))
  omega <- unname(design$omega)
  # the null lambda_max: the same design with noise-only responses
  null_design <- function(seed){
    set.seed(seed)
    d <- design
    d$y <- rbind(sqrt(omega[1])*qr.resid(hc$M, matrix(stats::rnorm(48, sd = sigma), 48)),
                 sqrt(omega[2])*qr.resid(hc$D, matrix(stats::rnorm(12, sd = sigma), 12)))
    d
  }
  null_lmax <- sapply(seq_len(400), function(r) lambda_max_prnode(null_design(10000 + r), mu = mu, bool_pooled = TRUE))
  lambda <- unname(stats::quantile(null_lmax, 0.95))
  norms <- sqrt(colSums(design$x^2))
  W <- c(S = mu[["S"]]*norms[1], M = mu[["M"]]*norms[2], D = mu[["D"]]*norms[3])
  a_D <- norms[3]
  upper <- min(W[["D"]], W[["S"]] - W[["M"]])
  expect_gt(upper, 0)
  leak_theory <- stats::pnorm(lambda*upper/a_D, lower.tail = FALSE) + stats::pnorm(lambda*W[["D"]]/a_D, lower.tail = FALSE)

  outcomes <- t(sapply(seq_len(400), function(r){
    sc <- .scalar_design(sigma = sigma, seed = 20000 + r)
    fit <- fit_prnode_path(sc$design, lambda = lambda, mu = mu)
    theta_M <- .path_mat(fit, "theta_M")[1, 1]
    score <- sign(theta_M)*sum(sc$design$x[, 3]*sc$design$y[, 1])
    excluded_by_rule <- (score >= -lambda*W[["D"]]) && (score <= lambda*upper)
    near_boundary <- min(abs(score + lambda*W[["D"]]), abs(score - lambda*upper)) < 1e-6*lambda*W[["D"]]
    c(leak = .path_net(fit, "D")[1, 1], rule_leak = !excluded_by_rule, boundary = near_boundary,
      mother_active = theta_M != 0)
  }))
  expect_true(all(outcomes[, "mother_active"]))
  clear <- !outcomes[, "boundary"]
  expect_identical(outcomes[clear, "leak"], outcomes[clear, "rule_leak"])
  rate <- mean(outcomes[, "leak"])
  se <- sqrt(leak_theory*(1 - leak_theory)/400)
  expect_lt(abs(rate - leak_theory), 2*se + 0.01)
  expect_lt(as.numeric(difftime(Sys.time(), started, units = "secs")), 60)
})

## R14 (Proposition 3 of the derivation notes `prop:threshold`): the discount read-out from
## the design's column norms, clamped at 1 (review pass 5): the
## formula 1 - min(W_D, W_S - W_M) / W_D exceeds 1 whenever W_S < W_M, i.e.
## mu > a_S / a_M, and the package reports min(eta, 1). On this fixture
## (n_D = 12, one sigma_j in both cells) a_S / a_M < 1.12, so mu = 1.2 is the
## clamped regime; the test checks that premise so the clamp is exercised.
test_that("R14: discount$D[k] = min(1, 1 - min(W_D, W_S - W_M) / W_D), in [0, 1], 0 when mu <= 1/sqrt(2), clamped at mu = 1.2", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda <- .lambda_path(design, num_lambda = 3)[2]
  p <- design$p
  ratio_SM <- sapply(seq_len(p), function(k) sqrt(sum(design$x[, k]^2))/sqrt(sum(design$x[, p + k]^2)))
  expect_true(all(ratio_SM < 1.2), info = paste("a_S / a_M:", paste(signif(ratio_SM, 4), collapse = ", ")))
  for(mu_val in c(0.5, 1/sqrt(2), 0.9, 1.2)){
    mu <- c(S = 1, M = mu_val, D = mu_val)
    fit <- fit_prnode_path(design, lambda = lambda, mu = mu)
    for(k in seq_len(p)){
      a_S <- sqrt(sum(design$x[, k]^2))
      a_M <- sqrt(sum(design$x[, p + k]^2))
      a_D <- sqrt(sum(design$x[, 2*p + k]^2))
      W_S <- a_S
      W_M <- mu_val*a_M
      W_D <- mu_val*a_D
      eta_D <- 1 - min(W_D, W_S - W_M)/W_D
      eta_M <- 1 - min(W_M, W_S - W_D)/W_M
      expect_equal(fit$discount$D[k], min(1, eta_D), tolerance = 1e-12, info = paste(mu_val, k))
      expect_equal(fit$discount$M[k], min(1, eta_M), tolerance = 1e-12, info = paste(mu_val, k))
      expect_gte(fit$discount$D[k], 0)
      expect_lte(fit$discount$D[k], 1)
      if(mu_val <= 1/sqrt(2)) expect_equal(fit$discount$D[k], 0, tolerance = 1e-12, info = paste(mu_val, k))
      if(mu_val == 1.2){
        expect_gt(eta_D, 1, label = paste("unclamped eta_D at mu = 1.2, regulator", k))
        expect_equal(fit$discount$D[k], 1, tolerance = 1e-12, info = paste("clamped", k))
      }
    }
  }
})

## R15/R17: guards and the y = 0 column, where the relative criterion is met
## at the first iterate with no NaN.
test_that("R15: increasing lambda warns and is sorted; bad mu, tol, design and NA each error; a zero response gives zeros", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda_path <- .lambda_path(design, num_lambda = 3)
  expect_warning(fit <- fit_prnode_path(design, lambda = rev(lambda_path)), "lambda")
  expect_equal(fit$lambda, lambda_path)
  expect_error(fit_prnode_path(design, lambda = lambda_path[2], mu = c(S = 1, M = -1, D = 1)), "mu")
  expect_error(fit_prnode_path(design, lambda = lambda_path[2], tol = 0), "tol")
  broken <- design
  broken$index <- NULL
  expect_error(fit_prnode_path(broken, lambda = lambda_path[2]), "index")
  na_design <- design
  na_design$y[3, 2] <- NA
  expect_error(fit_prnode_path(na_design, lambda = lambda_path[2]), "NA")
  zero_design <- design
  zero_design$y[, 4] <- 0
  fit0 <- fit_prnode_path(zero_design, lambda = lambda_path[2])
  expect_true(all(.path_beta(fit0, 4) == 0))
  expect_true(fit0$converged[4, 1])
  expect_equal(fit0$gap[4, 1], 0)
  expect_false(anyNA(fit0$objective))
  expect_false(anyNA(fit0$gap))
})

## R18: some random designs stall
## at the max_iter cap, and nothing showed where. `verbose` is a non-negative
## integer: 0 prints nothing, 1 names each protein as its path starts and
## ends, 2 adds one line per lambda (sweeps, certificate, active set and
## whether it converged), 3 adds the sweeps themselves. A stalled fit says
## "not converged" at level 2.
test_that("R18: verbose = 0 is silent; 1 names each protein, 2 each lambda, 3 each sweep; a stall is reported", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  lambda_path <- .lambda_path(design, num_lambda = 3)
  num_resp <- ncol(design$y)
  expect_silent(fit_prnode_path(design, lambda = lambda_path))
  expect_silent(fit_prnode_path(design, lambda = lambda_path, verbose = 0))
  out1 <- utils::capture.output(fit_prnode_path(design, lambda = lambda_path, verbose = 1))
  expect_true(any(grepl(paste0("protein 1 of ", num_resp), out1)))
  expect_true(any(grepl(paste0("protein ", num_resp, " of ", num_resp), out1)))
  expect_false(any(grepl("sweep [0-9]", out1)))
  out2 <- utils::capture.output(fit_prnode_path(design, lambda = lambda_path, verbose = 2))
  expect_equal(sum(grepl("lambda 3 of 3", out2)), num_resp)
  expect_true(any(grepl("converged", out2)))
  expect_false(any(grepl("sweep [0-9]", out2)))
  out3 <- utils::capture.output(fit_prnode_path(design, lambda = lambda_path, verbose = 3))
  expect_true(any(grepl("sweep [0-9]", out3)))
  expect_gt(length(out3), length(out2))
  stall <- utils::capture.output(suppressWarnings(fit_prnode_path(design, lambda = lambda_path[3], max_iter = 1,
                                                             verbose = 2)))
  expect_true(any(grepl("not converged", stall)))
  expect_error(fit_prnode_path(design, lambda = lambda_path, verbose = -1), "verbose")
  expect_error(fit_prnode_path(design, lambda = lambda_path, verbose = "yes"), "verbose")
})

## Z7 (tuning plan, option (b)): the finite-value and Matrix
## entry checks live in .design_blocks(), so a numeric-lambda fit names the
## protein (y) or the regulator and block (x) as tune_prnode() does (E12, E16),
## instead of passing Inf to the solver or calling a Matrix "not numeric".
test_that("Z7: a numeric-lambda fit and lambda_max_prnode() name the non-finite column and refuse a Matrix x", {
  .skip_unless_built("fit_prnode_path")
  design <- .fixture_dataset("regression", num_d = 12)$public
  pq <- design$p*design$q
  bad_y <- design
  bad_y$y[5, 3] <- Inf
  bad_x <- design
  bad_x$x[5, pq + 2] <- Inf
  for(call_name in c("fit_prnode", "fit_prnode_path", "lambda_max_prnode")){
    run <- switch(call_name,
                  fit_prnode = function(d) fit_prnode(d, lambda = 0.1),
                  fit_prnode_path = function(d) fit_prnode_path(d, lambda = 0.1),
                  lambda_max_prnode = function(d) lambda_max_prnode(d))
    message_y <- .error_message(run(bad_y))
    expect_match(message_y, "non-finite", info = call_name)
    expect_match(message_y, design$protein_names[3], info = call_name)
    message_x <- .error_message(run(bad_x))
    expect_match(message_x, "non-finite", info = call_name)
    expect_match(message_x, design$protein_names[2], info = call_name)
    expect_match(message_x, "\\bM\\b", info = call_name)
    sparse <- design
    sparse$x <- Matrix::Matrix(design$x, sparse = TRUE)
    expect_match(.error_message(run(sparse)), "sparse", info = call_name)
  }
})

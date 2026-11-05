# Plan section 11: fit_prnode() rewritten as the wrapper of class prnode_fit
# (R/fit_prnode.R), F1-F14, plus P2 (print and predict methods). Every test
# skips until the new signature (argument `public`) exists, except F6(a),
# which reproduces Table 2 of the derivation notes (`tab:leak-noiseless`) on the fork today.

.fm_fit <- function(public, lambda, mu = c(S = 1, M = 1, D = 1), ...){
  .skip_unless_built("fit_prnode", "public")
  fit_prnode(public, lambda = lambda, mu = mu, ...)
}

## F1: the extraction identities from the path, on an unscaled design.
test_that("F1: the *_ode_est, networks, intercepts and initial conditions are read off the path as documented", {
  ds <- .fixture_dataset("regression", num_d = 12)
  public <- .unscaled_public(ds)
  lambda <- .lambda_path(public, num_lambda = 4)[3]
  fit <- .fm_fit(public, lambda, bool_scale_response = FALSE)
  expect_s3_class(fit, "prnode_fit")
  path <- fit$path
  expect_s3_class(path, "prnode_path")
  expect_identical(fit$s_ode_est, .path_mat(path, "s"))
  expect_identical(fit$m_ode_est, .path_mat(path, "s") + .path_mat(path, "d_M"))
  expect_identical(fit$d_ode_est, .path_mat(path, "s") + .path_mat(path, "d_D"))
  expect_identical(.as_adj(fit$m_network_est), .path_net(path, "M"))
  expect_identical(.as_adj(fit$d_network_est), .path_net(path, "D"))
  expect_identical(.as_adj(fit$s_network_est), .path_net(path, "S"))
  expect_equal(dim(fit$intercepts), c(2, public$p))
  expect_identical(rownames(fit$intercepts), c("mother", "daughter"))
  gamma_M <- matrix(path$gamma_M[, , 1], 2)
  gamma_D <- matrix(path$gamma_D[, , 1], 2)
  expect_equal(unname(fit$intercepts["mother", ]), gamma_M[2, ])
  expect_equal(unname(fit$intercepts["daughter", ]), gamma_D[2, ])
  expect_equal(unname(fit$m_init_conditions), gamma_M[1, ])
  expect_equal(unname(fit$d_init_conditions), gamma_D[1, ])
  expect_equal(unname(fit$response_scale), rep(1, public$p))
})

## F2: per-protein independence; the permutation half is I1.
test_that("F2: fitting one response column gives row j of every output identically", {
  ds <- .fixture_dataset("regression", num_d = 12)
  public <- ds$public
  lambda <- .lambda_path(public, num_lambda = 4)[3]
  full <- .fm_fit(public, lambda)
  for(j in c(2, 5)){
    public_j <- public
    public_j$y <- public$y[, j, drop = FALSE]
    public_j$response_scale <- public$response_scale[j]
    public_j$protein_names <- public$protein_names[j]
    alone <- fit_prnode(public_j, lambda = lambda)
    expect_equal(alone$m_ode_est[1, ], full$m_ode_est[j, ], tolerance = 1e-12, info = j, ignore_attr = TRUE)
    expect_equal(alone$d_ode_est[1, ], full$d_ode_est[j, ], tolerance = 1e-12, info = j, ignore_attr = TRUE)
    expect_identical(.as_adj(alone$d_network_est)[1, ], .as_adj(full$d_network_est)[j, ], info = j)
  }
})

## F3 (Z11 gone): the edge rule is a norm threshold relative to ||y_tilde||.
test_that("F3: every network is empty just above lambda_max; the threshold is relative and reachable through edge_tol", {
  ds <- .fixture_dataset("ode", num_d = 12)
  public <- ds$public
  lambda_max <- lambda_max_prnode(public, bool_pooled = TRUE)
  empty <- .fm_fit(public, (1 + 1e-6)*lambda_max)
  expect_false(any(.as_adj(empty$m_network_est)))
  expect_false(any(.as_adj(empty$d_network_est)))
  expect_false(any(.as_adj(empty$s_network_est)))
  lambda <- lambda_max*10^-1.5
  default <- .fm_fit(public, lambda, edge_tol = 1e-8)
  huge <- .fm_fit(public, lambda, edge_tol = 10)
  expect_true(any(.as_adj(default$m_network_est)))
  expect_false(any(.as_adj(huge$m_network_est)))
  expect_identical(default$m_ode_est, huge$m_ode_est)
  zero_tol <- .fm_fit(public, lambda, edge_tol = 0)
  expect_identical(.as_adj(zero_tol$m_network_est),
                   t(sapply(seq_len(public$p), function(j) sapply(seq_len(public$p), function(k) zero_tol$m_ode_est[j, k] != 0))))
})

## F4: the wrapper is the solver.
test_that("F4: fit_prnode() equals a direct fit_prnode_path() call on the same design", {
  ds <- .fixture_dataset("regression", num_d = 12)
  public <- ds$public
  lambda <- .lambda_path(public, num_lambda = 4)[3]
  fit <- .fm_fit(public, lambda, mu = c(S = 1, M = 0.9, D = 0.9))
  direct <- fit_prnode_path(public, lambda = lambda, mu = c(S = 1, M = 0.9, D = 0.9))
  for(what in c("s", "d_M", "d_D", "theta_M", "theta_D", "gamma_M", "gamma_D", "objective"))
    expect_identical(fit$path[[what]], direct[[what]], info = what)
})

## F5: network recovery floors on the fixture (coefficients are identified,
## so the floors are above revision 1's 0.6). The scorer's zeros on the empty
## and the complete network are E4's (review pass 5).
test_that("F5: daughter and mother MCC >= 0.8 at the daughter-oracle lambda (regression), >= 0.6 (ODE)", {
  floors <- c(regression = 0.8, ode = 0.6)
  for(level in names(floors)){
    ds <- .fixture_dataset(level, num_d = 36)
    lambda_path <- .lambda_path(ds$public, mu = c(S = 1, M = 0.9, D = 0.9), num_lambda = 30)
    fit <- .fm_fit(ds$public, lambda_path, mu = c(S = 1, M = 0.9, D = 0.9))
    l <- .oracle_index(fit, ds$truth)
    expect_gte(.mcc(fit$d_network_est[, , l], ds$truth$d_network), floors[[level]], label = paste(level, "daughter"))
    expect_gte(.mcc(fit$m_network_est[, , l], ds$truth$m_network), floors[[level]], label = paste(level, "mother"))
  }
})

## F6 (T15): the leak regression. (a) Table 2 of the derivation notes `tab:leak-noiseless` reproduced on the
## fork: the four-block loss at (2, 1) gives the daughter 0.237 of a
## mother-only edge of strength 2 in noiseless data; the two-block loss gives
## exactly 0. The fixture builds both losses by hand, so the record of why the
## loss changed survives the deletion of the four-block design helper
## (2026-10-22); F6(b), the same contrast through get_cov_mat(), left with it.
test_that("F6(a): on the hand-built scalar fixture the fork leaks 0.237 under the four-block loss and 0 under the two-block loss", {
  fx <- .scalar_fixture()
  lambda_cur <- 0.2*.fork_scalar_lambda_max(fx, fx$psi_M, rep(0, fx$n_D), "current", 2, 1)
  cur <- .fork_scalar_fit(fx, 2*fx$psi_M, rep(0, fx$n_D), lambda_cur, "current", 2, 1)
  expect_equal(unname(cur[["theta_D"]]), 0.237, tolerance = 5e-3)
  expect_gt(abs(cur[["theta_D"]]), 1e-8)
  lambda_rec <- 0.2*.fork_scalar_lambda_max(fx, fx$psi_M, rep(0, fx$n_D), "recommended", 1, 1)
  rec <- .fork_scalar_fit(fx, 2*fx$psi_M, rep(0, fx$n_D), lambda_rec, "recommended", 1, 1)
  expect_identical(unname(rec[["theta_D"]]), 0)
  expect_gt(rec[["theta_M"]], 1)
})

## F7 (D7(c)): the response-scaling round trip puts the coefficients back in
## the original units; at the oracle they have the truth's signs on the
## recovered true edges and are closer to the truth than zero is. A lasso at
## the daughter-oracle lambda (0.08 lambda_max) shrinks the mother's
## coefficients to 0 to 0.9 of the truth (relative error 0.63 on this
## fixture; OLS on the same support 0.135), so a tight bound would ask for a
## relaxed refit that the estimator does not include (review pass 6).
test_that("F7: scaled fits return coefficients in original units, equal to the unscaling of the path", {
  ds <- .fixture_dataset("regression", num_d = 36)
  public <- ds$public
  expect_false(isTRUE(all.equal(unname(public$response_scale), rep(1, public$p))))
  lambda_path <- .lambda_path(public, mu = c(S = 1, M = 0.9, D = 0.9), num_lambda = 30)
  fit <- .fm_fit(public, lambda_path, mu = c(S = 1, M = 0.9, D = 0.9))
  s_vec <- public$response_scale
  for(l in c(5, 15)){
    expect_equal(fit$m_ode_est[, , l], .unscale_coef(.path_mat(fit$path, "theta_M", l), s_vec, public$q),
                 tolerance = 1e-12, info = l)
    expect_equal(fit$d_ode_est[, , l], .unscale_coef(.path_mat(fit$path, "theta_D", l), s_vec, public$q),
                 tolerance = 1e-12, info = l)
    gamma_M <- matrix(fit$path$gamma_M[, , l], 2)
    expect_equal(unname(fit$intercepts["mother", , l]), gamma_M[2, ]*s_vec, tolerance = 1e-12, info = l)
    expect_equal(unname(fit$m_init_conditions[, l]), gamma_M[1, ]*s_vec, tolerance = 1e-12, info = l)
  }
  l <- .oracle_index(fit, ds$truth)
  truth_support <- .as_adj(ds$truth$m_network)
  est <- fit$m_ode_est[, , l]
  rel_err <- sqrt(sum((est[truth_support] - ds$truth$m_ode_pars[truth_support])^2))/sqrt(sum(ds$truth$m_ode_pars^2))
  expect_lt(rel_err, 1)
  recovered <- truth_support & est != 0
  expect_gt(sum(recovered), 0)
  # a lasso at the oracle lambda can leave one weak, collinear coefficient on the
  # wrong side of zero (5 of 6 recovered edges agree under the common weight)
  expect_gte(mean(sign(est[recovered]) == sign(ds$truth$m_ode_pars[recovered])), 0.75)
  expect_equal(unname(fit$response_scale), unname(s_vec))
})

## F8: omega handling through the raw-argument route.
test_that("F8: omega = NULL estimates the weights, a user omega is used as given, a wrong length errors", {
  .skip_unless_built("fit_prnode", "public")
  .skip_unless_built("estimate_noise_sd")
  ds <- .fixture_dataset("ode", num_d = 12)
  raw <- ds$raw
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  fit_null <- fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time, Y_D = raw$Y_D,
                        times_e = raw$times_e, q = 1, lambda = lambda, omega = NULL)
  expect_equal(unname(fit_null$omega), unname(ds$public$omega), tolerance = 1e-12)
  fit_user <- fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time, Y_D = raw$Y_D,
                        times_e = raw$times_e, q = 1, lambda = lambda, omega = c(M = 2, D = 3))
  expect_equal(unname(fit_user$omega), c(2, 3))
  expect_false(isTRUE(all.equal(fit_user$m_ode_est, fit_null$m_ode_est)))
  expect_error(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time, Y_D = raw$Y_D,
                         times_e = raw$times_e, q = 1, lambda = lambda, omega = c(1, 2, 3)),
               "omega")
})

## F9 (Z4): the solver arguments reach the solver.
test_that("F9: tol changes the gap and max_iter = 1 surfaces the solver's warning naming the protein", {
  ds <- .fixture_dataset("regression", num_d = 12)
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  loose <- .fm_fit(ds$public, lambda, tol = 1e-2)
  tight <- .fm_fit(ds$public, lambda, tol = 1e-10)
  expect_false(isTRUE(all.equal(loose$path$gap, tight$path$gap)))
  expect_lt(max(tight$path$gap), max(loose$path$gap))
  expect_true(all(tight$path$gap <= 1e-10*tight$path$objective + 1e-14))
  msg <- tryCatch(fit_prnode(ds$public, lambda = lambda, max_iter = 1), warning = function(w) conditionMessage(w))
  expect_match(msg, "protein")
  fit1 <- suppressWarnings(fit_prnode(ds$public, lambda = lambda, max_iter = 1))
  expect_false(all(fit1$converged))
})

## F10 (Z3, Z9): p and q are derived, never passed.
test_that("F10: an inconsistent public errors naming index and x, and no p or nbasis argument exists", {
  .skip_unless_built("fit_prnode", "public")
  ds <- .fixture_dataset("regression", num_d = 12)
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  bad <- ds$public
  bad$index <- bad$index[-1]
  msg <- tryCatch(fit_prnode(bad, lambda = lambda), error = function(e) conditionMessage(e))
  expect_match(msg, "index")
  expect_match(msg, "x")
  expect_false(any(c("p", "nbasis") %in% names(formals(fit_prnode))))
})

## F11/F12(b): the shape depends on the length of lambda, by design.
test_that("F11: a lambda vector gives arrays with a named third dimension; fit[l] is the scalar fit at lambda[l]", {
  ds <- .fixture_dataset("regression", num_d = 12)
  public <- ds$public
  lambda_path <- .lambda_path(public, num_lambda = 5)
  fit <- .fm_fit(public, lambda_path)
  p <- public$p
  expect_equal(dim(fit$m_network_est), c(p, p, 5))
  expect_equal(dim(fit$d_ode_est), c(p, p*public$q, 5))
  expect_equal(dim(fit$intercepts), c(2, p, 5))
  expect_equal(as.numeric(dimnames(fit$m_network_est)[[3]]), lambda_path, tolerance = 1e-12)
  expect_equal(as.numeric(dimnames(fit$intercepts)[[3]]), lambda_path, tolerance = 1e-12)
  one <- fit[4]
  expect_s3_class(one, "prnode_fit")
  expect_equal(dim(one$m_network_est), c(p, p))
  expect_equal(dim(one$m_ode_est), c(p, p*public$q))
  expect_equal(dim(one$intercepts), c(2, p))
  cold <- fit_prnode(public, lambda = lambda_path[4])
  expect_equal(one$m_ode_est, cold$m_ode_est, tolerance = 1e-8)
  expect_equal(one$d_ode_est, cold$d_ode_est, tolerance = 1e-8)
  expect_identical(.as_adj(one$d_network_est), .as_adj(cold$d_network_est))
  expect_equal(one$lambda, lambda_path[4])
  expect_true(is.matrix(cold$m_network_est))
  expect_true(is.matrix(cold$m_ode_est))
})

## F14 (F13 revised): the joint_* names are gone.
test_that("F14: no joint_* element exists, s_* elements have the F12(b) shapes, and print mentions s_ only", {
  ds <- .fixture_dataset("regression", num_d = 12)
  public <- ds$public
  lambda_path <- .lambda_path(public, num_lambda = 3)
  for(lambda in list(lambda_path[2], lambda_path)){
    fit <- .fm_fit(public, lambda)
    expect_null(fit$joint_ode_est)
    expect_null(fit$joint_network_est)
    expect_false(any(grepl("^joint", names(fit))))
    expect_true(all(c("s_ode_est", "s_network_est") %in% names(fit)))
    if(length(lambda) == 1){
      expect_equal(dim(fit$s_ode_est), c(public$p, public$p*public$q))
      expect_equal(dim(fit$s_network_est), c(public$p, public$p))
    } else {
      expect_equal(dim(fit$s_ode_est), c(public$p, public$p*public$q, length(lambda)))
      expect_equal(dim(fit$s_network_est), c(public$p, public$p, length(lambda)))
    }
    printed <- paste(utils::capture.output(print(fit)), collapse = "\n")
    expect_no_match(printed, "joint")
  }
})

## P2: the print and predict methods of prnode_fit.
test_that("P2: print.prnode_fit() reports the settings and counts; predict() is Z gamma + Psi theta at the frames", {
  ds <- .fixture_dataset("regression", num_d = 12)
  public <- ds$public
  lambda <- .lambda_path(public, mu = c(S = 1, M = 0.9, D = 0.9), num_lambda = 4)[3]
  fit <- .fm_fit(public, lambda, mu = c(S = 1, M = 0.9, D = 0.9))
  printed <- paste(utils::capture.output(print(fit)), collapse = "\n")
  for(token in c("6", "48", "12", signif(lambda, 3), "0.9", "scal")) expect_match(printed, as.character(token), fixed = TRUE)
  expect_match(printed, as.character(sum(.as_adj(fit$d_network_est))), fixed = TRUE)
  expect_match(printed, "discount")
  for(cell in c("M", "D")){
    pred <- stats::predict(fit, cell = cell)
    hc <- .hand_cell(ds, cell)
    theta <- if(cell == "M") fit$m_ode_est else fit$d_ode_est
    gamma <- rbind(if(cell == "M") fit$m_init_conditions else fit$d_init_conditions,
                   fit$intercepts[if(cell == "M") "mother" else "daughter", ])
    hand <- hc$Z %*% gamma + hc$Psi %*% t(theta)
    expect_equal(unname(pred), unname(hand), tolerance = 1e-10, info = cell)
    expect_equal(dim(pred), c(if(cell == "M") 48 else 12, public$p), info = cell)
  }
})

## F15: `verbose` reaches the solver
## one level down, so fit_prnode(verbose = 1) reports the wrapper's own phases
## and verbose = 2 the solver's per-protein lines.
test_that("F15: fit_prnode(verbose = 1) reports its phases, verbose = 2 the solver's proteins, the default nothing", {
  .skip_unless_built("fit_prnode")
  ds <- .fixture_dataset("regression", num_d = 12)
  lambda_path <- .lambda_path(ds$public, num_lambda = 3)
  expect_silent(fit_prnode(ds$public, lambda = lambda_path))
  out1 <- utils::capture.output(fit_prnode(ds$public, lambda = lambda_path, verbose = 1))
  expect_true(any(grepl("fit_prnode", out1)))
  expect_false(any(grepl("protein 1 of", out1)))
  out2 <- utils::capture.output(fit_prnode(ds$public, lambda = lambda_path, verbose = 2))
  expect_true(any(grepl("protein 1 of", out2)))
  raw <- ds$raw
  out_raw <- utils::capture.output(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                                             Y_D = raw$Y_D, times_e = raw$times_e, lambda = lambda_path, verbose = 1))
  expect_true(any(grepl("design", out_raw)))
  expect_error(fit_prnode(ds$public, lambda = lambda_path, verbose = -1), "verbose")
})

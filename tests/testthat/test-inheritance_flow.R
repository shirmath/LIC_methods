# The inheritance measure (R/inheritance.R): the one-frame flow-map
# similarity (the definition), its estimator on a tuned fit (since 0.9.0 the
# standardized plug-in at the known-variance BIC pick, no window; the pipeline comparison),
# the permutation p-value and the wrong-mother test, and the plot's numbers.
#
# Oracles: a hand-coded field and central differences for the Jacobian at
# q = 3; the first-order (Frobenius) form of the similarity at a short
# horizon; the fit's own unscaled coefficients at the pick
# (`fit$m_ode_est`, `fit$d_ode_est`, the second implementation of the
# unscaling in the package) conjugated to standardized units by hand; the
# p-value counted by hand.

# --- fixtures -----------------------------------------------------------------

# the p = 4 examples' dataset tuned at the defaults (cached in .fx by the
# tuning helper); the fit at that tuning carries `tuning`
.flow_fit_small <- function(){
  key <- "flow_fit_small"
  if(!exists(key, envir = .fx)){
    tn <- .fixture_tuning_default_small()
    assign(key, .with_warnings(fit_prnode(.fixture_dataset_small()$public, lambda = tn))$value, envir = .fx)
  }
  get(key, envir = .fx)
}

# a q = 3 dataset: p = 4, 36 mother frames, an 18-frame daughter, cubic
# truth; tuned on a 20-point, 2-decade
# path so that the fit stays quick
.flow_dataset_q3 <- function(){
  key <- "flow_dataset_q3"
  if(!exists(key, envir = .fx)){
    assign(key,
           build_md_dataset(seed = 3, num_d = 18, noise_level_ratio = 0.1,
                            truth_args = list(p = 4, q = 3, num_shared = 4, num_mother_only = 1,
                                              num_daughter_only = 1),
                            obs_time = seq(0, 1, length.out = 36), steps_per_frame = 200, q = 3,
                            level = "ode"),
           envir = .fx)
  }
  get(key, envir = .fx)
}

.flow_fit_q3 <- function(){
  key <- "flow_fit_q3"
  if(!exists(key, envir = .fx)){
    ds <- .flow_dataset_q3()
    tn <- .with_warnings(tune_prnode(ds$public, num_lambda = 20, num_decades = 2))$value
    assign(key, .with_warnings(fit_prnode(ds$public, lambda = tn))$value, envir = .fx)
  }
  get(key, envir = .fx)
}

# the field f_j(x) = sum_k sum_m pars[j, (k - 1) q + m] x_k^m, as a function
.poly_field <- function(pars, q){
  p <- nrow(pars)
  function(x){
    vapply(seq_len(p), function(j){
      sum(vapply(seq_len(p), function(k) sum(pars[j, (k - 1)*q + seq_len(q)]*x[k]^seq_len(q)), numeric(1)))
    }, numeric(1))
  }
}

# central-difference Jacobian of a vector field
.numeric_jacobian <- function(f, x, eps = 1e-6){
  p <- length(x)
  jac <- matrix(0, p, p)
  for(k in seq_len(p)){
    e_k <- replace(numeric(p), k, eps)
    jac[, k] <- (f(x + e_k) - f(x - e_k))/(2*eps)
  }
  jac
}

# the first-order form of the similarity: 1 - ||A - B|| / (||A|| + ||B||)
.frobenius_sim <- function(A, B){
  1 - norm(A - B, "F")/(norm(A, "F") + norm(B, "F"))
}

# --- the definition -----------------------------------------------------------

## FM1: identities of the definition at q = 1 on random systems: in [0, 1],
## symmetric in the two cells, 1 on a clone, 1 when both are empty (the
## convention), and 0 when the two one-frame displacements are opposite
## (A_D = -A_M at a short horizon is opposite only to first order, so the
## bound is loose there).
test_that("FM1: flow_map_similarity() lies in [0, 1], is symmetric, is 1 on a clone and on two empty systems, and near 0 on opposite systems", {
  grid <- expand.grid(seed = 1:5, p = c(3, 6))
  for(r in seq_len(nrow(grid))){
    set.seed(grid$seed[r])
    p <- grid$p[r]
    A <- matrix(stats::rnorm(p*p), p, p)
    B <- matrix(stats::rnorm(p*p), p, p)
    label <- paste("seed", grid$seed[r], "p", p)
    tau <- flow_map_similarity(A, B, h = 0.05)
    expect_true(is.finite(tau) && tau >= 0 && tau <= 1, info = label)
    expect_equal(tau, flow_map_similarity(B, A, h = 0.05), info = label)
    expect_equal(flow_map_similarity(A, A, h = 0.05), 1, info = label)
    expect_equal(flow_map_similarity(A, A*1.0000001, h = 0.05), 1, tolerance = 1e-5, info = label)
    expect_lt(flow_map_similarity(A, -A, h = 1e-3), 0.01, label = label)
  }
  expect_identical(flow_map_similarity(matrix(0, 3, 3), matrix(0, 3, 3), h = 0.1), 1)
  expect_equal(flow_map_similarity(matrix(0, 3, 3), diag(3), h = 0.1), 0)
})

## FM2 [oracle]: to first order in h the similarity is the Frobenius form
## 1 - ||A - B|| / (||A|| + ||B||), so at a short horizon the two agree to
## O(h), and a daughter who keeps every edge at a fraction rho reads
## 1 - |1 - rho| / (1 + rho) (the appendix's second property).
test_that("FM2: at a short horizon flow_map_similarity() equals the Frobenius similarity to O(h), and a uniformly reweighted daughter reads 1 - |1 - rho| / (1 + rho)", {
  set.seed(2)
  A <- matrix(stats::rnorm(36), 6, 6)
  B <- matrix(stats::rnorm(36), 6, 6)
  h_vec <- c(1e-1, 1e-2, 1e-3, 1e-4)
  gap_vec <- vapply(h_vec, function(h) abs(flow_map_similarity(A, B, h = h) - .frobenius_sim(A, B)), numeric(1))
  expect_lt(gap_vec[4], 1e-4)
  # the gap shrinks with h (it is O(h))
  expect_true(all(diff(gap_vec) < 0))
  for(rho in c(0.25, 0.5, 2)){
    expect_equal(flow_map_similarity(A, rho*A, h = 1e-4), 1 - abs(1 - rho)/(1 + rho), tolerance = 1e-3,
                 info = paste("rho", rho))
  }
})

## FM3 [oracle]: at q = 3 each system matrix is the Jacobian of that cell's
## polynomial field at the cell's own state (`x_ref$M` for the mother,
## `x_ref$D` for the daughter; package 0.13.0), checked against central differences of the
## hand-coded field; at 0 for both cells it is the first-power coefficients,
## so the q = 3 similarity at 0 equals the q = 1 similarity of the
## first-power matrices. The old single state (a bare vector, both cells at
## the mother's mean) is refused, as is a list missing a cell: the two
## readings differ, so a silent fallback would change the measure.
test_that("FM3: at q = 3 flow_map_similarity() linearizes each cell at its own x_ref: the Jacobians match central differences and at 0 they are the first powers", {
  ds <- .flow_dataset_q3()
  truth <- ds$truth
  # a dataset's truth carries no `p` or `q`; both are read off the matrix
  p <- nrow(truth$m_ode_pars)
  q <- ncol(truth$m_ode_pars)/p
  expect_identical(q, 3)
  x_ref <- list(M = colMeans(ds$raw$Y_M), D = colMeans(ds$raw$Y_D))
  jac_M <- LICmethods:::.jacobian_at(truth$m_ode_pars, q, x_ref$M)
  expect_equal(jac_M, .numeric_jacobian(.poly_field(truth$m_ode_pars, q), x_ref$M), tolerance = 1e-5)
  jac_D <- LICmethods:::.jacobian_at(truth$d_ode_pars, q, x_ref$D)
  expect_equal(jac_D, .numeric_jacobian(.poly_field(truth$d_ode_pars, q), x_ref$D), tolerance = 1e-5)
  h <- 1/35
  tau_q3 <- flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = h, q = q, x_ref = x_ref)
  expect_equal(tau_q3, LICmethods:::.flow_map_similarity(jac_M, jac_D, h))
  # both cells at the mother's mean (the reading until 0.12.0) is another number
  x_mother <- list(M = x_ref$M, D = x_ref$M)
  tau_mother <- flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = h, q = q, x_ref = x_mother)
  expect_equal(tau_mother, LICmethods:::.flow_map_similarity(jac_M, LICmethods:::.jacobian_at(truth$d_ode_pars, q, x_ref$M), h))
  expect_false(isTRUE(all.equal(tau_q3, tau_mother)))
  first_cols <- (seq_len(p) - 1)*q + 1
  expect_equal(flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = h, q = q,
                                   x_ref = list(M = rep(0, p), D = rep(0, p))),
               flow_map_similarity(truth$m_ode_pars[, first_cols], truth$d_ode_pars[, first_cols], h = h))
  # the clone and the range hold at q = 3 too: a clone shares the mother's
  # frames, so both states are hers and the two Jacobians coincide
  expect_equal(flow_map_similarity(truth$m_ode_pars, truth$m_ode_pars, h = h, q = q, x_ref = x_mother), 1)
  expect_true(tau_q3 >= 0 && tau_q3 <= 1)
  # x_ref is required at q > 1 as a list with both cells, and the dimensions must agree
  expect_error(flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = h, q = q), "x_ref")
  expect_error(flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = h, q = q, x_ref = x_ref$M), "x_ref")
  expect_error(flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = h, q = q, x_ref = x_ref["M"]), "x_ref")
  expect_error(flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = h, q = q,
                                   x_ref = list(M = x_ref$M, D = 1:3)), "x_ref")
  expect_error(flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars[, first_cols], h = h, q = q), "dimensions")
  expect_error(flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = 0, q = q, x_ref = x_ref))
})

## FM4: overflow returns NA rather than an error or a wrong number (a q = 3
## path's dense end reaches raw coefficients of 1e5 and more; the flow-map study).
test_that("FM4: flow_map_similarity() is NA when the exponential overflows", {
  A <- matrix(1e6, 6, 6)
  expect_true(is.na(flow_map_similarity(A, diag(6), h = 1)))
  expect_false(is.na(flow_map_similarity(A, diag(6), h = 1e-9)))
})

# --- the estimator ------------------------------------------------------------

## FM5 [oracle]: the plug-in at the pick equals the similarity of the fit's own
## unscaled coefficients (`m_ode_est`, `d_ode_est`, the package's other
## implementation of the per-power unscaling) taken to the standardized units
## by hand (the raw system, or at q = 3 each cell's Jacobian at the mean of
## its own observed frames, the design's `mean_state`, conjugated as
## S^-1 A S with S the diagonal of the design's response scale), at q = 1 and
## at q = 3; the path profile is the same hand computation at every point,
## starts at 1 (both systems empty at lambda_max), and the estimate is its
## value at the pick. At q = 3 no `x_ref` is passed: the default is the
## design's own means, an explicit list of the same means agrees, and the
## reading at the mother's mean for both cells (until 0.12.0) differs.
test_that("FM5: inheritance_flow() is the standardized plug-in at the pick: it matches the fit's unscaled coefficients conjugated by hand, at q = 1 and q = 3", {
  cases <- list(q1 = list(fit = .flow_fit_small(), raw = .fixture_dataset_small()$raw),
                q3 = list(fit = .flow_fit_q3(), raw = .flow_dataset_q3()$raw))
  for(nm in names(cases)){
    fit <- cases[[nm]]$fit
    raw <- cases[[nm]]$raw
    tn <- fit$tuning
    q <- fit$design$q
    p <- fit$design$p
    s_vec <- fit$design$response_scale
    # the design stores each cell's mean observed state in observation units
    mean_state <- list(M = colMeans(raw$Y_M), D = colMeans(raw$Y_D))
    expect_equal(fit$design$mean_state, mean_state, info = nm)
    expect_equal(tn$path$design$mean_state, mean_state, info = nm)
    x_ref <- if(q > 1) mean_state else list(M = NULL, D = NULL)
    out <- inheritance_flow(fit, bool_path = TRUE)
    h <- stats::median(diff(cases[[nm]]$raw$m_obs_time))
    expect_equal(out$h, h, info = nm)
    L <- length(tn$lambda_path)
    expect_equal(nrow(out$path_df), L, info = nm)
    expect_equal(out$path_df$lambda, tn$lambda_path, info = nm)
    expect_equal(out$pick_idx, tn$pick_idx, info = nm)
    expect_equal(out$path_df$bool_pick, seq_len(L) == tn$pick_idx, info = nm)
    expect_equal(out$path_df$similarity[1], 1, info = nm)
    to_std <- function(A) diag(1/s_vec) %*% A %*% diag(s_vec)
    # the pick: the fit's own read-out in raw units, conjugated to standardized
    at_pick <- flow_map_similarity(to_std(.jacobian_at(fit$m_ode_est, q, x_ref$M)),
                                   to_std(.jacobian_at(fit$d_ode_est, q, x_ref$D)), h = h)
    expect_equal(out$path_df$similarity[tn$pick_idx], at_pick, tolerance = 1e-10, info = nm)
    expect_equal(out$estimate, at_pick, tolerance = 1e-10, info = nm)
    # every path point: the hand unscaling of the tuned path, then the same conjugation
    hand <- vapply(seq_len(L), function(l){
      theta_M <- .unscale_coef(matrix(tn$path$theta_M[, , l], p), s_vec, q)
      theta_D <- .unscale_coef(matrix(tn$path$theta_D[, , l], p), s_vec, q)
      flow_map_similarity(to_std(.jacobian_at(theta_M, q, x_ref$M)), to_std(.jacobian_at(theta_D, q, x_ref$D)), h = h)
    }, numeric(1))
    expect_equal(out$path_df$similarity, hand, tolerance = 1e-10, info = nm)
    expect_equal(inheritance_flow(fit), out$estimate, info = nm)
    expect_true(out$estimate >= 0 && out$estimate <= 1, info = nm)
    # the same from the tuning object alone
    expect_equal(inheritance_flow(tn), out$estimate, info = nm)
    if(q > 1){
      # an explicit x_ref equal to the design's means is the default; the
      # mother's mean for both cells is another number
      expect_equal(inheritance_flow(fit, x_ref = mean_state), out$estimate, info = nm)
      at_mother <- inheritance_flow(fit, x_ref = list(M = mean_state$M, D = mean_state$M))
      expect_false(isTRUE(all.equal(at_mother, out$estimate)), info = nm)
      expect_equal(at_mother,
                   flow_map_similarity(to_std(.jacobian_at(fit$m_ode_est, q, mean_state$M)),
                                       to_std(.jacobian_at(fit$d_ode_est, q, mean_state$M)), h = h),
                   tolerance = 1e-10, info = nm)
    }
  }
})

## FM6 [invariant]: a daughter who is a copy of her mother (the same frames
## and observations) is read as a clone: the two fitted systems coincide at
## every path point, so the estimate is 1.
test_that("FM6: a daughter identical to her mother reads 1 at every path point", {
  raw <- .fixture_dataset("regression", num_d = 12)$raw
  public_same <- make_prnode_design(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$m_obs_time,
                                    Y_D = raw$Y_M, times_e = raw$times_e, q = 1)
  fit_same <- .with_warnings(fit_prnode(public_same))$value
  out <- inheritance_flow(fit_same, bool_path = TRUE)
  expect_equal(out$estimate, 1, tolerance = 1e-6)
  expect_true(all(abs(out$path_df$similarity[seq_len(out$pick_idx)] - 1) < 1e-6))
})

## FM7: the horizon argument: a user `h` replaces the default, the retired
## window argument is gone, and a non-positive `h` errors.
test_that("FM7: h is honoured and the window argument is gone", {
  fit <- .flow_fit_small()
  tn <- fit$tuning
  out_default <- inheritance_flow(fit, bool_path = TRUE)
  out_h <- inheritance_flow(fit, bool_path = TRUE, h = 0.5)
  expect_equal(out_h$h, 0.5)
  expect_false(isTRUE(all.equal(out_h$path_df$similarity[tn$pick_idx], out_default$path_df$similarity[tn$pick_idx])))
  expect_equal(out_h$estimate, out_h$path_df$similarity[tn$pick_idx])
  expect_error(inheritance_flow(fit, num_past_anchor = 8), "unused argument")
  expect_error(inheritance_flow(fit, h = -1))
})

## FM8: the entry guards: a fit without a tuning and a plain path are
## refused with a message naming the remedy; at q = 3 a malformed `x_ref`
## (a bare vector, the single state of 0.12.0; a list missing a cell; a
## wrong length) is refused, and a design without `mean_state` (built by
## hand or before 0.13.0) is refused unless `x_ref` supplies the states.
test_that("FM8: inheritance_flow() refuses a fit without a tuning, a non-fit, a malformed x_ref and a q = 3 design without mean_state", {
  fit <- .flow_fit_small()
  plain <- .with_warnings(fit_prnode(.fixture_dataset_small()$public, lambda = fit$tuning$lambda))$value
  expect_null(plain$tuning)
  expect_error(inheritance_flow(plain), "tuning")
  expect_error(inheritance_flow(fit$path), "prnode_fit")
  expect_error(inheritance_flow(list()), "prnode_fit")
  fit_q3 <- .flow_fit_q3()
  expect_true(is.finite(inheritance_flow(fit_q3)))
  expect_error(inheritance_flow(fit_q3, x_ref = 1:4), "x_ref")
  expect_error(inheritance_flow(fit_q3, x_ref = list(M = rep(1, 4))), "x_ref")
  expect_error(inheritance_flow(fit_q3, x_ref = list(M = rep(1, 4), D = 1:3)), "x_ref")
  fit_bare <- fit_q3
  fit_bare$tuning$path$design$mean_state <- NULL
  expect_error(inheritance_flow(fit_bare), "mean_state")
  expect_equal(inheritance_flow(fit_bare, x_ref = fit_q3$design$mean_state), inheritance_flow(fit_q3))
  # at q = 1 no state is read, so a stripped design still works
  fit_bare1 <- fit
  fit_bare1$tuning$path$design$mean_state <- NULL
  expect_equal(inheritance_flow(fit_bare1), inheritance_flow(fit))
})

# --- the p-value and the test -------------------------------------------------

## FM9 [oracle]: the permutation p-value counted by hand, its range
## {1/(K+1), ..., 1}, ties counted as at-or-above, and its guards.
test_that("FM9: inheritance_pvalue() is (1 + #{null >= estimate}) / (1 + K)", {
  expect_equal(inheritance_pvalue(0.8, c(0.5, 0.6, 0.9, 0.4)), 2/5)
  expect_equal(inheritance_pvalue(0.95, c(0.5, 0.6, 0.9, 0.4)), 1/5)
  expect_equal(inheritance_pvalue(0.1, c(0.5, 0.6, 0.9, 0.4)), 1)
  expect_equal(inheritance_pvalue(0.6, c(0.5, 0.6, 0.9, 0.4)), 3/5)
  set.seed(9)
  for(K in c(1, 7, 24)){
    null_vec <- stats::runif(K)
    est <- stats::runif(1)
    pv <- inheritance_pvalue(est, null_vec)
    expect_equal(pv, (1 + sum(null_vec >= est))/(1 + K), info = paste("K", K))
    expect_true(pv >= 1/(1 + K) && pv <= 1, info = paste("K", K))
  }
  expect_error(inheritance_pvalue(c(0.1, 0.2), 0.5))
  expect_error(inheritance_pvalue(0.5, numeric(0)))
  expect_error(inheritance_pvalue(0.5, c(0.2, NA)))
})

## FM10: inheritance_test() on the small dataset with two wrong mothers (the
## mothers of two other seeds): the true pair's estimate is
## inheritance_flow() of its fit, each null estimate is inheritance_flow()
## of the wrong pairing's own tuned fit (rebuilt here by hand, at the true
## pair's fixed spar so that the rule is the same for every pairing), the
## p-value is inheritance_pvalue() of the two, a passed `fit` gives the same
## numbers as fitting inside, and the matrix and list forms of a wrong
## mother agree. The test's exchangeability: pairing the daughter with her
## own mother as a "wrong" mother ties the estimate and is counted.
test_that("FM10: inheritance_test() refits every wrong pairing at the tuning rule and reads it by inheritance_flow()", {
  ds <- .fixture_dataset_small()
  raw <- ds$raw
  fit <- .flow_fit_small()
  others <- lapply(2:3, function(s){
    build_md_dataset(seed = s, num_d = 18, noise_level_ratio = 0.1, truth_args = list(p = 4),
                     obs_time = seq(0, 1, length.out = 36), steps_per_frame = 100)$raw$Y_M
  })
  names(others) <- c("seed2", "seed3")
  # the fixture's design was built on raw$times_e (100 steps per frame), so
  # every refit is on that grid: a derived grid changes the integral by 1e-7
  test <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = others,
                                          fit = fit, bool_keep_fits = TRUE, times_e = raw$times_e,
                                          design_args = list(spar = 0.5)))$value
  expect_s3_class(test, "prnode_inheritance_test")
  expect_equal(names(test), c("estimate", "fit", "fits", "h", "null", "num_null", "p_value"))
  expect_equal(test$estimate, inheritance_flow(fit))
  expect_equal(test$num_null, 2)
  expect_identical(names(test$null), c("seed2", "seed3"))
  expect_equal(test$p_value, inheritance_pvalue(test$estimate, test$null))
  # each null estimate is the wrong pairing's own tuned fit, read the same way
  for(i in 1:2){
    public_i <- make_prnode_design(m_obs_time = raw$m_obs_time, Y_M = others[[i]], d_obs_time = raw$d_obs_time,
                                   Y_D = raw$Y_D, q = 1, times_e = raw$times_e, spar = 0.5)
    fit_i <- .with_warnings(fit_prnode(public_i))$value
    expect_equal(unname(test$null[i]), inheritance_flow(fit_i), info = names(others)[i])
    expect_equal(test$fits[[i]]$lambda, fit_i$lambda, info = names(others)[i])
  }
  # fitting the true pair inside gives the same numbers as passing `fit`
  test_inside <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                                 null_mothers = others[1], times_e = raw$times_e,
                                                 design_args = list(spar = 0.5)))$value
  expect_equal(test_inside$estimate, test$estimate)
  expect_equal(unname(test_inside$null), unname(test$null[1]))
  expect_null(test_inside$fits)
  # the list form of a wrong mother
  test_list <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                               null_mothers = list(list(m_obs_time = raw$m_obs_time, Y_M = others[[1]])),
                                               fit = fit, times_e = raw$times_e,
                                               design_args = list(spar = 0.5)))$value
  expect_equal(unname(test_list$null), unname(test$null[1]))
  expect_identical(names(test_list$null), "null_1")
  # the own mother as a wrong mother ties and is counted at or above
  test_own <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                              null_mothers = list(own = raw$Y_M), fit = fit,
                                              times_e = raw$times_e, design_args = list(spar = 0.5)))$value
  expect_equal(unname(test_own$null), test_own$estimate)
  expect_equal(test_own$p_value, 1)
  # guards
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list()), "null_mothers")
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(others[[1]][, 1:2])),
               "null_mothers")
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = others,
                                fit = fit, mu = c(S = 1, M = 0.9, D = 0.9)), "mu")
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = others,
                                fit = fit, design_args = list(spar = 0.25)), "spar")
  expect_output(print(test), "p-value")
})

## FM11: power on a clone: a daughter who is a copy of the mother's frames,
## against three unrelated wrong mothers, reaches the smallest attainable
## p-value 1/4 (she reads 1 with her own mother, FM6, and below 1 with any
## other).
test_that("FM11: a clone daughter against three unrelated mothers gets the smallest attainable p-value", {
  ds <- .fixture_dataset("regression", num_d = 12)
  raw <- ds$raw
  others <- lapply(1:3, function(s){
    build_md_dataset(seed = 200 + s, num_d = 12, noise_level_ratio = 0.1, truth_args = .fx_truth_args,
                     obs_time = .fx_obs_time, times_e = .fx_times_coarse, level = "regression")$raw$Y_M
  })
  test <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$m_obs_time, raw$Y_M, null_mothers = others,
                                          times_e = raw$times_e, design_args = list(spar = 0.5)))$value
  expect_equal(test$estimate, 1, tolerance = 1e-6)
  expect_true(all(test$null < 1 - 1e-3))
  expect_equal(test$p_value, 1/4)
})

## FM12: at q = 3 the test linearizes every pairing's two fields at their own
## cells' mean observed states (the design's `mean_state`, so the true pair's
## estimate is inheritance_flow(fit)) unless `x_ref` fixes two states for all
## pairings; with a fixed `x_ref` the true pair's estimate is
## inheritance_flow(fit, x_ref = ), and the single state of 0.12.0 is refused.
test_that("FM12: inheritance_test() at q = 3 reads each pairing at its own means by default and honours a fixed list", {
  ds <- .flow_dataset_q3()
  raw <- ds$raw
  fit <- .flow_fit_q3()
  # seed 6: seeds 4 and 5 draw an unstable cubic system (deSolve stops)
  other <- build_md_dataset(seed = 6, num_d = 18, noise_level_ratio = 0.1,
                            truth_args = list(p = 4, q = 3, num_shared = 4, num_mother_only = 1,
                                              num_daughter_only = 1),
                            obs_time = seq(0, 1, length.out = 36), steps_per_frame = 200, q = 3,
                            level = "ode")$raw$Y_M
  x_fixed <- list(M = rep(1, 4), D = rep(1, 4))
  test_fixed <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                                null_mothers = list(other), fit = fit, q = 3, x_ref = x_fixed,
                                                times_e = raw$times_e, design_args = list(spar = 0.5)))$value
  expect_equal(test_fixed$estimate, inheritance_flow(fit, x_ref = x_fixed))
  test_own <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                              null_mothers = list(other), fit = fit, q = 3,
                                              times_e = raw$times_e, design_args = list(spar = 0.5)))$value
  expect_equal(test_own$estimate, inheritance_flow(fit))
  expect_equal(test_own$estimate, inheritance_flow(fit, x_ref = list(M = colMeans(raw$Y_M), D = colMeans(raw$Y_D))))
  expect_false(isTRUE(all.equal(test_own$estimate, test_fixed$estimate)))
  expect_true(test_own$p_value %in% c(0.5, 1))
  # the single state of 0.12.0 is refused before any wrong pairing is fitted
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(other),
                                fit = fit, q = 3, x_ref = rep(1, 4), times_e = raw$times_e,
                                design_args = list(spar = 0.5)), "x_ref")
  # a passed fit must have been built on the data passed: the true pair is
  # read at the fit's stored means, the wrong pairings at the data's, so a
  # rescaled Y_M of the same dimensions is refused (with a fixed x_ref the
  # states are the caller's and the check does not apply)
  expect_error(inheritance_test(raw$m_obs_time, 2*raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(other),
                                fit = fit, q = 3, times_e = raw$times_e, design_args = list(spar = 0.5)),
               "mean_state")
})

# --- the plot -----------------------------------------------------------------

## FM13 (replaces L14): plot_md_inheritance() returns inheritance_flow()'s
## numbers; `truth` sets the reference to the truth's standardized similarity
## (`truth$inheritance_std`) at the same horizon and widens the y-range to
## include it.
test_that("FM13: plot_md_inheritance() returns inheritance_flow()'s path and the truth's similarity", {
  ds <- .fixture_dataset_small()
  fit <- .flow_fit_small()
  flow <- inheritance_flow(fit, bool_path = TRUE)
  out <- .on_null_device(plot_md_inheritance(fit))
  expect_equal(names(out), c("estimate", "path_df", "reference", "pick_idx", "ylim"))
  expect_equal(out$path_df, flow$path_df)
  expect_equal(out$estimate, flow$estimate)
  expect_equal(out$pick_idx, flow$pick_idx)
  expect_true(is.na(out$reference))
  out_truth <- .on_null_device(plot_md_inheritance(fit, truth = ds$truth))
  expect_equal(out_truth$reference, ds$truth$inheritance_std, tolerance = 1e-12)
  expect_true(out_truth$ylim[1] <= out_truth$reference && out_truth$ylim[2] >= out_truth$reference)
  expect_equal(out_truth$path_df, out$path_df)
  # at q = 3 the plot needs no x_ref: the design's means, and the truth's
  # reference is the dataset's own standardized similarity at its own means
  out_q3 <- .on_null_device(plot_md_inheritance(.flow_fit_q3(), truth = .flow_dataset_q3()$truth))
  expect_true(is.finite(out_q3$reference))
  expect_equal(out_q3$reference, .flow_dataset_q3()$truth$inheritance_std, tolerance = 1e-12)
  expect_equal(out_q3$estimate, inheritance_flow(.flow_fit_q3()))
  expect_error(.on_null_device(plot_md_inheritance(fit, truth = list(p = 4))), "m_ode_pars")
  expect_error(.on_null_device(plot_md_inheritance(list())), "prnode")
})

## FM14 (review fixes of 2026-10-28): with `fit` given and no `times_e`, the
## wrong pairings are built on the fit's own grid, so the own mother as a
## wrong mother ties exactly (before the fix the derived grid moved the
## estimate by 1e-7 and the tie could read as "below"); `design_args` and
## `fit_args` reach every pairing (a refit with `spar = 0.7` by hand equals
## the null estimate) and are checked for names; a fit of a subset of
## responses is refused by inheritance_flow().
test_that("FM14: inheritance_test() takes the fit's grid by default, passes design_args and fit_args to every pairing, and inheritance_flow() guards the square system", {
  ds <- .fixture_dataset_small()
  raw <- ds$raw
  fit <- .flow_fit_small()
  test_own <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                              null_mothers = list(own = raw$Y_M), fit = fit,
                                              design_args = list(spar = 0.5)))$value
  expect_identical(unname(test_own$null), test_own$estimate)
  expect_equal(test_own$p_value, 1)
  # design_args and fit_args reach the wrong pairing's design and fit
  other <- build_md_dataset(seed = 2, num_d = 18, noise_level_ratio = 0.1, truth_args = list(p = 4),
                            obs_time = seq(0, 1, length.out = 36), steps_per_frame = 100)$raw$Y_M
  test_spar <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                               null_mothers = list(seed2 = other), times_e = raw$times_e,
                                               design_args = list(spar = 0.7), fit_args = list(edge_tol = 1e-6),
                                               bool_keep_fits = TRUE))$value
  public_hand <- make_prnode_design(m_obs_time = raw$m_obs_time, Y_M = other, d_obs_time = raw$d_obs_time,
                                    Y_D = raw$Y_D, q = 1, times_e = raw$times_e, spar = 0.7)
  fit_hand <- .with_warnings(fit_prnode(public_hand, edge_tol = 1e-6))$value
  expect_equal(unname(test_spar$null), inheritance_flow(fit_hand))
  expect_equal(test_spar$fits[[1]]$tuning$path$edge_tol, 1e-6)
  expect_equal(test_spar$fit$tuning$path$edge_tol, 1e-6)
  public_default <- make_prnode_design(m_obs_time = raw$m_obs_time, Y_M = other, d_obs_time = raw$d_obs_time,
                                       Y_D = raw$Y_D, q = 1, times_e = raw$times_e)
  fit_default <- .with_warnings(fit_prnode(public_default))$value
  expect_false(isTRUE(all.equal(inheritance_flow(fit_default), inheritance_flow(fit_hand))))
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(other),
                                design_args = list(0.7)), "named")
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(other),
                                design_args = list(q = 2)), "design_args")
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(other),
                                fit_args = list(lambda = 0.1)), "fit_args")
  # a subset of responses is not a square system
  public_sub <- ds$public
  public_sub$y <- public_sub$y[, 1:2, drop = FALSE]
  public_sub$response_scale <- public_sub$response_scale[1:2]
  fit_sub <- .with_warnings(fit_prnode(public_sub))$value
  expect_error(inheritance_flow(fit_sub), "every protein as a response")
})

## FM15 (0.9.0): with no `design_args$spar`, every pairing chooses its own
## smoothing level and penalty (tune_prnode_spar()): the null estimate is
## inheritance_flow() of fit_prnode() on the wrong pairing's raw observations
## by that rule, and the true pair fitted inside is read the same way.
test_that("FM15: inheritance_test() without a fixed spar refits every pairing by the spar search", {
  ds <- .fixture_dataset_small()
  raw <- ds$raw
  other <- build_md_dataset(seed = 2, num_d = 18, noise_level_ratio = 0.1, truth_args = list(p = 4),
                            obs_time = seq(0, 1, length.out = 36), steps_per_frame = 100)$raw$Y_M
  test <- .with_warnings(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D,
                                          null_mothers = list(seed2 = other), times_e = raw$times_e,
                                          bool_keep_fits = TRUE))$value
  fit_true <- .with_warnings(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                                        Y_D = raw$Y_D, times_e = raw$times_e))$value
  fit_null <- .with_warnings(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = other, d_obs_time = raw$d_obs_time,
                                        Y_D = raw$Y_D, times_e = raw$times_e))$value
  expect_equal(test$estimate, inheritance_flow(fit_true))
  expect_equal(unname(test$null), inheritance_flow(fit_null))
  expect_equal(test$fit$tuning$spar, fit_true$tuning$spar)
  expect_true(test$fits[[1]]$tuning$spar %in% c(0, 0.25, 0.5, 0.75))
})

## FM16: the rule of the true pair's fit and of the wrong pairings must agree:
## a fit from the spar search with a fixed `design_args$spar`, or a fit at a
## fixed spar with none, would make the null come from another rule.
test_that("FM16: inheritance_test() refuses a fit whose spar rule differs from design_args", {
  raw <- .fixture_dataset_small()$raw
  other <- raw$Y_M
  fit_searched <- .with_warnings(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                                            Y_D = raw$Y_D, times_e = raw$times_e, spar_grid = c(0.25, 0.5)))$value
  fit_fixed <- .with_warnings(fit_prnode(m_obs_time = raw$m_obs_time, Y_M = raw$Y_M, d_obs_time = raw$d_obs_time,
                                         Y_D = raw$Y_D, times_e = raw$times_e, spar = 0.5))$value
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(a = other),
                                times_e = raw$times_e, fit = fit_searched, design_args = list(spar = 0.5)),
               "search")
  expect_error(inheritance_test(raw$m_obs_time, raw$Y_M, raw$d_obs_time, raw$Y_D, null_mothers = list(a = other),
                                times_e = raw$times_e, fit = fit_fixed), "NULL")
})

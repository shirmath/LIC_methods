# Plan section 14: end-to-end recovery (E1, E2) and the scorer (E4, E5); E3
# left with the measure it tested in package 0.7.0. E6, the leak regression against the retired
# four-block loss on the generator, left with the four-block helper on
# 2026-10-22; F6(a) in test-fit_prnode.R keeps the leak record on a
# hand-built fixture.

## E1: recovery floors along a path at the ODE level, with the regression
## level as the higher reference (Z16).
test_that("E1: along a 30-point path the ODE-level fit reaches MCC >= 0.6 in both cells and starts empty; regression >= 0.8", {
  .skip_unless_built("fit_prnode", "public")
  mu <- c(S = 1, M = 0.9, D = 0.9)
  for(level in c("ode", "regression")){
    ds <- .fixture_dataset(level, num_d = 36)
    lambda_path <- .lambda_path(ds$public, mu = mu, num_lambda = 30)
    fit <- fit_prnode(ds$public, lambda = lambda_path, mu = mu)
    mcc_d <- sapply(seq_along(lambda_path), function(l) .mcc(fit$d_network_est[, , l], ds$truth$d_network))
    mcc_m <- sapply(seq_along(lambda_path), function(l) .mcc(fit$m_network_est[, , l], ds$truth$m_network))
    floor_val <- if(level == "ode") 0.6 else 0.8
    expect_gte(max(mcc_d), floor_val, label = paste(level, "daughter"))
    expect_gte(max(mcc_m), floor_val, label = paste(level, "mother"))
    expect_false(any(.as_adj(fit$m_network_est[, , 1])), info = level)
    expect_false(any(.as_adj(fit$d_network_est[, , 1])), info = level)
  }
})

## E2 (Proposition 5 of the derivation notes `prop:split`): identical networks on noiseless regression-level data
## are recovered as identical: at the daughter-oracle lambda the two networks
## coincide and the deviations are exactly zero. The proposition's scalar
## condition, a_M / a_S <= mu for every regulator, is reported with any
## failure but is not sufficient here (review pass 6): with six regulators on
## different windows the fit at mu = 0.9 and at 0.95 has d_M != 0 and drops a
## mother edge from the daughter although every ratio is at most 0.924 (the
## KKT certificate at 1e-10, so this is the optimum, not the solver); from
## mu = 1 on, d = 0 exactly and the networks coincide (checked at 1, 1.2, 1.5,
## 2). The test runs at mu = 1; the condition the proposition needs in the
## multi-regulator case is an open item for the derivation notes.
test_that("E2: with d_delta = 0 and no noise, the fused fit returns equal networks and d = 0 at the oracle lambda", {
  .skip_unless_new_pipeline()
  .skip_unless_built("fit_prnode", "public")
  ds <- build_md_dataset(seed = 22, num_d = 36, noise_level_ratio = 0,
                         truth_args = list(p = 6, q = 1, num_shared = 6, num_mother_only = 0, num_daughter_only = 0,
                                           frac_unequal = 0),
                         obs_time = .fx_obs_time, times_e = .fx_times_coarse, level = "regression")
  expect_true(all(ds$truth$d_delta == 0))
  public <- ds$public
  mu <- c(S = 1, M = 1, D = 1)
  lambda_path <- .lambda_path(public, mu = mu, num_lambda = 30, num_decades = 3)
  fit <- fit_prnode(public, lambda = lambda_path, mu = mu)
  l <- .oracle_index(fit, ds$truth)
  p <- public$p
  ratio_vec <- sapply(seq_len(p), function(k){
    a_M <- sqrt(sum(public$x[, p + k]^2))
    a_S <- sqrt(sum(public$x[, k]^2))
    a_M/a_S
  })
  info <- paste("a_M / a_S per regulator:", paste(round(ratio_vec, 3), collapse = " "))
  expect_identical(.as_adj(fit$d_network_est[, , l]), .as_adj(fit$m_network_est[, , l]), info = info)
  expect_true(all(.path_mat(fit$path, "d_M", l) == 0), info = info)
  expect_true(all(.path_mat(fit$path, "d_D", l) == 0), info = info)
  expect_gte(.mcc(fit$d_network_est[, , l], ds$truth$d_network), 0.8)
})

## E3 left in package 0.7.0; the inheritance measure is tested in
## test-inheritance_flow.R.

## E4 (Z17, D0d(a)): the scorer against hand-built confusion matrices.
test_that("E4: score_md_networks() gives MCC 1 on the truth, 0 on the empty and the complete network, and counts by edge type", {
  .skip_unless_built("score_md_networks")
  .skip_unless_built("sim_md_network")
  truth <- sim_md_network(6, 1, num_shared = 6, num_mother_only = 2, num_daughter_only = 2, frac_unequal = 0.5,
                          seed_number = 44)
  perfect <- list(m_network_est = truth$m_network, d_network_est = truth$d_network, s_network_est = truth$s_network)
  sc <- score_md_networks(perfect, truth)
  expect_equal(unname(sc$mcc[c("M", "D", "S")]), c(1, 1, 1))
  expect_equal(unname(sc$f1[c("M", "D", "S")]), c(1, 1, 1))
  expect_equal(unname(sc$leak_rate[c("M", "D")]), c(0, 0))
  expect_equal(unname(sc$borrow_rate[c("M", "D")]), c(1, 1))
  conf_D <- sc$confusion$D
  expect_true(all(c("TP", "FP", "FN", "TN") %in% colnames(conf_D)))
  expect_true(all(c("shared_equal", "shared_unequal", "mother_only", "daughter_only", "absent") %in% rownames(conf_D)))
  expect_equal(sum(conf_D[, c("FP", "FN")]), 0)
  expect_equal(conf_D["daughter_only", "TP"], 2)
  expect_equal(conf_D["mother_only", "TN"], 2)
  expect_equal(sum(conf_D[, c("TP", "TN")]), 30)
  empty <- list(m_network_est = matrix(FALSE, 6, 6), d_network_est = matrix(FALSE, 6, 6), s_network_est = matrix(FALSE, 6, 6))
  sc0 <- score_md_networks(empty, truth)
  expect_equal(unname(sc0$mcc[c("M", "D", "S")]), c(0, 0, 0))
  expect_equal(unname(sc0$leak_rate[c("M", "D")]), c(0, 0))
  expect_equal(sc0$confusion$D["daughter_only", "FN"], 2)
  # the complete network has no true negatives, so the MCC denominator
  # vanishes and the scorer returns 0 (F5's former helper check, now on the
  # package)
  complete <- list(m_network_est = matrix(TRUE, 6, 6), d_network_est = matrix(TRUE, 6, 6), s_network_est = matrix(TRUE, 6, 6))
  sc1 <- score_md_networks(complete, truth)
  expect_equal(unname(sc1$mcc[c("M", "D", "S")]), c(0, 0, 0))
  expect_equal(unname(sc1$leak_rate[c("M", "D")]), c(1, 1))
  with_diag <- score_md_networks(perfect, truth, bool_score_diagonal = TRUE)
  expect_equal(sum(with_diag$confusion$D[, c("TP", "TN", "FP", "FN")]), 36)
  expect_equal(sum(sc$confusion$D[, c("TP", "TN", "FP", "FN")]), 30)
})

## E5 (section 8.6 of the derivation notes): the scalar MCC hides the leak.
test_that("E5: a daughter estimate equal to the mother's network has high MCC_D and a leak rate of 1", {
  .skip_unless_built("score_md_networks")
  .skip_unless_built("sim_md_network")
  truth <- sim_md_network(6, 1, num_shared = 6, num_mother_only = 2, num_daughter_only = 2, frac_unequal = 0.5,
                          seed_number = 45)
  copied <- list(m_network_est = truth$m_network, d_network_est = truth$m_network, s_network_est = truth$s_network)
  sc <- score_md_networks(copied, truth)
  expect_gt(sc$mcc[["D"]], 0.6)
  expect_equal(sc$leak_rate[["D"]], 1)
  expect_equal(sc$leak_rate[["M"]], 0)
  expect_equal(.mcc(truth$m_network, truth$d_network), sc$mcc[["D"]], tolerance = 1e-12)
})

# Plan, E1, restated for the
# known-variance BIC rule of 0.9.0: the package's criterion and pick against
# plain lm.fit() on its own supports and against a glmnet path, with the
# helpers in helper-tuning.R. glmnet is used here only; the file skips
# without it.

## E1 (ii): the criterion recomputed by oracle_kv_bic() on the package's own
## supports shares nothing with tune_prnode() but lm.fit(), so it must agree
## to rounding, on both fixtures, at the defaults and with mu moved.
test_that("E1 (ii): the criterion on the package's supports equals the oracle's, and the pick is its first minimum", {
  skip_if_not_installed("glmnet")
  small <- .fixture_dataset_small()$public
  mu_09 <- c(S = 1, M = 0.9, D = 0.9)
  cases <- list(
    list(label = "p4 defaults 80/4", d = small, tn = .fixture_tuning_default_small()),
    list(label = "p4 20/2", d = small, tn = .fixture_tuning_small()),
    list(label = "p6 40/3", d = .tuning_public(), tn = .fixture_tuning()),
    list(label = "p4 mu = (1, 0.9, 0.9)", d = small, tn = .with_warnings(tune_prnode(small, mu = mu_09))$value))
  for(case in cases){
    d <- case$d
    tn <- case$tn
    pkg_path <- lapply(seq_len(d$p), function(j) .path_coef_mat(tn$path, j))
    crit <- oracle_kv_bic(d, pkg_path)
    expect_equal(crit, tn$criterion, tolerance = 1e-10, info = case$label)
    expect_identical(which.min(crit), tn$pick_idx, info = case$label)
  }
})

## E1 (i): the same rule on a glmnet path. The two solvers agree exactly only in
## the limit of their tolerances, so their supports differ on a fraction of a
## percent of the (column, lambda) entries and the picks can differ by a step
## on a coarse grid; the comparison is made on the 80-point, 4-decade grid.
test_that("E1 (i): on the 80-point, 4-decade grid the glmnet path gives the same pick as the package's", {
  skip_if_not_installed("glmnet")
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_default_small()
  expect_length(tn$lambda_path, 80)
  expect_equal(tn$lambda_path[1], lambda_max_prnode(small, mu = c(S = 1, M = 1, D = 1), bool_pooled = TRUE),
               tolerance = 1e-12)
  glmnet_fit <- oracle_mu_path(small, c(S = 1, M = 1, D = 1), tn$lambda_path)
  crit <- oracle_kv_bic(small, glmnet_fit$path)
  expect_lte(abs(which.min(crit) - tn$pick_idx), 1)
  pkg_path <- lapply(seq_len(small$p), function(j) .path_coef_mat(tn$path, j))
  agreement <- mean(unlist(lapply(seq_len(small$p), function(j){
    (glmnet_fit$path[[j]] != 0) == (pkg_path[[j]] != 0)
  })))
  expect_gt(agreement, 0.99)
})

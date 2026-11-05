# Plan on tune_prnode()
# itself, restated for the known-variance BIC rule of package 0.9.0: the path
# (§1, T1-T8 and Z1 (c)'s E8) and the criterion along the path (§3, C4-C7).
# The anchor's autocorrelation correction and the cross-validation prune (the
# plan's V tests) left with the rule. The helpers' unit tests are in
# test-tune_prnode_helpers.R, the fit_prnode() branches, printing and
# invariances in test-tune_prnode_fit.R, and the glmnet oracle in
# test-tune_prnode_oracle.R.

# --- §1 the path --------------------------------------------------------------------

## T1: the grid is num_lambda points log-spaced over num_decades decades below
## the pooled lambda_max, and the chosen lambda is one of them.
test_that("T1: the default grid runs from the pooled lambda_max over num_decades decades; lambda is on the path", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  L <- 40
  expect_length(tn$lambda_path, L)
  expect_true(all(diff(tn$lambda_path) < 0))
  expect_equal(tn$lambda_path[1], lambda_max_prnode(public, mu = c(S = 1, M = 1, D = 1), bool_pooled = TRUE),
               tolerance = 1e-12)
  expect_equal(tn$lambda_path[L], tn$lambda_path[1]*10^-3, tolerance = 1e-12)
  expect_identical(tn$lambda, tn$lambda_path[tn$pick_idx])
  expect_identical(tn$path$lambda, tn$lambda_path)

  small <- .fixture_dataset_small()$public
  two <- .with_warnings(tune_prnode(small, num_lambda = 2, num_decades = 0.5))$value
  expect_length(two$lambda_path, 2)
  expect_equal(two$lambda_path[1], lambda_max_prnode(small, bool_pooled = TRUE), tolerance = 1e-12)
  expect_equal(two$lambda_path[2], two$lambda_path[1]*10^-0.5, tolerance = 1e-12)
  expect_true(two$pick_idx %in% 1:2)
})

## T2: the first point is the pooled lambda_max, where every protein's KKT
## condition holds at zero (the H-series pins lambda_max_prnode()); this is
## what makes the two cap definitions of Z4 agree on the default grid (V5).
test_that("T2: the first path point is the null model for every protein", {
  for(tn in list(.fixture_tuning(), .fixture_tuning_small())){
    p <- length(tn$rho)
    for(j in seq_len(p)){
      expect_true(all(.path_coef_mat(tn$path, j)[, 1] == 0), info = paste("protein", j))
    }
    expect_true(any(.path_coef_mat(tn$path, 1)[, 2:length(tn$lambda_path)] != 0))
  }
})

## T3: a user path is used as given (any length, including one point); an
## increasing one reaches fit_prnode_path()'s sorting warning and the tuning
## is on the sorted path.
test_that("T3: a user lambda_path is used as given, including L = 1; an increasing one is sorted with the solver's warning", {
  small <- .fixture_dataset_small()$public
  lambda_max <- lambda_max_prnode(small, bool_pooled = TRUE)
  user <- lambda_max*c(0.7, 0.2, 0.05)
  tn <- .with_warnings(tune_prnode(small, lambda_path = user))$value
  expect_identical(tn$path$lambda, user)
  expect_identical(tn$lambda_path, user)
  expect_length(tn$criterion, 3)

  one <- .with_warnings(tune_prnode(small, lambda_path = lambda_max/10))$value
  expect_identical(one$lambda_path, lambda_max/10)
  expect_identical(one$pick_idx, 1L)
  expect_identical(one$lambda, lambda_max/10)

  expect_warning(inc <- tune_prnode(small, lambda_path = rev(user)), regexp = "sorting it in decreasing order")
  expect_identical(inc$lambda_path, user)
  expect_equal(inc$criterion, tn$criterion, tolerance = 1e-10)
  expect_identical(inc$pick_idx, tn$pick_idx)
})

## T4: the `path` argument. (a) no refit: the tuning of a pre-fitted path is the
## tuning that fitted it; (b) `mu` is read from the path, a disagreeing `mu`
## errors, an agreeing one (named or not) passes; (c) the arguments documented
## as ignored with `path` are ignored (plan row 36, 87); (d) a path of another
## design or a non-path errors.
test_that("T4: a pre-fitted path is tuned without refitting, carries its mu, ignores the solver arguments, and must match the design", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  from_path <- tune_prnode(public, path = tn$path)
  for(nm in c("criterion", "k", "rss", "pick_idx", "lambda", "lambda_path", "mu", "spar")){
    expect_identical(from_path[[nm]], tn[[nm]], info = nm)
  }
  expect_identical(from_path$path, tn$path)

  small <- .fixture_dataset_small()$public
  mu_09 <- c(S = 1, M = 0.9, D = 0.9)
  path_09 <- .with_warnings(fit_prnode_path(small, lambda = .lambda_path(small, mu = mu_09, num_lambda = 5,
                                                                          num_decades = 1), mu = mu_09))$value
  tn_09 <- tune_prnode(small, path = path_09)
  expect_identical(tn_09$mu, mu_09)
  expect_error(tune_prnode(small, path = path_09, mu = c(1, 1, 1)), regexp = "`path` was fitted at mu")
  expect_identical(tune_prnode(small, path = path_09, mu = c(1, 0.9, 0.9)), tn_09)
  expect_identical(tune_prnode(small, path = path_09, mu = mu_09), tn_09)

  ignored <- tune_prnode(public, path = tn$path, edge_tol = 1, max_iter = 1, tol = 1, lambda_path = c(5, 1),
                         num_lambda = 3, num_decades = 1)
  expect_identical(ignored, from_path)

  expect_error(tune_prnode(public, path = .fixture_tuning_small()$path), regexp = "does not match")
  expect_error(tune_prnode(small, path = .fixture_tuning_q2()$path), regexp = "does not match")
  expect_error(tune_prnode(public, path = list(lambda = 1)), regexp = "must be a prnode_path")
})

## T5(a): .path_coef_mat() stacks s, d_M, d_D in the design's column order. The
## solver's objective is computed from its coefficients in the orthonormalized
## groups, so the fixtures' .objective() on the stacked vector reproduces it
## only if the stacking is in the order `x` expects.
test_that("T5(a): the stacked coefficients reproduce the solver's objective on the original design", {
  small <- .fixture_dataset_small()
  tn <- .fixture_tuning_small()
  path <- tn$path
  mu <- c(S = 1, M = 1, D = 1)
  for(j in seq_len(4)) for(l in seq_along(path$lambda)){
    beta <- .path_coef_mat(path, j)[, l]
    expect_equal(.objective(small$public, j, beta, path$lambda[l], mu), path$objective[j, l],
                 tolerance = 1e-10, info = paste(j, l))
  }
  public <- .tuning_public()
  tn6 <- .fixture_tuning()
  for(j in c(1, 6)) for(l in c(1, 10, 25, 40)){
    beta <- .path_coef_mat(tn6$path, j)[, l]
    expect_equal(.objective(public, j, beta, tn6$path$lambda[l], mu), tn6$path$objective[j, l],
                 tolerance = 1e-10, info = paste(j, l))
  }
})

## T5(b): the `drop = FALSE` corners (plan row 43): pq = 1 with L > 1, and L = 1.
test_that("T5(b): .path_coef_mat() stays a matrix at pq = 1 and at L = 1", {
  ds <- .fixture_dataset_small()
  one <- .one_protein_public(ds, 1)
  path_1 <- .with_warnings(fit_prnode_path(one, lambda = .lambda_path(one, num_lambda = 3, num_decades = 1)))$value
  m <- .path_coef_mat(path_1, 1)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(3, 3))
  small <- ds$public
  path_L1 <- .with_warnings(fit_prnode_path(small, lambda = lambda_max_prnode(small, bool_pooled = TRUE)/5))$value
  m2 <- .path_coef_mat(path_L1, 2)
  expect_true(is.matrix(m2))
  expect_equal(dim(m2), c(12, 1))
  expect_equal(c(m2), c(path_L1$s[2, , 1], path_L1$d_M[2, , 1], path_L1$d_D[2, , 1]))
})

## T5(c) (Phase 4, Z3 (b)): one definition of "support" (plan row 79). The
## print method reports the active groups at the anchor and at the pick, the
## count .path_coef_mat() gives when each regulator's q columns of one block
## are one group; at q = 2 that is half the nonzero-column count the method
## printed before Phase 4.
test_that("T5(c): print() reports the active groups at the pick, counted as groups", {
  tn <- .fixture_tuning_small()
  p <- dim(tn$path$s)[1]
  groups_at <- function(tn, l) sum(vapply(seq_len(p), function(j) .active_groups(tn$path, j, l), numeric(1)))
  out <- paste(utils::capture.output(print(tn)), collapse = "\n")
  expect_match(out, paste0("active groups[ =:]+", groups_at(tn, tn$pick_idx), "\\b"))
  expect_equal(tn$k[tn$pick_idx], groups_at(tn, tn$pick_idx))
  tq2 <- .fixture_tuning_q2()
  out_q2 <- paste(utils::capture.output(print(tq2)), collapse = "\n")
  expect_match(out_q2, paste0("active groups[ =:]+", groups_at(tq2, tq2$pick_idx), "\\b"))
})

## T6: the solver arguments reach fit_prnode_path() when the path is fitted
## here (plan row 36). With one outer sweep the solver still runs its inner
## active-set sweeps, so a one-group solution can converge; the test asserts
## the warning names max_iter = 1, the null point converges, and more fits are
## left unconverged than at the default, not that every point fails.
test_that("T6: edge_tol, max_iter and tol reach fit_prnode_path() when the path is fitted here", {
  small <- .fixture_dataset_small()$public
  base <- .fixture_tuning_small()
  rough <- .with_warnings(tune_prnode(small, num_lambda = 20, num_decades = 2, max_iter = 1))
  expect_true(any(grepl("max_iter = 1 ", rough$warnings, fixed = TRUE)))
  converged <- rough$value$path$converged
  expect_true(all(converged[, 1]))
  expect_gt(sum(!converged), sum(!base$path$converged))

  loose <- .with_warnings(tune_prnode(small, num_lambda = 20, num_decades = 2, edge_tol = 10))$value
  expect_false(any(loose$path$networks$D))
  expect_true(any(base$path$networks$D))
  for(nm in c("criterion", "k", "rss", "pick_idx", "lambda", "spar")){
    expect_identical(loose[[nm]], base[[nm]], info = nm)
  }

  tol_set <- .with_warnings(tune_prnode(small, num_lambda = 20, num_decades = 2, tol = 1e-4))$value
  expect_equal(tol_set$path$tol, 1e-4)
  expect_equal(base$path$tol, 1e-10)
})

## T7: every protein shares the pooled path (the tuning study suggestion 5 rejects a
## per-protein lambda): restricting `public` to the protein with the smallest
## own lambda_max moves the path's first point down to that value.
test_that("T7: the path is shared by the proteins; a one-response public starts at its own lambda_max", {
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_small()
  expect_null(dim(tn$path$lambda))
  expect_length(tn$path$lambda, 20)
  per_protein <- lambda_max_prnode(small, bool_pooled = FALSE)
  expect_equal(tn$lambda_path[1], max(per_protein), tolerance = 1e-12)
  j <- which.min(per_protein)
  public_j <- small
  public_j$y <- small$y[, j, drop = FALSE]
  public_j$response_scale <- small$response_scale[j]
  public_j$protein_names <- small$protein_names[j]
  tn_j <- .with_warnings(tune_prnode(public_j, num_lambda = 2, num_decades = 1))$value
  expect_equal(tn_j$lambda_path[1], per_protein[j], tolerance = 1e-12)
  expect_lt(tn_j$lambda_path[1], tn$lambda_path[1])
  expect_length(tn_j$criterion, 2)
})

## T8 (Z1): non-convergence propagates as the solver's warning through
## tune_prnode() and fit_prnode(). At Phase 3 the examples' 20-point path
## left one (protein, lambda) at max_iter (protein 1, point 14) and the test
## asserted it; that non-convergence was the strong rule's re-admission gate
## (estimator plan R7b, fixed 2026-10-23), so the fixture
## now converges everywhere and the propagation is asserted through
## max_iter = 1 alone. Restated in place; recorded in both plans.
test_that("T8: the solver's max_iter warning propagates through tune_prnode() and fit_prnode(public), and `converged` records the fit", {
  tn <- .fixture_tuning_small()
  warnings <- .fixture_tuning_small_warnings()
  expect_false(any(grepl("max_iter", warnings)))
  expect_true(all(tn$path$converged))
  expect_equal(dim(tn$path$converged), c(4, 20))

  small <- .fixture_dataset_small()$public
  through_fit <- .with_warnings(fit_prnode(small, max_iter = 1))
  expect_true(any(grepl("max_iter = 1 ", through_fit$warnings, fixed = TRUE)))
  expect_true(any(!through_fit$value$tuning$path$converged))
})

## E8 (Phase 4, Z1 (c)): a second, specific warning from tune_prnode() when the
## anchor or the pick sits on a (protein, lambda) that did not converge,
## naming the index and the protein(s); absent when the non-converged points
## lie elsewhere. The non-converged path is passed pre-fitted, so the solver's
## own warning cannot fire and the new warning is isolated.
test_that("E8: tune_prnode() warns, naming index and protein, only when the pick did not converge", {
  small <- .fixture_dataset_small()$public
  rough <- .with_warnings(tune_prnode(small, num_lambda = 20, num_decades = 2, max_iter = 1))$value
  bad_pick <- which(!rough$path$converged[, rough$pick_idx])
  stopifnot(length(bad_pick) > 0)
  res <- .with_warnings(tune_prnode(small, path = rough$path))
  expect_length(res$warnings, 1)
  message <- paste(res$warnings, collapse = " ")
  expect_match(message, "pick")
  expect_match(message, paste0("\\b", rough$pick_idx, "\\b"))
  expect_match(message, paste0("protein[^0-9]*", bad_pick[1], "\\b"))
  through_fit <- .with_warnings(fit_prnode(small, max_iter = 1))
  expect_true(any(grepl("pick", through_fit$warnings)))

  base <- .fixture_tuning_small()
  expect_true(all(base$path$converged[, base$pick_idx]))
  expect_no_warning(tune_prnode(small, path = base$path))
})

# --- §3 the criterion along the path ---------------------------------------------

## C4: at the null model k = 0 and RSS = ||y_j||^2, so the criterion is the
## weighted total sum of squares (plan row 81).
test_that("C4: at the null model the criterion is sum_j ||y_j||^2", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  expect_equal(tn$criterion[1], sum(public$y^2), tolerance = 1e-10)
  expect_identical(tn$k[1], 0)
})

## C5: the criterion recomputed with nothing from the package but lm.fit(): per
## path point and protein the RSS of the least-squares refit on the support
## plus log(n) k; the pick is its first minimum, and `spar` is the design's.
test_that("C5: the criterion along the path equals the inline lm.fit() recomputation, and the pick is its minimum", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  x <- public$x
  y <- public$y
  n <- nrow(x)
  p <- ncol(y)
  L <- length(tn$lambda_path)
  crit <- numeric(L)
  k_vec <- numeric(L)
  for(j in seq_len(p)){
    support <- .path_coef_mat(tn$path, j) != 0
    for(l in seq_len(L)){
      S <- which(support[, l])
      rss <- if(length(S) == 0) sum(y[, j]^2) else sum(stats::lm.fit(x[, S, drop = FALSE], y[, j])$residuals^2)
      crit[l] <- crit[l] + rss + log(n)*length(S)
      k_vec[l] <- k_vec[l] + length(S)
    }
  }
  expect_equal(tn$criterion, crit, tolerance = 1e-10)
  expect_equal(tn$k, k_vec)
  expect_equal(tn$rss, crit - log(n)*k_vec, tolerance = 1e-10)
  expect_identical(tn$pick_idx, which.min(crit))
  expect_equal(tn$lambda, tn$lambda_path[tn$pick_idx])
  expect_equal(tn$spar, public$spar)
})

## C6: ties and the first minimum (plan row 71): a path point repeated three
## times gives identical supports, hence identical criterion values, and
## `which.min()` returns the first, the largest lambda.
test_that("C6: among tied minima the pick is the first (sparsest) index", {
  small <- .fixture_dataset_small()$public
  base <- .fixture_tuning_small()
  a <- base$pick_idx
  lp <- base$lambda_path
  tied <- c(lp[seq_len(a)], lp[a], lp[a], lp[seq(a + 1, length(lp))])
  expect_false(is.unsorted(rev(tied)))
  tn <- .with_warnings(tune_prnode(small, lambda_path = tied))$value
  expect_identical(tn$criterion[a], tn$criterion[a + 1])
  expect_identical(tn$criterion[a + 1], tn$criterion[a + 2])
  expect_identical(tn$pick_idx, a)
  expect_identical(tn$lambda_path[a + 2], lp[a])
})

## C7: the known-variance criterion is not invariant to a common rescaling of
## the loss weights: scaling x and y by sqrt(2) (omega doubled) doubles every
## RSS and leaves the supports unchanged, so the criterion becomes
## 2 RSS + log(n) k. The supports are unchanged because the penalty is on the
## fitted contribution ||x_g beta_g|| (Definition 1), of the derivation notes so the
## objective at (sqrt(2) x, sqrt(2) y, sqrt(2) lambda) is twice the objective
## at (x, y, lambda) and the grid, from lambda_max, scales by sqrt(2) with it.
## The profile-likelihood BIC was invariant; this test documents the
## difference, not an invariance, and why the design's weights matter.
test_that("C7: doubling omega doubles the criterion's fit term and leaves its penalty term", {
  small <- .fixture_dataset_small()$public
  base <- .fixture_tuning_small()
  scaled <- small
  scaled$x <- sqrt(2)*small$x
  scaled$y <- sqrt(2)*small$y
  scaled$omega <- 2*small$omega
  tn <- .with_warnings(tune_prnode(scaled, num_lambda = 20, num_decades = 2))$value
  expect_equal(tn$lambda_path, sqrt(2)*base$lambda_path, tolerance = 1e-12)
  for(j in seq_len(4)){
    expect_identical(.path_coef_mat(tn$path, j) != 0, .path_coef_mat(base$path, j) != 0, info = paste("protein", j))
  }
  n <- small$n_M + small$n_D
  expect_equal(tn$criterion, 2*base$rss + log(n)*base$k, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(tn$criterion, base$criterion)))
})

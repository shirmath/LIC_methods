# Plan, §2-§4 on the
# internal helpers of R/tune_prnode.R: the least-squares refit
# .refit_support() (R1-R4) and the known-variance BIC (C1-C3). The helpers are called bare: the tests run with the package namespace
# as their parent, under devtools::test() and R CMD check alike.

# --- §2 the refit ---------------------------------------------------------------

## R1: the coefficients are placed by `active_idx` position (plan row 44), so a
## gapped, non-consecutive support is the one that would expose a positional
## mix-up. No support here contains a regulator's S, M and D columns together,
## so lm.fit() has no aliased column and the comparison is exact.
test_that("R1: .refit_support() equals lm.fit() on each support, placed by position", {
  public <- .tuning_public()
  x <- public$x
  y <- public$y[, 1]
  num_col <- ncol(x)
  supports <- list(c(2, 7, 13), c(1, 2, 3), c(5, 9, 10, 18), 7)
  support_mat <- sapply(supports, function(s) seq_len(num_col) %in% s)
  beta <- .refit_support(x, y, support_mat)
  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(num_col, length(supports)))
  for(i in seq_along(supports)){
    label <- paste(supports[[i]], collapse = ",")
    expected <- numeric(num_col)
    expected[supports[[i]]] <- stats::lm.fit(x[, supports[[i]], drop = FALSE], y)$coefficients
    expect_equal(beta[, i], expected, tolerance = 1e-12, info = label)
    expect_true(all(beta[-supports[[i]], i] == 0), info = label)
  }
})

## R2: `x_S = x_M + x_D` column by column, so a support holding all three
## columns of one regulator is rank deficient by construction (plan row 42):
## lm.fit() pivots one column out (NA), the helper sets it to zero, and the
## fitted values equal those of the two-column support that spans the same
## space.
test_that("R2: an aliased scaffold column gets zero and the fitted values equal the two-column support's", {
  public <- .tuning_public()
  x <- public$x
  y <- public$y[, 2]
  num_col <- ncol(x)
  pq <- public$p*public$q
  k <- 1
  triple <- c(k, pq + k, 2*pq + k)
  support_mat <- cbind(seq_len(num_col) %in% triple, seq_len(num_col) %in% triple[2:3])
  beta <- .refit_support(x, y, support_mat)
  expect_false(anyNA(beta))
  expect_equal(sum(beta[triple, 1] == 0), 1)
  expect_false(all(beta[, 1] == 0))
  expect_equal(c(x %*% beta[, 1]), c(x %*% beta[, 2]), tolerance = 1e-10)
  expect_equal(x[, k], x[, pq + k] + x[, 2*pq + k])
})

## R3: the dedupe by `match()` on the pasted support keys (plan row 51): equal
## supports share one refit, the empty support is the zero column, and the
## output keeps the input's column order.
test_that("R3: the empty support is zero, duplicate supports are identical, columns stay in input order", {
  public <- .tuning_public()
  x <- public$x
  y <- public$y[, 3]
  num_col <- ncol(x)
  in_set <- function(s) seq_len(num_col) %in% s
  support_mat <- cbind(rep(FALSE, num_col), in_set(c(2, 7)), rep(FALSE, num_col), in_set(c(2, 7)), in_set(4))
  beta <- .refit_support(x, y, support_mat)
  expect_equal(ncol(beta), 5)
  expect_true(all(beta[, 1] == 0))
  expect_true(all(beta[, 3] == 0))
  expect_identical(beta[, 2], beta[, 4])
  expect_false(identical(beta[, 2], beta[, 5]))
  expect_equal(beta[4, 5], unname(stats::lm.fit(x[, 4, drop = FALSE], y)$coefficients))
  expect_true(all(beta[-4, 5] == 0))
})

## R4: `apply()` over one row and the final `drop = FALSE` (plan rows 43, 56):
## a one-column support and a one-row support both come back as matrices.
test_that("R4: a 3pq x 1 support returns 3pq x 1 and a 1 x L support returns 1 x L", {
  public <- .tuning_public()
  x <- public$x
  y <- public$y[, 1]
  num_col <- ncol(x)
  one <- .refit_support(x, y, matrix(seq_len(num_col) %in% c(2, 7), num_col, 1))
  expect_true(is.matrix(one))
  expect_equal(dim(one), c(num_col, 1))
  row <- .refit_support(x[, 2, drop = FALSE], y, matrix(c(TRUE, FALSE, TRUE), 1, 3))
  expect_true(is.matrix(row))
  expect_equal(dim(row), c(1, 3))
  expect_identical(row[1, 2], 0)
  expect_identical(row[1, 1], row[1, 3])
  expect_equal(row[1, 1], unname(stats::lm.fit(x[, 2, drop = FALSE], y)$coefficients))
})

# --- §3 the known-variance BIC -------------------------------------------------

## C1: the criterion is sum_j [RSS_jl + log(n) k_jl] on the least-squares refit
## of each support, with k the active-group count; computed here inline with
## lm.fit() on a small random design and two hand-made supports per protein
## (the empty one included, where the criterion is the null RSS).
test_that("C1: .known_variance_bic() is sum_j RSS_jl + log(n) k_jl on the refit, with k and rss returned summed over proteins", {
  set.seed(3)
  n <- 40
  x <- matrix(stats::rnorm(n*6), n, 6)
  y <- cbind(x[, 1] - 2*x[, 4] + stats::rnorm(n), stats::rnorm(n))
  support_1 <- cbind(rep(FALSE, 6), c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE), rep(TRUE, 6))
  support_2 <- cbind(rep(FALSE, 6), c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  out <- .known_variance_bic(x, y, list(support_1, support_2), q = 1)
  rss_of <- function(j, support){
    vapply(seq_len(ncol(support)), function(l){
      cols <- which(support[, l])
      if(length(cols) == 0) return(sum(y[, j]^2))
      sum(stats::lm.fit(x[, cols, drop = FALSE], y[, j])$residuals^2)
    }, numeric(1))
  }
  rss <- rss_of(1, support_1) + rss_of(2, support_2)
  k <- colSums(support_1) + colSums(support_2)
  expect_equal(out$rss, rss, tolerance = 1e-10)
  expect_equal(out$k, k)
  expect_equal(out$criterion, rss + log(n)*k, tolerance = 1e-10)
})

## C2: at the null model the criterion is the sum of the squared responses,
## and a zero response (RSS exactly 0) is floored so that nothing is NaN.
test_that("C2: at the empty support the criterion is sum_j ||y_j||^2, and a zero response stays finite", {
  public <- .tuning_public()
  n <- public$n_M + public$n_D
  empty <- matrix(FALSE, ncol(public$x), 1)
  out <- .known_variance_bic(public$x, public$y, rep(list(empty), ncol(public$y)), q = 1)
  expect_equal(out$criterion, sum(public$y^2), tolerance = 1e-10)
  expect_identical(out$k, 0)
  zero <- .known_variance_bic(public$x, matrix(0, n, 1), list(cbind(empty, !empty)), q = 1)
  expect_true(all(is.finite(zero$criterion)))
  expect_equal(zero$k, c(0, ncol(public$x)))
})

## C3: the criterion's penalty is log(n) with n the rows of the whole design
## (no effective-sample-size correction): doubling the rows of x and y adds
## log(2) k and doubles the fit term.
test_that("C3: the penalty is log(n) k with n the number of rows", {
  set.seed(8)
  n <- 30
  x <- matrix(stats::rnorm(n*4), n, 4)
  y <- matrix(x[, 1] + stats::rnorm(n), n, 1)
  support <- cbind(rep(FALSE, 4), c(TRUE, TRUE, FALSE, FALSE))
  one <- .known_variance_bic(x, y, list(support), q = 1)
  two <- .known_variance_bic(rbind(x, x), rbind(y, y), list(support), q = 1)
  expect_equal(two$rss, 2*one$rss, tolerance = 1e-10)
  expect_equal(two$criterion - 2*one$rss, log(2*n)*one$k, tolerance = 1e-10)
})

## E21: .expand_groups() widens a support to whole groups, so
## the refit is on "all columns of the active groups" as the roxygen says. The
## solver zeroes a group whole, so the partial group is hand-built.
test_that("E21: .expand_groups() turns one TRUE row of a q = 2 group into both, and the refit uses both columns", {
  support_mat <- matrix(FALSE, 6, 3)
  support_mat[2, 1] <- TRUE
  support_mat[c(3, 6), 2] <- TRUE
  expanded <- .expand_groups(support_mat, q = 2)
  expect_identical(expanded[, 1], c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(expanded[, 2], c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_identical(expanded[, 3], rep(FALSE, 6))
  expect_identical(.expand_groups(expanded, q = 2), expanded)

  public <- .fixture_public_q2()
  x <- public$x
  y <- public$y[, 1]
  partial <- matrix(FALSE, ncol(x), 1)
  partial[c(4, 9), 1] <- TRUE
  beta <- .refit_support(x, y, .expand_groups(partial, q = 2))
  cols <- c(3, 4, 9, 10)
  expected <- numeric(ncol(x))
  expected[cols] <- stats::lm.fit(x[, cols, drop = FALSE], y)$coefficients
  expect_equal(beta[, 1], expected, tolerance = 1e-12)
})

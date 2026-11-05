# Plan, §5-§7: the return
# contract of tune_prnode() and the two new branches of fit_prnode() (E2-E18;
# E1 is in test-tune_prnode_oracle.R, E8 beside T8), the print
# methods (P1-P3) and the invariances (I1-I3). Tests marked "Phase 4" assert
# the resolved decisions (E12, E16) and are red until the code is edited.

# --- §5 tune_prnode() end to end and the fit_prnode() branches ----------------------

## E2: the return contract, with the lengths asserted before the names (plan
## row 29) and `mu` renamed deliberately (row 12).
test_that("E2: a prnode_tuning carries the documented elements, in alphabetical order, with the documented shapes", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  L <- length(tn$lambda_path)
  expect_s3_class(tn, "prnode_tuning")
  expect_length(tn, 9)
  expect_identical(names(tn), c("criterion", "k", "lambda", "lambda_path", "mu", "path", "pick_idx", "rss", "spar"))
  expect_length(tn$criterion, L)
  expect_length(tn$k, L)
  expect_length(tn$rss, L)
  expect_length(tn$lambda, 1)
  expect_length(tn$lambda_path, L)
  expect_length(tn$mu, 3)
  expect_identical(names(tn$mu), c("S", "M", "D"))
  expect_s3_class(tn$path, "prnode_path")
  expect_length(tn$pick_idx, 1)
  expect_true(tn$pick_idx == round(tn$pick_idx) && tn$pick_idx >= 1 && tn$pick_idx <= L)
  expect_length(tn$spar, 1)
  expect_equal(tn$spar, public$spar)
  expect_true(all(tn$k >= 0 & tn$k <= 3*ncol(public$y)^2))

  small <- .fixture_dataset_small()$public
  unnamed <- .with_warnings(tune_prnode(small, mu = c(1, 1, 1), num_lambda = 2, num_decades = 1))$value
  expect_identical(unnamed$mu, c(S = 1, M = 1, D = 1))
})

## E3: a tuning passed as `lambda` fits at its lambda and mu and is carried as
## `tuning`; the single-lambda subset drops it.
test_that("E3: fit_prnode(public, lambda = tuning) equals the numeric-lambda fit plus the tuning; fit[1] drops it", {
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_small()
  tuned <- .with_warnings(fit_prnode(small, lambda = tn))$value
  plain <- .with_warnings(fit_prnode(small, lambda = tn$lambda))$value
  expect_identical(tuned$tuning, tn)
  expect_null(plain$tuning)
  expect_identical(tuned$lambda, tn$lambda)
  expect_identical(tuned$mu, tn$mu)
  tuned$tuning <- NULL
  expect_identical(tuned, plain)
  expect_null(.with_warnings(fit_prnode(small, lambda = tn))$value[1]$tuning)
})

## E4: the lambda = NULL route is tune_prnode() at the defaults, and the
## deterministic pair is identical (plan rows 65, 85).
test_that("E4: fit_prnode(public) is identical to fit_prnode(public, lambda = tune_prnode(public))", {
  small <- .fixture_dataset_small()$public
  default_fit <- .fixture_fit_default_small()$fit
  tn <- .fixture_tuning_default_small()
  expect_identical(default_fit$tuning, tn)
  expect_identical(default_fit, .with_warnings(fit_prnode(small, lambda = tn))$value)
  expect_length(tn$lambda_path, 80)
  expect_equal(tn$lambda_path[80], tn$lambda_path[1]*1e-4, tolerance = 1e-12)
})

## E5: the mu handshake and the entry guards of both functions (plan row 80).
test_that("E5: a disagreeing mu errors, an agreeing one passes, a bare list errors in both functions, a bad verbose errors", {
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_small()
  expect_error(fit_prnode(small, lambda = tn, mu = c(1, 0.9, 0.9)), regexp = "disagrees with the tuning's mu")
  named <- .with_warnings(fit_prnode(small, lambda = tn, mu = c(S = 1, M = 1, D = 1)))$value
  unnamed <- .with_warnings(fit_prnode(small, lambda = tn, mu = c(1, 1, 1)))$value
  expect_identical(named, unnamed)
  expect_identical(named$lambda, tn$lambda)
  expect_error(tune_prnode(list(a = 1)), regexp = "must be the design list")
  expect_error(fit_prnode(list(a = 1)), regexp = "must be the design list")
  expect_error(tune_prnode(small, verbose = -1), regexp = "`verbose`")
  expect_error(tune_prnode(small, verbose = "a"), regexp = "`verbose`")
})

## E6: the tuning sees the reweighted design (plan row 4): a user omega is
## applied before the penalty is chosen, so the tuning's path carries the new
## omega and its lambda is tune_prnode()'s on the manually reweighted public.
test_that("E6: fit_prnode(public, omega = w) tunes on the reweighted design", {
  small <- .fixture_dataset_small()$public
  w <- unname(c(2, 0.5)*small$omega)
  fit_w <- .with_warnings(fit_prnode(small, omega = w))$value
  expect_equal(fit_w$tuning$path$design$omega, c(M = w[1], D = w[2]))
  expect_equal(fit_w$omega, c(M = w[1], D = w[2]))
  m_rows <- seq_len(small$n_M)
  d_rows <- small$n_M + seq_len(small$n_D)
  reweighted <- small
  reweighted$x[m_rows, ] <- sqrt(w[1]/small$omega[[1]])*small$x[m_rows, ]
  reweighted$x[d_rows, ] <- sqrt(w[2]/small$omega[[2]])*small$x[d_rows, ]
  reweighted$y[m_rows, ] <- sqrt(w[1]/small$omega[[1]])*small$y[m_rows, ]
  reweighted$y[d_rows, ] <- sqrt(w[2]/small$omega[[2]])*small$y[d_rows, ]
  reweighted$omega <- c(M = w[1], D = w[2])
  tn_w <- .with_warnings(tune_prnode(reweighted))$value
  expect_equal(fit_w$tuning$lambda, tn_w$lambda, tolerance = 1e-12)
  expect_identical(fit_w$tuning$pick_idx, tn_w$pick_idx)
  expect_identical(fit_w$tuning$pick_idx, tn_w$pick_idx)
  expect_equal(fit_w$tuning$criterion, tn_w$criterion, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(tn_w$lambda_path[1], .fixture_tuning_default_small()$lambda_path[1])))
})

## E7: every verbose level runs and prints what it says (plan row 48): 0 is
## silent, 1 names the phases, 2 adds the solver's per-protein lines, and
## through fit_prnode(public, verbose = 1) the tuning is decremented to silent.
test_that("E7: verbose 0 is silent, 1 reports the phases, 2 adds the solver's lines, and fit_prnode() decrements it", {
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_small()
  expect_silent(tune_prnode(small, path = tn$path))
  phases <- utils::capture.output(tune_prnode(small, path = tn$path, verbose = 1))
  expect_true(any(grepl("known-variance BIC", phases)))
  expect_true(any(grepl("pick at path index", phases)))
  expect_false(any(grepl("fitting the path", phases)))
  fitted <- utils::capture.output(.with_warnings(tune_prnode(small, num_lambda = 2, num_decades = 1, verbose = 1)))
  expect_true(any(grepl("fitting the path of 2 lambda values for 4 proteins", fitted)))
  expect_false(any(grepl("fit_prnode_path", fitted)))
  solver <- utils::capture.output(.with_warnings(tune_prnode(small, num_lambda = 2, num_decades = 1, verbose = 2)))
  expect_true(any(grepl("fit_prnode_path: protein", solver)))
  wrapper <- .fixture_fit_default_small()$output
  expect_true(any(grepl("no `lambda` given", wrapper, fixed = TRUE)))
  expect_false(any(grepl("tune_prnode:", wrapper, fixed = TRUE)))
  expect_false(any(grepl("fit_prnode_path", wrapper, fixed = TRUE)))
})

## E9 (Phase 4, Z3 (b)): at q = 2, k counts active groups, not nonzero
## columns. The solver zeroes a group whole, so every active group has both
## columns nonzero and the group count is half the column count; the refit
## still uses every column of the active groups (the RSS of lm.fit() on them).
test_that("E9: at q = 2 the criterion's k is the number of active groups and the refit uses all their columns", {
  public <- .fixture_public_q2()
  tq2 <- .fixture_tuning_q2()
  x <- public$x
  y <- public$y
  n_M <- public$n_M
  n_D <- public$n_D
  n <- n_M + n_D
  p <- ncol(y)
  L <- length(tq2$lambda_path)
  expect_equal(public$q, 2)
  expect_equal(ncol(x), 3*p*2)
  crit_groups <- numeric(L)
  crit_columns <- numeric(L)
  for(j in seq_len(p)){
    support <- .path_coef_mat(tq2$path, j) != 0
    odd <- seq(1, nrow(support), by = 2)
    expect_identical(support[odd, ], support[odd + 1, ], info = paste("protein", j))
    for(l in seq_len(L)){
      S <- which(support[, l])
      rss <- if(length(S) == 0) sum(y[, j]^2) else sum(stats::lm.fit(x[, S, drop = FALSE], y[, j])$residuals^2)
      crit_groups[l] <- crit_groups[l] + rss + log(n)*.active_groups(tq2$path, j, l)
      crit_columns[l] <- crit_columns[l] + rss + log(n)*length(S)
    }
  }
  expect_false(isTRUE(all.equal(crit_groups, crit_columns)))
  expect_equal(tq2$criterion, crit_groups, tolerance = 1e-10)
  expect_identical(tq2$pick_idx, which.min(crit_groups))
})

## E10: the chosen lambda is a point of the path and the single-lambda fit at
## it (a cold solve) reproduces the warm-started path solution.
test_that("E10: the fit at the tuning's lambda reproduces the path's coefficients at pick_idx", {
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_small()
  fit <- .with_warnings(fit_prnode(small, lambda = tn))$value
  l <- tn$pick_idx
  expect_identical(fit$lambda, tn$lambda_path[l])
  expect_equal(fit$path$theta_D[, , 1], tn$path$theta_D[, , l], tolerance = 1e-8)
  expect_equal(fit$path$theta_M[, , 1], tn$path$theta_M[, , l], tolerance = 1e-8)
  expect_equal(fit$d_ode_est, .unscale_coef(matrix(tn$path$theta_D[, , l], 4, 4), small$response_scale, 1),
               tolerance = 1e-8)
  expect_true(any(fit$d_ode_est != 0))
})

## E11: a deliberately weak sanity check (plan row 67): the pick is strictly
## inside the path and the daughter's network is neither empty nor full. No
## MCC tolerance; recovery is the bank's business.
test_that("E11: on the low-noise fixture the pick is strictly inside the path and the daughter network is non-trivial", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  L <- length(tn$lambda_path)
  expect_gt(tn$pick_idx, 1)
  expect_lt(tn$pick_idx, L)
  fit <- .with_warnings(fit_prnode(public, lambda = tn))$value
  edges <- sum(.as_adj(fit$d_network_est))
  expect_gte(edges, 1)
  expect_lt(edges, public$p^2)
})

## E12 (Phase 4): a non-finite entry in `public$y` or `public$x` errors at the
## top of tune_prnode(), naming the protein (y) or the regulator and block (x,
## mapped through the S, M, D column order of .path_coef_mat()), before any
## path is fitted and through fit_prnode(public).
test_that("E12: a non-finite entry in public$y names the protein; in public$x the regulator and block; before the path and through fit_prnode()", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  for(bad in list(NA, NaN, Inf)){
    label <- paste("y", bad)
    broken <- public
    broken$y[5, 3] <- bad
    message <- .error_message(tune_prnode(broken, path = tn$path))
    expect_match(message, "public\\$y", info = label)
    expect_match(message, "protein_3", info = label)
    expect_no_match(message, "protein_[1245]", info = label)
  }
  two <- public
  two$y[5, 3] <- NA
  two$y[7, 5] <- Inf
  message <- .error_message(tune_prnode(two, num_lambda = 2, num_decades = 1))
  expect_match(message, "protein_3")
  expect_match(message, "protein_5")
  expect_match(.error_message(fit_prnode(two)), "protein_3")
  unnamed <- two
  unnamed$protein_names <- NULL
  expect_match(.error_message(tune_prnode(unnamed, path = tn$path)), "\\b3\\b")

  pq <- public$p*public$q
  for(bad in list(NA, NaN, Inf)){
    label <- paste("x", bad)
    broken <- public
    broken$x[5, pq + 2] <- bad
    message <- .error_message(tune_prnode(broken, path = tn$path))
    expect_match(message, "public\\$x", info = label)
    expect_match(message, "protein_2", info = label)
    expect_match(message, "\\bM\\b", info = label)
  }
  two_x <- public
  two_x$x[5, pq + 2] <- NA
  two_x$x[6, 2*pq + 3] <- Inf
  message <- .error_message(tune_prnode(two_x, num_lambda = 2, num_decades = 1))
  expect_match(message, "protein_2")
  expect_match(message, "protein_3")
  expect_match(message, "\\bD\\b")
  expect_match(.error_message(fit_prnode(two_x)), "protein_2")
})

## E13: p = 1 end to end (plan row 7, degenerate sizes): one response, three
## design columns.
test_that("E13: a one-protein design tunes, with one criterion value per path point", {
  ds <- .fixture_dataset_small()
  one <- .one_protein_public(ds, 2)
  expect_equal(ncol(one$x), 3)
  expect_equal(ncol(one$y), 1)
  tn <- .with_warnings(tune_prnode(one, num_lambda = 20, num_decades = 2))$value
  expect_length(tn$criterion, 20)
  expect_lte(max(tn$k), 3)
  fit <- .with_warnings(fit_prnode(one, lambda = tn))$value
  expect_equal(dim(fit$d_ode_est), c(1, 1))
  expect_identical(fit$tuning, tn)
})

## E14: no RNG anywhere in the rule (plan row 65).
test_that("E14: two tunings of the same public are identical", {
  small <- .fixture_dataset_small()$public
  again <- .with_warnings(tune_prnode(small, num_lambda = 20, num_decades = 2))$value
  expect_identical(again, .fixture_tuning_small())
})

## E15: the blast radius of the new element (plan row 91): the downstream
## readers of a prnode_fit ignore `tuning` (inheritance_flow() reads
## `tuning` by design, test-inheritance_flow.R). plot_md_adjacency() takes the
## network matrix, so the daughter's network is what is passed.
test_that("E15: score_md_networks(), predict() and plot_md_adjacency() give the same values with and without `tuning`", {
  ds <- .fixture_dataset_small()
  public <- ds$public
  tn <- .fixture_tuning_small()
  tuned <- .with_warnings(fit_prnode(public, lambda = tn))$value
  plain <- .with_warnings(fit_prnode(public, lambda = tn$lambda))$value
  expect_identical(score_md_networks(tuned, ds$truth), score_md_networks(plain, ds$truth))
  expect_identical(stats::predict(tuned, cell = "D"), stats::predict(plain, cell = "D"))
  expect_identical(stats::predict(tuned, cell = "M"), stats::predict(plain, cell = "M"))
  drawn <- .on_null_device(list(plot_md_adjacency(tuned$d_network_est, protein_names = public$protein_names),
                                plot_md_adjacency(plain$d_network_est, protein_names = public$protein_names)))
  expect_identical(drawn[[1]], drawn[[2]])
})

## E16 (Phase 4): a sparse design is never expected: any
## Matrix-class `public$x` errors at entry with a named message, in both
## functions; the base matrix passes.
test_that("E16: a Matrix-class public$x errors naming sparse formats in tune_prnode() and fit_prnode()", {
  public <- .tuning_public()
  tn <- .fixture_tuning()
  sparse <- public
  sparse$x <- Matrix::Matrix(public$x, sparse = TRUE)
  expect_s4_class(sparse$x, "dgCMatrix")
  expect_match(.error_message(tune_prnode(sparse, path = tn$path)), "sparse")
  expect_match(.error_message(tune_prnode(sparse, num_lambda = 2, num_decades = 1)), "sparse")
  expect_match(.error_message(fit_prnode(sparse)), "sparse")
  dense <- public
  dense$x <- Matrix::Matrix(public$x, sparse = FALSE)
  expect_match(.error_message(tune_prnode(dense, path = tn$path)), "sparse")
  expect_s3_class(tune_prnode(public, path = tn$path), "prnode_tuning")
})

## E19: `y` may hold a subset of the proteins (as T7 builds),
## so `q` is read from the design's `q` or `index`, never from
## ncol(x) / (3 ncol(y)), which is 2 on the 2-of-4 subset below. The tuning
## and print() run, the anchor's printed active groups are its nonzero
## columns (q = 1), the messages name the regulator through protein_names and
## the response by index (protein_names has four names for two columns), and
## a design whose `index` and `x` disagree is refused.
test_that("E19: a subset-y design reads q from the design, tunes, prints, and names its columns", {
  small <- .fixture_dataset_small()$public
  subset <- small
  subset$y <- small$y[, c(1, 3), drop = FALSE]
  colnames(subset$y) <- NULL
  subset$response_scale <- small$response_scale[c(1, 3)]
  expect_equal(ncol(subset$x)/(3*ncol(subset$y)), 2)
  expect_identical(.design_q(subset, pq = ncol(subset$x)/3, p_y = ncol(subset$y)), 1L)
  tn <- .with_warnings(tune_prnode(subset, num_lambda = 20, num_decades = 2))$value
  expect_s3_class(tn, "prnode_tuning")
  printed <- utils::capture.output(print(tn))
  num_nonzero <- sum(vapply(seq_len(2), function(j) sum(.path_coef_mat(tn$path, j)[, tn$pick_idx] != 0),
                            numeric(1)))
  pick_line <- grep("^  pick", printed, value = TRUE)
  expect_match(pick_line, paste0("active groups = ", num_nonzero, "\\b"))

  pq <- ncol(subset$x)/3
  bad_x <- subset
  bad_x$x[5, pq + 2] <- Inf
  expect_match(.error_message(tune_prnode(bad_x, path = tn$path)), "protein_2 \\(M block\\)")
  bad_y <- subset
  bad_y$y[5, 2] <- Inf
  message_y <- .error_message(tune_prnode(bad_y, path = tn$path))
  expect_match(message_y, "protein\\(s\\) 2$")

  mismatch <- subset
  mismatch$q <- NULL
  mismatch$index <- rep(seq_len(9), length.out = ncol(subset$x))
  expect_match(.error_message(tune_prnode(mismatch, num_lambda = 2, num_decades = 1)), "cannot read the group size")
})

## E20: a non-numeric `y` is its own error, not reported as
## non-finite in every protein.
test_that("E20: a character public$y errors 'must be numeric'", {
  public <- .tuning_public()
  character_y <- public
  character_y$y <- matrix(as.character(public$y), nrow(public$y))
  message <- .error_message(tune_prnode(character_y, num_lambda = 2, num_decades = 1))
  expect_match(message, "public\\$y` must be numeric")
  expect_no_match(message, "non-finite")
})

## E17: perturbation (plan regimes row): noise of 1e-9 times the response's
## scale leaves the anchor and the pick where they were; the criterion's
## margin at the anchor is measured first so the fixture is known to be off
## a tie.
test_that("E17: a 1e-9 relative perturbation of y leaves pick_idx unchanged", {
  small <- .fixture_dataset_small()$public
  base <- .fixture_tuning_small()
  expect_gt(sort(base$criterion)[2] - min(base$criterion), 1e-6)
  perturbed <- small
  set.seed(17)
  perturbed$y <- small$y + 1e-9*max(abs(small$y))*matrix(stats::rnorm(length(small$y)), nrow(small$y))
  tn <- .with_warnings(tune_prnode(perturbed, num_lambda = 20, num_decades = 2))$value
  expect_identical(tn$pick_idx, base$pick_idx)
  expect_equal(tn$criterion, base$criterion, tolerance = 1e-6)
})

## E18: serialization (plan regimes row): a prnode_tuning round-trips through
## saveRDS()/readRDS() and still drives fit_prnode().
test_that("E18: a prnode_tuning round-trips through RDS and fits identically", {
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_small()
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(tn, path)
  restored <- readRDS(path)
  expect_identical(restored, tn)
  expect_identical(.with_warnings(fit_prnode(small, lambda = restored))$value,
                   .with_warnings(fit_prnode(small, lambda = tn))$value)
})

# --- §6 printing --------------------------------------------------------------------

## P1: print.prnode_tuning().
test_that("P1: print.prnode_tuning() reports the path, the pick, its active groups and the smoothing level, invisibly", {
  tn <- .fixture_tuning_small()
  L <- length(tn$lambda_path)
  lines <- utils::capture.output(res <- withVisible(print(tn)))
  expect_false(res$visible)
  expect_identical(res$value, tn)
  out <- paste(lines, collapse = "\n")
  expect_match(out, "PRNODE tuning (class prnode_tuning)", fixed = TRUE)
  expect_match(out, paste0("path of ", L, " lambda values from ", signif(tn$lambda_path[1], 3), " to ",
                           signif(tn$lambda_path[L], 3)), fixed = TRUE)
  expect_match(out, paste0("pick (known-variance BIC): index ", tn$pick_idx), fixed = TRUE)
  expect_match(out, paste0("active groups = ", tn$k[tn$pick_idx]), fixed = TRUE)
  expect_match(out, paste0("smoothing level spar = ", tn$spar), fixed = TRUE)
})

## P2: print.prnode_fit() adds its tuning line only when `tuning` exists.
test_that("P2: print.prnode_fit() reports the pick only for a tuned fit", {
  small <- .fixture_dataset_small()$public
  tn <- .fixture_tuning_small()
  tuned <- .with_warnings(fit_prnode(small, lambda = tn))$value
  plain <- .with_warnings(fit_prnode(small, lambda = tn$lambda))$value
  out_tuned <- paste(utils::capture.output(print(tuned)), collapse = "\n")
  expect_match(out_tuned, paste0("lambda chosen by tune_prnode(): pick at path index ", tn$pick_idx, " of ",
                                 length(tn$lambda_path)), fixed = TRUE)
  out_plain <- paste(utils::capture.output(print(plain)), collapse = "\n")
  expect_no_match(out_plain, "tune_prnode", fixed = TRUE)
  expect_match(out_plain, paste0("lambda = ", signif(tn$lambda, 3)), fixed = TRUE)
})

## P3: the roxygen examples run as documented (plan regimes row R15), read from
## man/ so the test follows the documentation rather than a copy of it; the
## tune_prnode() example completes in under 5 s (1 s at Phase 3). Under
## R CMD check man/ is not reachable from tests/ and the check runs the
## examples itself, so the test skips there.
test_that("P3: the examples of tune_prnode() and fit_prnode() run, the former in under 5 s", {
  man_dir <- test_path("..", "..", "man")
  skip_if_not(dir.exists(man_dir), "man/ is not reachable from the test directory (R CMD check runs the examples)")
  for(topic in c("tune_prnode", "fit_prnode")){
    rd_file <- file.path(man_dir, paste0(topic, ".Rd"))
    skip_if_not(file.exists(rd_file), paste(rd_file, "is missing"))
    example_file <- withr::local_tempfile(fileext = ".R")
    tools::Rd2ex(rd_file, example_file)
    example_env <- new.env(parent = environment())
    seconds <- system.time({
      utils::capture.output(res <- .with_warnings(source(example_file, local = example_env, echo = FALSE)))
    })[["elapsed"]]
    expect_true(exists("tuning", envir = example_env, inherits = FALSE), info = topic)
    expect_s3_class(get("tuning", envir = example_env), "prnode_tuning")
    if(topic == "tune_prnode") expect_lt(seconds, 5)
  }
})

# --- §7 invariances -------------------------------------------------------------------

## I1: protein relabelling (plan row 6, 83): the design is rebuilt from the
## permuted raw (the I-series' construction), so the response columns and the
## regulator groups move together. rho and n_eff are permuted; the criterion
## and the cv error are sums over proteins and the anchor, pick and lambda are
## unchanged. The sums run in a different order, so the comparison is to
## 1e-8 rather than identical.
test_that("I1: permuting the proteins leaves the criterion, the pick and lambda unchanged", {
  ds <- .fixture_dataset_small()
  raw <- ds$raw
  base <- .fixture_tuning_small()
  expect_identical(.rebuild_public(raw)$x, ds$public$x)
  for(perm in list(c(3, 1, 4, 2), c(2, 4, 1, 3))){
    label <- paste(perm, collapse = "")
    raw_p <- raw
    raw_p$Y_M <- raw$Y_M[, perm]
    raw_p$Y_D <- raw$Y_D[, perm]
    moved <- .with_warnings(tune_prnode(.rebuild_public(raw_p), num_lambda = 20, num_decades = 2))$value
    expect_equal(moved$lambda_path, base$lambda_path, tolerance = 1e-12, info = label)
    expect_equal(moved$criterion, base$criterion, tolerance = 1e-8, info = label)
    expect_identical(moved$pick_idx, base$pick_idx, info = label)
    expect_equal(moved$lambda, base$lambda, tolerance = 1e-12, info = label)
  }
})

## I2: row reversal within each cell leaves every residual sum of squares and
## the supports unchanged (the loss is a sum over rows and the design has no
## time order of its own), so the criterion and the pick are unchanged.
test_that("I2: reversing the rows within each cell leaves the criterion and the pick unchanged", {
  small <- .fixture_dataset_small()$public
  base <- .fixture_tuning_small()
  rows <- c(rev(seq_len(small$n_M)), small$n_M + rev(seq_len(small$n_D)))
  flipped <- small
  flipped$x <- small$x[rows, , drop = FALSE]
  flipped$y <- small$y[rows, , drop = FALSE]
  tn <- .with_warnings(tune_prnode(flipped, num_lambda = 20, num_decades = 2))$value
  expect_equal(tn$criterion, base$criterion, tolerance = 1e-8)
  expect_identical(tn$pick_idx, base$pick_idx)
})

## I3: names are never read (plan regimes row R6): `mu` named or not, and a
## public without protein names, tune identically (the path carries the
## design's settings, so it is compared after dropping them).
test_that("I3: a named and an unnamed mu, and a public without protein names, tune identically", {
  small <- .fixture_dataset_small()$public
  named <- .with_warnings(tune_prnode(small, mu = c(S = 1, M = 1, D = 1), num_lambda = 5, num_decades = 1))$value
  unnamed <- .with_warnings(tune_prnode(small, mu = c(1, 1, 1), num_lambda = 5, num_decades = 1))$value
  expect_identical(named, unnamed)
  nameless <- small
  nameless$protein_names <- NULL
  tn <- .with_warnings(tune_prnode(nameless, num_lambda = 5, num_decades = 1))$value
  expect_null(tn$path$design$protein_names)
  tn$path <- NULL
  named$path <- NULL
  expect_identical(tn, named)
})

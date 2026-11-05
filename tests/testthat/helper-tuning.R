# Fixtures and the glmnet oracle for the tuning-rule suite.
# The estimator's fixtures stay in helper-fixtures.R; this file adds
# the cached tunings, the examples' p = 4 dataset, two hand-built designs
# that reach guards make_prnode_design() never lets through, and a copy of
# the scored rule (the 0.5.0 rule's glmnet path solver) under `oracle_`
# names for E1.
#
# Every tuning is built once per run and cached in `.fx` together with the
# warnings it raised, so that the solver's max_iter warning (plan Z1) is
# captured at the fixture and never escapes a test (pitfalls row 27).

# --- warnings and errors as data ----------------------------------------------

# run `expr`; return its value with the messages of the warnings it raised
.with_warnings <- function(expr){
  messages <- character(0)
  value <- withCallingHandlers(expr, warning = function(cnd){
    messages <<- c(messages, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = messages)
}

# the message of the error `expr` raises, or a sentinel string when it raises
# none (so an expect_match() on the result fails with a readable message)
.error_message <- function(expr){
  tryCatch({
    force(expr)
    "<no error raised>"
  }, error = function(cnd) conditionMessage(cnd))
}

# --- datasets -----------------------------------------------------------------

# the examples' dataset of ?tune_prnode and ?fit_prnode: p = 4, 36 mother
# frames, 18 daughter frames, low noise, seed 1, 100 steps per frame
.fixture_dataset_small <- function(){
  .skip_unless_new_pipeline()
  key <- "dataset_small_p4"
  if(!exists(key, envir = .fx)){
    assign(key,
           build_md_dataset(seed = 1,
                            num_d = 18,
                            noise_level_ratio = 0.1,
                            truth_args = list(p = 4),
                            obs_time = seq(0, 1, length.out = 36),
                            steps_per_frame = 100),
           envir = .fx)
  }
  get(key, envir = .fx)
}

# the p = 6 fixture at the ODE level: 48 mother frames, `num_d` daughter frames
.tuning_public <- function(num_d = 12){
  .fixture_dataset("ode", num_d = num_d)$public
}

# the design rebuilt from a dataset's `raw` with the arguments
# build_md_dataset() used (reproduces `public` exactly; asserted in I1)
.rebuild_public <- function(raw,
                            q = 1,
                            x_list = NULL){
  make_prnode_design(m_obs_time = raw$m_obs_time,
                     Y_M = raw$Y_M,
                     d_obs_time = raw$d_obs_time,
                     Y_D = raw$Y_D,
                     times_e = raw$times_e,
                     q = q,
                     x_list = x_list)
}

# a one-protein design (p = 1, so pq = q): the raw of `ds` restricted to
# protein `j`, through make_prnode_design()
.one_protein_public <- function(ds,
                                j = 1){
  raw <- ds$raw
  make_prnode_design(m_obs_time = raw$m_obs_time,
                     Y_M = raw$Y_M[, j, drop = FALSE],
                     d_obs_time = raw$d_obs_time,
                     Y_D = raw$Y_D[, j, drop = FALSE],
                     times_e = raw$times_e,
                     q = ds$public$q)
}

# a `public` with the daughter cut to her first `n_D` rows of `x` and `y`
# (the profile is left as built). Not a valid profiled design, but a numeric
# design the tuning and the solver accept; it is the only way to reach the
# daughter-length guard of plan V3, since make_prnode_design() refuses fewer
# than 3 frames per cell and its smoother fewer than 4.
.truncated_public <- function(public,
                              n_D){
  stopifnot(n_D <= public$n_D)
  keep <- c(seq_len(public$n_M), public$n_M + seq_len(n_D))
  out <- public
  out$x <- public$x[keep, , drop = FALSE]
  out$y <- public$y[keep, , drop = FALSE]
  out$n_D <- n_D
  out$d_obs_time <- public$d_obs_time[seq_len(n_D)]
  out
}

# --- cached tunings -----------------------------------------------------------

# one tune_prnode() per key, with its warnings; `...` are tune_prnode()'s
.cached_tuning <- function(key,
                           public,
                           ...){
  .skip_unless_built("tune_prnode")
  key <- paste0("tuning_", key)
  if(!exists(key, envir = .fx)){
    assign(key, .with_warnings(tune_prnode(public, ...)), envir = .fx)
  }
  get(key, envir = .fx)
}

# the p = 6 fixture tuning: a 40-point, 3-decade path (1.7 s at Phase 3,
# no non-converged fit; the 30-point path of the same span hit max_iter once)
.fixture_tuning <- function(num_d = 12,
                            num_lambda = 40,
                            num_decades = 3){
  .cached_tuning(paste("p6", num_d, num_lambda, num_decades, sep = "_"), .tuning_public(num_d),
                 num_lambda = num_lambda, num_decades = num_decades)$value
}

.fixture_tuning_warnings <- function(num_d = 12,
                                     num_lambda = 40,
                                     num_decades = 3){
  .cached_tuning(paste("p6", num_d, num_lambda, num_decades, sep = "_"), .tuning_public(num_d),
                 num_lambda = num_lambda, num_decades = num_decades)$warnings
}

# the examples' tuning: the p = 4 dataset with the roxygen example's 20-point,
# 2-decade path (1 s; converges everywhere since the strong rule's re-admission
# gate was fixed, estimator plan R7b; before that one (protein, lambda) warned)
.fixture_tuning_small <- function(){
  .cached_tuning("small_20_2", .fixture_dataset_small()$public, num_lambda = 20, num_decades = 2)$value
}

.fixture_tuning_small_warnings <- function(){
  .cached_tuning("small_20_2", .fixture_dataset_small()$public, num_lambda = 20, num_decades = 2)$warnings
}

# the examples' dataset at the defaults (80 points, 4 decades; 0.5 s)
.fixture_tuning_default_small <- function(){
  .cached_tuning("small_default", .fixture_dataset_small()$public)$value
}

# fit_prnode(public) at the defaults on the examples' dataset, with its
# verbose = 1 output and warnings kept beside it (E4, E7, T8)
.fixture_fit_default_small <- function(){
  .skip_unless_built("fit_prnode", "public")
  key <- "fit_small_default"
  if(!exists(key, envir = .fx)){
    public <- .fixture_dataset_small()$public
    output <- utils::capture.output(res <- .with_warnings(fit_prnode(public, verbose = 1)))
    assign(key, list(fit = res$value, warnings = res$warnings, output = output), envir = .fx)
  }
  get(key, envir = .fx)
}

# a q = 2 design on the examples' raw with the 20-point, 2-decade path (E9, Z3)
.fixture_public_q2 <- function(){
  key <- "public_small_q2"
  if(!exists(key, envir = .fx)){
    assign(key, .rebuild_public(.fixture_dataset_small()$raw, q = 2), envir = .fx)
  }
  get(key, envir = .fx)
}

.fixture_tuning_q2 <- function(){
  .cached_tuning("small_q2_20_2", .fixture_public_q2(), num_lambda = 20, num_decades = 2)$value
}

# --- support arithmetic shared by the tests ----------------------------------

# the number of active groups of protein `j` at path point `l`: a group is
# the `q` consecutive columns of one regulator in one block, in the column
# order of .path_coef_mat() (plan Z3 (b), T5(c))
.active_groups <- function(path,
                           j,
                           l){
  coef_vec <- .path_coef_mat(path, j)[, l]
  q <- dim(path$s)[2]/dim(path$s)[1]
  sum(colSums(matrix(coef_vec != 0, nrow = q)) > 0)
}

# --- the glmnet oracle: the path solver of the 0.5.0 tuning rule ---
#
# A glmnet path at the package's group penalty (which is exact at q = 1) and a
# least-squares refit, under `oracle_` names, so that E1 can read the
# known-variance BIC on a path of another solver against the package's
# tune_prnode(). The edits relative to the original are the `oracle_` prefix, a `lambdas` argument to oracle_mu_path(),
# and the qualified `glmnet::glmnet()` and `stats::lm.fit()` calls (E1 skips
# without glmnet). The scored rule's anchor and prune (oracle_neff_ic(),
# oracle_block_cv(), oracle_prune_pick()) left with the rule in 0.9.0.

oracle_path_fit <- function(x, y, index, M_groups, D_groups, lambdas, mu){
  stopifnot(all(table(index) == 1))
  n <- nrow(x)
  w <- ifelse(index %in% M_groups, mu[["M"]], ifelse(index %in% D_groups, mu[["D"]], mu[["S"]]))
  pf <- w*sqrt(colSums(x^2))
  lambda_g <- lambdas*sum(pf)/(n*ncol(x))
  tol_args <- if(utils::packageVersion("glmnet") >= "5.0"){
    list(control = list(thresh = 1e-10, maxit = 1e6))
  } else {
    list(thresh = 1e-10, maxit = 1e6)
  }
  fit <- suppressWarnings(do.call(glmnet::glmnet,
                                  c(list(x = x, y = y, lambda = lambda_g, penalty.factor = pf,
                                         standardize = FALSE, intercept = FALSE), tol_args)))
  coefs <- as.matrix(stats::coef(fit))[-1, , drop = FALSE]
  if(ncol(coefs) < length(lambdas)){
    coefs <- coefs[, c(seq_len(ncol(coefs)), rep(ncol(coefs), length(lambdas) - ncol(coefs))), drop = FALSE]
  }
  coefs
}

oracle_refit_coefs <- function(x_train, y_train, coefs){
  key <- apply(coefs != 0, 2, paste, collapse = "")
  first <- !duplicated(key)
  beta <- sapply(which(first), function(i){
    S <- which(coefs[, i] != 0)
    b <- rep(0, nrow(coefs))
    if(length(S) > 0){
      fit <- stats::lm.fit(x_train[, S, drop = FALSE], y_train)$coefficients
      fit[is.na(fit)] <- 0
      b[S] <- fit
    }
    b
  })
  beta <- matrix(beta, nrow = nrow(coefs))
  beta[, match(key, key[first]), drop = FALSE]
}

oracle_mu_path <- function(d, mu, lambdas){
  path <- lapply(1:d$p, function(j){
    oracle_path_fit(d$x, d$y[, j], d$index, d$M_groups, d$D_groups, lambdas, mu = mu)
  })
  for(j in 1:d$p) path[[j]][, 1] <- 0
  list(lambdas = lambdas, path = path)
}

# the known-variance BIC of a glmnet-style support list, in plain lm.fit()
oracle_kv_bic <- function(d, path){
  n <- d$n_M + d$n_D
  Reduce(`+`, lapply(1:d$p, function(j){
    rss <- colSums((d$y[, j] - d$x %*% oracle_refit_coefs(d$x, d$y[, j], path[[j]]))^2)
    rss + log(n)*colSums(path[[j]] != 0)
  }))
}

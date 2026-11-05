# Fixtures and accessors for the PRNODE test suite.
#
# Every fixture is built in code from a fixed seed and cached in `.fx` for the
# run; nothing is read from disk. The builders skip (never fail) when the
# function they need has not been built yet, so that the suite is green with
# skips until Phase 4 lands and the skip count says what is missing.
#
# ---------------------------------------------------------------------------
# ASSUMED SHAPES OF THE UNBUILT FUNCTIONS (Phase 4 builds to these, or edits
# the accessor below once). Everything the plan leaves unspecified is decided
# here and nowhere else in the suite:
#
#   fit_prnode_path(design, lambda, ...) returns class "prnode_path" whose
#   coefficient elements always carry the lambda dimension (it is a path
#   object; `fit_prnode()` drops that dimension for a scalar lambda):
#     $s, $d_M, $d_D, $theta_M, $theta_D   p x pq x L arrays
#     $gamma_M, $gamma_D                   2 x p x L arrays (intercept, slope)
#     $networks$S, $M, $D                  p x p x L logical arrays
#     $objective, $gap, $kkt_max, $iterations, $converged   p x L matrices
#     $discount$M, $discount$D             length-p vectors (per regulator)
#     $lambda, $mu, $tol, $edge_tol        the settings
#
#   make_prnode_design(...) returns the `public` list of plan section 5 with
#   exactly these names: D_groups, d_obs_time, index, M_groups, m_obs_time,
#   n_D, n_M, omega, p, profile, protein_names, q, response_scale, times_e,
#   x, y. `omega` is a length-2 numeric in the order (M, D); `profile$M` and
#   `profile$D` are `profile_cell()` outputs (Z, qr_Z, y_tilde, Psi_tilde, n).
#
#   fit_prnode(public = NULL, lambda, ...) returns class "prnode_fit"; when
#   `public` is NULL the design is built from the named arguments
#   `m_obs_time, Y_M, d_obs_time, Y_D, times_e, q`. Its elements (F1, F7, F8,
#   F9, F11, F14, P2, I2, I3):
#     $s_ode_est, $m_ode_est, $d_ode_est     p x pq (scalar lambda) or
#                                            p x pq x L, dimnames[[3]] = lambda
#     $s_network_est, $m_network_est,
#     $d_network_est                         p x p logical, or p x p x L
#     $intercepts                            2 x p with rownames "mother",
#                                            "daughter" (row 2 of gamma_c, in
#                                            original units), or 2 x p x L
#     $m_init_conditions, $d_init_conditions length-p (row 1 of gamma_c, in
#                                            original units), or p x L
#     $path                                  the prnode_path (scaled units)
#     $lambda                                the lambda vector; fit[l]$lambda
#                                            is the scalar lambda[l]
#     $mu, $omega (length 2, order M, D), $response_scale (length p),
#     $converged (p x L logical), $discount   carried from the path/design
#     no element named joint_* (F14)
#   `[.prnode_fit` (fit[l]) returns the scalar shapes above; `print()` names
#   s_ (never joint_) and says whether the response was scaled; `predict(
#   object, cell)` is Z_c gamma_c + Psi_c theta_c at the frames (P2).
#
#   inheritance_flow(fit, bool_path = FALSE) (package 0.7.0) returns one
#   number from a fit carrying its tuning (test-inheritance_flow.R);
#   a dataset's truth$inheritance is the same measure on the true systems.
#
#   sim_md_trajectory(truth, ...) returns list(mother = list(X, obs_time,
#   Y_noiseless, Y, sigma), daughter = list(X, obs_time, Y_noiseless, Y,
#   sigma, birth_time, birth_idx), truth), where `birth_idx` is the index of
#   the birth frame in `obs_time` (M1, L7), not in `times_e`.
#
#   draw_md_truth(num_d, obs_time, truth_args, times_e, seed_number,
#   max_attempts, max_growth, ...) returns list(sim, num_attempts): `sim` is
#   the sim_md_trajectory() object above (sigma = 0, so `Y_noiseless` is what
#   B7 reads) and `num_attempts` the accepted attempt's index (B6, B7).
#
#   build_md_dataset() returns list(meta, raw, truth, public) with exactly the
#   names B5 lists; `raw$birth_time == raw$d_obs_time[1]`; `meta$level` is the
#   level label (Q3 exempts it from the identity between levels).
#
#   plot_md_pca3d() returns `scores`, `projection`, `colours` as lists of two
#   (mother, daughter) in `pair` order.
#
#   score_md_networks(fit, truth, ...) reads `fit$m_network_est`,
#   `fit$d_network_est`, `fit$s_network_est` (so a bare list with those three
#   elements is accepted) and returns `$mcc`, `$f1` (named M, D, S),
#   `$leak_rate`, `$borrow_rate` (named M, D) and `$confusion` (per cell, a
#   matrix with one row per edge type and columns TP, FP, FN, TN).
#
#   estimate_noise_sd(obs_time, Y, method) errors below four frames naming the
#   minimum (H8); make_prnode_design() with omega = NULL therefore refuses a
#   three-frame cell naming the cell and the minimum of 4 (D6), and accepts it
#   with a user omega (D1).
#
#   fit_prnode_path()$discount$M/$D are clamped at 1: min(1, 1 - min{W_c, W_S -
#   W_c'}/W_c) (R14, review pass 5).
# ---------------------------------------------------------------------------

.fx <- new.env(parent = emptyenv())

# The 48 frames sit at multiples of 1/47, which a grid of step 1e-5 does not
# contain, and the new pipeline reads frames off the grid by exact index (plan
# C2, M7). The fixture grids therefore have 2000 (fine) or 200 (coarse) steps
# per frame interval: step 1.06e-5 and 1.06e-4. The smoother and grid tests
# (B3, B4) keep the old banks' seq(0, 1, by = 1e-5) with nearest-neighbour
# lookup, because those two functions still serve that grid.
.fx_obs_time <- seq(0, 1, length.out = 48)
.fx_times_e <- seq(0, 1, by = 1/(47*2000))
.fx_times_coarse <- seq(0, 1, by = 1/(47*200))
.fx_seed <- 101
.fx_truth_args <- list(p = 6,
                       q = 1,
                       num_shared = 6,
                       num_mother_only = 2,
                       num_daughter_only = 2)
.fx_noise_low <- 0.1

# Is a function built, and (optionally) does it carry the argument that marks
# the rewritten signature? `fit_prnode()`, `build_md_dataset()`,
# `sim_md_trajectory()` and `add_md_noise()` exist today with the old
# signatures, so their tests look for an argument the rewrite adds.
.built <- function(fn_name,
                   marker = NULL){
  if(!exists(fn_name, mode = "function")) return(FALSE)
  if(is.null(marker)) return(TRUE)
  marker %in% names(formals(get(fn_name, mode = "function")))
}

.skip_unless_built <- function(fn_name,
                               marker = NULL){
  if(!.built(fn_name, marker)){
    skip(paste0(fn_name, "() is not built yet (test plan, Phase 4)",
                if(!is.null(marker)) paste0("; marker argument `", marker, "`")))
  }
}

.skip_unless_new_pipeline <- function(){
  .skip_unless_built("sim_md_network")
  .skip_unless_built("sim_ode_trajectory")
  .skip_unless_built("sim_md_trajectory", "truth")
  .skip_unless_built("make_prnode_design")
  .skip_unless_built("build_md_dataset", "level")
}

# One dataset of the new pipeline at p = 6, q = 1, n_M = 48, low noise; cached
# by (level, num_d, noise, seed, q, grid).
.fixture_dataset <- function(level = c("regression", "ode"),
                             num_d = 12,
                             noise_level_ratio = .fx_noise_low,
                             q = 1,
                             seed = .fx_seed,
                             times_e = .fx_times_e){
  level <- match.arg(level)
  .skip_unless_new_pipeline()
  key <- paste("dataset", level, num_d, noise_level_ratio, q, seed, length(times_e), sep = "_")
  if(!exists(key, envir = .fx)){
    assign(key,
           build_md_dataset(seed = seed,
                            num_d = num_d,
                            noise_level_ratio = noise_level_ratio,
                            truth_args = .fx_truth_args,
                            obs_time = .fx_obs_time,
                            times_e = times_e,
                            q = q,
                            level = level),
           envir = .fx)
  }
  get(key, envir = .fx)
}

# A q = 3 design on the q = 1 truth (plan S5: the misspecified-basis regime;
# identities only, never recovery).
.fixture_design_q3 <- function(num_d = 36,
                               level = c("regression", "ode")){
  level <- match.arg(level)
  ds <- .fixture_dataset(level = level, num_d = num_d)
  key <- paste("design_q3", level, num_d, sep = "_")
  if(!exists(key, envir = .fx)){
    x_list <- if(level == "regression") list(M = ds$truth$X_M, D = ds$truth$X_D) else NULL
    assign(key,
           make_prnode_design(m_obs_time = ds$raw$m_obs_time,
                              Y_M = ds$raw$Y_M,
                              d_obs_time = ds$raw$d_obs_time,
                              Y_D = ds$raw$Y_D,
                              times_e = ds$raw$times_e,
                              q = 3,
                              x_list = x_list),
           envir = .fx)
  }
  get(key, envir = .fx)
}

# A fit of the reference solver on a fixture, cached by its settings.
.fixture_fit <- function(design,
                         lambda,
                         mu = c(S = 1, M = 1, D = 1),
                         key = NULL,
                         ...){
  .skip_unless_built("fit_prnode_path")
  if(is.null(key)) return(fit_prnode_path(design, lambda = lambda, mu = mu, ...))
  key <- paste("fit", key, paste(signif(lambda, 8), collapse = ","), paste(mu, collapse = ","), sep = "_")
  if(!exists(key, envir = .fx)){
    assign(key, fit_prnode_path(design, lambda = lambda, mu = mu, ...), envir = .fx)
  }
  get(key, envir = .fx)
}

# A decreasing lambda path from the pooled lambda_max of a design.
.lambda_path <- function(design,
                         mu = c(S = 1, M = 1, D = 1),
                         num_lambda = 30,
                         num_decades = 2){
  .skip_unless_built("lambda_max_prnode")
  lambda_max <- lambda_max_prnode(design, mu = mu, bool_pooled = TRUE)
  lambda_max*10^seq(0, -num_decades, length.out = num_lambda)
}

# --- accessors on a prnode_path ---------------------------------------------

.path_mat <- function(path,
                      what,
                      l = 1){
  arr <- path[[what]]
  stopifnot(length(dim(arr)) == 3)
  matrix(arr[, , l], nrow = dim(arr)[1], ncol = dim(arr)[2])
}

.path_net <- function(path,
                      cell,
                      l = 1){
  arr <- path$networks[[cell]]
  stopifnot(length(dim(arr)) == 3)
  matrix(as.logical(arr[, , l]), nrow = dim(arr)[1], ncol = dim(arr)[2])
}

# the stacked coefficient vector of protein j in the column order of the design
.path_beta <- function(path,
                       j,
                       l = 1){
  c(.path_mat(path, "s", l)[j, ], .path_mat(path, "d_M", l)[j, ], .path_mat(path, "d_D", l)[j, ])
}

# an adjacency of any storage as logical p x p
.as_adj <- function(net_mat){
  matrix(as.logical(net_mat != 0), nrow = nrow(net_mat), ncol = ncol(net_mat))
}

# --- the design's groups ------------------------------------------------------

# columns of each group in `index` order, with its thin QR
.design_groups <- function(design){
  labels <- unique(design$index)
  lapply(labels, function(g){
    cols <- which(design$index == g)
    decomp <- qr(design$x[, cols, drop = FALSE])
    list(label = g,
         cols = cols,
         Q = qr.Q(decomp),
         R = qr.R(decomp),
         type = if(g %in% design$M_groups) "M" else if(g %in% design$D_groups) "D" else "S")
  })
}

.group_weights <- function(design,
                           mu = c(S = 1, M = 1, D = 1)){
  mu <- .mu_vec(mu)
  q <- design$q
  sapply(.design_groups(design), function(gr) sqrt(q)*mu[[gr$type]])
}

.mu_vec <- function(mu){
  if(is.null(names(mu))) names(mu) <- c("S", "M", "D")
  mu
}

# the reference objective of protein j at the stacked coefficient `beta`
# (original coordinates): 1/2 ||y_j - x beta||^2 + lambda sqrt(q) sum_g mu_g ||x_g beta_g||
.objective <- function(design,
                       j,
                       beta,
                       lambda,
                       mu = c(S = 1, M = 1, D = 1)){
  mu <- .mu_vec(mu)
  resid <- design$y[, j] - c(design$x %*% beta)
  groups <- .design_groups(design)
  penalty <- sum(sapply(groups, function(gr){
    mu[[gr$type]]*sqrt(sum((design$x[, gr$cols, drop = FALSE] %*% beta[gr$cols])^2))
  }))
  0.5*sum(resid^2) + lambda*sqrt(design$q)*penalty
}

# KKT violation of protein j at `beta`, in the orthonormalized coordinates
.kkt_violation <- function(design,
                           j,
                           beta,
                           lambda,
                           mu = c(S = 1, M = 1, D = 1)){
  mu <- .mu_vec(mu)
  resid <- design$y[, j] - c(design$x %*% beta)
  groups <- .design_groups(design)
  max(sapply(groups, function(gr){
    w <- sqrt(design$q)*mu[[gr$type]]
    b <- c(gr$R %*% beta[gr$cols])
    z <- c(crossprod(gr$Q, resid))
    if(sqrt(sum(b^2)) > 0){
      sqrt(sum((z - lambda*w*b/sqrt(sum(b^2)))^2))
    } else {
      max(0, sqrt(sum(z^2)) - lambda*w)
    }
  }))
}

# --- independent solvers (oracles) ------------------------------------------

# FISTA on the orthonormalized problem min 1/2 ||y - sum_g Q_g b_g||^2 +
# lambda sum_g w_g ||b_g||, an implementation unrelated to cyclic descent.
# Returns beta in the original coordinates.
.fista_group_lasso <- function(design,
                               j,
                               lambda,
                               mu = c(S = 1, M = 1, D = 1),
                               max_iter = 20000,
                               tol = 1e-13){
  groups <- .design_groups(design)
  w_vec <- .group_weights(design, mu)
  Q <- do.call(cbind, lapply(groups, function(gr) gr$Q))
  sizes <- sapply(groups, function(gr) length(gr$cols))
  starts <- cumsum(c(0, sizes[-length(sizes)]))
  y <- design$y[, j]
  lip <- max(svd(Q, nu = 0, nv = 0)$d)^2
  b <- rep(0, ncol(Q))
  z <- b
  t_old <- 1
  prox <- function(v){
    for(g in seq_along(groups)){
      idx <- starts[g] + seq_len(sizes[g])
      nv <- sqrt(sum(v[idx]^2))
      shrink <- max(0, 1 - lambda*w_vec[g]/(lip*nv))
      v[idx] <- if(nv > 0) shrink*v[idx] else 0
    }
    v
  }
  for(it in seq_len(max_iter)){
    grad <- -c(crossprod(Q, y - Q %*% z))
    b_new <- prox(z - grad/lip)
    t_new <- (1 + sqrt(1 + 4*t_old^2))/2
    z <- b_new + ((t_old - 1)/t_new)*(b_new - b)
    done <- sqrt(sum((b_new - b)^2)) <= tol*(1 + sqrt(sum(b_new^2)))
    b <- b_new
    t_old <- t_new
    if(done && it > 50) break
  }
  beta <- numeric(ncol(design$x))
  for(g in seq_along(groups)){
    idx <- starts[g] + seq_len(sizes[g])
    beta[groups[[g]]$cols] <- backsolve(groups[[g]]$R, b[idx])
  }
  beta
}

# the single-cell standardized group lasso of one cell (plan R9's oracle):
# min 1/2 ||y_c - sum_k Psi_ck theta_k||^2 + weight sum_k ||Psi_ck theta_k||,
# solved by FISTA in the orthonormalized coordinates. `x_cell` is the cell's
# rows of the deviation block (already carrying sqrt(omega)), `y_cell` its rows
# of y; `index_cell` the regulator label of each column.
.single_cell_group_lasso <- function(x_cell,
                                     y_cell,
                                     index_cell,
                                     weight,
                                     max_iter = 20000,
                                     tol = 1e-13){
  labels <- unique(index_cell)
  groups <- lapply(labels, function(g){
    cols <- which(index_cell == g)
    decomp <- qr(x_cell[, cols, drop = FALSE])
    list(cols = cols, Q = qr.Q(decomp), R = qr.R(decomp))
  })
  Q <- do.call(cbind, lapply(groups, function(gr) gr$Q))
  sizes <- sapply(groups, function(gr) length(gr$cols))
  starts <- cumsum(c(0, sizes[-length(sizes)]))
  lip <- max(svd(Q, nu = 0, nv = 0)$d)^2
  b <- rep(0, ncol(Q))
  z <- b
  t_old <- 1
  for(it in seq_len(max_iter)){
    grad <- -c(crossprod(Q, y_cell - Q %*% z))
    v <- z - grad/lip
    for(g in seq_along(groups)){
      idx <- starts[g] + seq_len(sizes[g])
      nv <- sqrt(sum(v[idx]^2))
      v[idx] <- if(nv > 0) max(0, 1 - weight/(lip*nv))*v[idx] else 0
    }
    t_new <- (1 + sqrt(1 + 4*t_old^2))/2
    z <- v + ((t_old - 1)/t_new)*(v - b)
    done <- sqrt(sum((v - b)^2)) <= tol*(1 + sqrt(sum(v^2)))
    b <- v
    t_old <- t_new
    if(done && it > 50) break
  }
  theta <- numeric(ncol(x_cell))
  for(g in seq_along(groups)){
    idx <- starts[g] + seq_len(sizes[g])
    theta[groups[[g]]$cols] <- backsolve(groups[[g]]$R, b[idx])
  }
  theta
}

# --- hand reconstruction of one cell from `raw` ------------------------------

# Rebuilds what `make_prnode_design()` sees for one cell from the dataset's
# `raw`, `truth` and `meta` alone: the trajectory the design was built from
# (true at the regression level, smoothed at the ODE level), the integrated
# basis, Z and the response in the design's units (divided by
# `response_scale`). Used wherever a test needs the unprofiled Psi or Z.
# With `bool_scale_regulators = TRUE` the trajectory is divided by
# `response_scale` as well before the basis is built, so that regulator k's
# m-th power is scaled by s_k^m: the design's units (plan D7(c), F7). The
# default keeps Psi in raw units for the hand computations in raw units
# (Q1, the smoother tests).
.hand_cell <- function(ds,
                       cell = c("M", "D"),
                       bool_scale_regulators = FALSE,
                       q = ds$public$q,
                       spar = 0.5){
  cell <- match.arg(cell)
  .skip_unless_built("make_integrated_basis")
  raw <- ds$raw
  bool_regression <- identical(ds$meta$level, "regression")
  if(cell == "M"){
    obs_time <- raw$m_obs_time
    Y <- raw$Y_M
    grid <- raw$times_e
    X <- if(bool_regression) ds$truth$X_M else smooth_md_trajectory(obs_time, Y, grid, spar = spar)
  } else {
    obs_time <- raw$d_obs_time
    Y <- raw$Y_D
    grid <- make_daughter_grid(raw$birth_time, raw$times_e)
    X <- if(bool_regression) ds$truth$X_D else smooth_md_trajectory(obs_time, Y, grid, spar = spar)
  }
  scale_vec <- ds$public$response_scale
  if(is.null(scale_vec)) scale_vec <- rep(1, ncol(Y))
  if(bool_scale_regulators) X <- sweep(X, 2, scale_vec, "/")
  Psi <- make_integrated_basis(X, grid, obs_time, q)
  Z <- cbind(1, obs_time - obs_time[1])
  list(grid = grid,
       obs_time = obs_time,
       Psi = Psi,
       qr_Z = qr(Z),
       X = X,
       Y = Y,
       y = sweep(Y, 2, scale_vec, "/"),
       Z = Z)
}

# --- scoring ------------------------------------------------------------------

# Matthews correlation over the off-diagonal positions (plan N3(b)); 0 when
# the denominator vanishes, as the pilot and evaluate.R score it.
.mcc <- function(est_mat,
                 true_mat){
  est <- .as_adj(est_mat)
  tru <- .as_adj(true_mat)
  off <- row(est) != col(est)
  tp <- sum(est[off] & tru[off])
  fp <- sum(est[off] & !tru[off])
  fn <- sum(!est[off] & tru[off])
  tn <- sum(!est[off] & !tru[off])
  denom <- sqrt(as.numeric(tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
  if(denom == 0) return(0)
  (tp*tn - fp*fn)/denom
}

# --- the scalar fixture of the derivation notes (plan F6, R13) -------------------------------

# the one-regulator fixture of window_check.R: 48 mother frames over
# three sine cycles and the daughter's last 12 frames (the last quarter cycle
# and a half)
.scalar_fixture <- function(){
  n_M <- 48
  n_D <- 12
  list(n_M = n_M,
       n_D = n_D,
       psi_M = sin(seq(0, 6*pi, length.out = n_M)) + 0.5,
       psi_D = sin(seq(4.5*pi, 6*pi, length.out = n_D)) + 0.5)
}

# the fork on the hand-built scalar design (window_check.R's
# fit_scalar), four-block "current" or two-block "recommended"
.fork_scalar_fit <- function(fx,
                             y_M,
                             y_D,
                             lambda,
                             design = c("current", "recommended"),
                             mscale = 1,
                             dscale = 1){
  design <- match.arg(design)
  n_M <- fx$n_M
  n_D <- fx$n_D
  N <- n_M + n_D
  if(design == "current"){
    w <- c(rep(1/sqrt(2*N), N), rep(1/sqrt(2*n_M), n_M), rep(1/sqrt(2*n_D), n_D))
    x_S <- c(fx$psi_M, fx$psi_D, fx$psi_M, fx$psi_D)
    x_M <- c(rep(0, N), fx$psi_M, rep(0, n_D))
    x_D <- c(rep(0, N), rep(0, n_M), fx$psi_D)
    y <- c(y_M, y_D, y_M, y_D)
  } else {
    w <- rep(1, N)
    x_S <- c(fx$psi_M, fx$psi_D)
    x_M <- c(fx$psi_M, rep(0, n_D))
    x_D <- c(rep(0, n_M), fx$psi_D)
    y <- c(y_M, y_D)
  }
  x <- cbind(x_S, x_M, x_D)*w
  fit <- grplasso.default(x, y*w, index = 1:3, M_groups = 2, D_groups = 3,
                          center = FALSE, lambda = lambda, model = LinReg(),
                          mscale = mscale, dscale = dscale, residualize = FALSE,
                          control = grpl.control(tol = 1e-12, max.iter = 5000, trace = 0))
  b <- c(fit$coefficients)
  c(s = b[1], d_M = b[2], d_D = b[3], theta_M = b[1] + b[2], theta_D = b[1] + b[3])
}

# window_check's grid search for the first lambda with a nonzero coefficient
.fork_scalar_lambda_max <- function(fx,
                                    y_M,
                                    y_D,
                                    design,
                                    mscale = 1,
                                    dscale = 1){
  grid <- 10^seq(4, -4, length.out = 400)
  for(l in grid){
    f <- .fork_scalar_fit(fx, y_M, y_D, l, design, mscale, dscale)
    if(any(abs(f[1:3]) > 1e-10)) return(l)
  }
  NA_real_
}

# --- the fork called as section 6.3 of the derivation notes (plan G1, G2) -------------------

# the four-column fork design of the recommended loss: Z_M and Z_D zero-padded
# as unpenalized columns, then the 3p groups of `design$x`
.fork_design <- function(design){
  n_M <- design$n_M
  n_D <- design$n_D
  t_M <- design$m_obs_time
  t_D <- design$d_obs_time
  omega <- unname(design$omega)
  Z <- rbind(cbind(sqrt(omega[1])*cbind(1, t_M - t_M[1]), 0, 0),
             cbind(0, 0, sqrt(omega[2])*cbind(1, t_D - t_D[1])))
  list(x = cbind(Z, design$x),
       index = c(rep(NA, 4), design$index),
       M_groups = design$M_groups,
       D_groups = design$D_groups,
       n = n_M + n_D)
}

# one protein through the fork at the settings of section 6.3 of the derivation notes; returns
# the stacked (s, d_M, d_D) coefficients in the design's column order
.fork_fit_protein <- function(design,
                              j,
                              lambda,
                              mu = c(S = 1, M = 1, D = 1),
                              tol = 1e-12){
  mu <- .mu_vec(mu)
  fd <- .fork_design(design)
  lambda_fork <- 2*sqrt(fd$n)*mu[["S"]]*lambda
  fit <- suppressWarnings(
    grplasso.default(x = fd$x, y = design$y[, j], index = fd$index,
                     M_groups = fd$M_groups, D_groups = fd$D_groups,
                     lambda = lambda_fork, model = LinReg(), center = FALSE,
                     standardize = TRUE, residualize = TRUE,
                     mscale = mu[["M"]]/mu[["S"]], dscale = mu[["D"]]/mu[["S"]],
                     control = grpl.control(tol = tol, max.iter = 20000, trace = 0)))
  coef_vec <- c(fit$coefficients[, ncol(fit$coefficients)])
  list(gamma = coef_vec[1:4], beta = coef_vec[-(1:4)], fit = fit)
}

# --- the p = 1 regression-level fixture and its exact solver (plan R4, R13, F6) ---

# the exact scalar solution of the three-group problem for p = 1, q = 1 by
# enumeration of the sign patterns (Proposition 5 of the derivation notes `prop:split`: the optimal
# scaffold is one of {0, u, v}, so at most two of s, d_M, d_D are nonzero)
.scalar_exact <- function(design, lambda, mu = c(S = 1, M = 1, D = 1)){
  mu <- .mu_vec(mu)
  x <- design$x
  y <- design$y[, 1]
  norms <- sqrt(colSums(x^2))
  W <- lambda*c(mu[["S"]]*norms[1], mu[["M"]]*norms[2], mu[["D"]]*norms[3])
  objective <- function(beta) 0.5*sum((y - x %*% beta)^2) + sum(W*abs(beta))
  patterns <- list(integer(0), 1, 2, 3, c(1, 2), c(1, 3), c(2, 3))
  best <- list(value = objective(c(0, 0, 0)), beta = c(0, 0, 0))
  for(free in patterns){
    if(length(free) == 0) next
    signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(free))))
    for(r in seq_len(nrow(signs))){
      sg <- signs[r, ]
      xf <- x[, free, drop = FALSE]
      rhs <- c(crossprod(xf, y)) - W[free]*sg
      b <- solve(crossprod(xf), rhs)
      if(all(sign(b) == sg)){
        beta <- c(0, 0, 0)
        beta[free] <- b
        val <- objective(beta)
        if(val < best$value) best <- list(value = val, beta = beta)
      }
    }
  }
  best
}

# a hand-built p = 1 regression-level design: mother with a self-edge, daughter
# without (plan R13's fixture), on the coarse grid. x0 = 10 (review pass 6):
# at x0 = 3 the profiled mother column had norm 0.88 against a threshold of
# 2.6 for the edge of strength 3, so the mother's edge was inactive in 149 of
# R13's 400 draws; at x0 = 10 the norm is 3.9 and the edge is always active.
.scalar_design <- function(sigma = 0.5, seed = 1, a = -3, b_M = 3, b_D = 2, x0 = 10, omega = c(1/sigma^2, 1/sigma^2)){
  .skip_unless_built("sim_ode_trajectory")
  .skip_unless_built("make_prnode_design")
  grid <- .fx_times_coarse
  obs_time <- .fx_obs_time
  birth_time <- obs_time[37]
  X_M <- sim_ode_trajectory(matrix(a, 1, 1), b_M, x0, grid, q = 1)
  d_grid <- grid[grid >= birth_time - 1e-12]
  X_D <- sim_ode_trajectory(matrix(0, 1, 1), b_D, X_M[which(abs(grid - birth_time) < 1e-10), ], d_grid, q = 1)
  m_idx <- sapply(obs_time, function(tt) which(abs(grid - tt) < 1e-10))
  d_obs <- obs_time[37:48]
  d_idx <- sapply(d_obs, function(tt) which(abs(d_grid - tt) < 1e-10))
  set.seed(seed)
  Y_M <- X_M[m_idx, , drop = FALSE] + stats::rnorm(48, sd = sigma)
  Y_D <- X_D[d_idx, , drop = FALSE] + stats::rnorm(12, sd = sigma)
  design <- make_prnode_design(m_obs_time = obs_time, Y_M = Y_M, d_obs_time = d_obs, Y_D = Y_D, times_e = grid,
                               q = 1, omega = omega, x_list = list(M = X_M, D = X_D), bool_scale_response = FALSE)
  list(design = design, X_M = X_M, X_D = X_D, Y_M = Y_M, Y_D = Y_D, obs_time = obs_time, d_obs = d_obs, grid = grid)
}


# --- fit_prnode() helpers (plan sections 11-12) --------------------------------

# Further assumed shapes (see the header): `predict.prnode_fit(object, cell)`
# needs no design argument, so a `prnode_fit` must keep what it needs (the
# plan's section 9 lists "the design settings" among the path's elements);
# `[.prnode_fit` takes one lambda index.

# coefficients fitted on the scaled response back in the original units: the
# coefficient of regulator k's m-th power in equation j is multiplied by
# s_j / s_k^m (plan D7(c), F7)
.unscale_coef <- function(theta_scaled, s_vec, q){
  p <- nrow(theta_scaled)
  out <- theta_scaled
  for(j in seq_len(p)) for(k in seq_len(p)) for(m in seq_len(q)){
    out[j, (k - 1)*q + m] <- theta_scaled[j, (k - 1)*q + m]*s_vec[j]/s_vec[k]^m
  }
  out
}

# an unscaled public from a dataset's raw (bool_scale_response = FALSE)
.unscaled_public <- function(ds, ...){
  .skip_unless_built("make_prnode_design")
  x_list <- if(identical(ds$meta$level, "regression")) list(M = ds$truth$X_M, D = ds$truth$X_D) else NULL
  make_prnode_design(m_obs_time = ds$raw$m_obs_time, Y_M = ds$raw$Y_M, d_obs_time = ds$raw$d_obs_time,
                     Y_D = ds$raw$Y_D, times_e = ds$raw$times_e, q = ds$public$q, x_list = x_list,
                     bool_scale_response = FALSE, ...)
}

# the daughter-oracle lambda index of a path fit: the point of highest
# daughter MCC (ties to the sparser, larger lambda)
.oracle_index <- function(fit, truth){
  mcc_vec <- sapply(seq_len(dim(fit$d_network_est)[3]), function(l) .mcc(fit$d_network_est[, , l], truth$d_network))
  which.max(mcc_vec)
}

# run a plotting call on a null device
.on_null_device <- function(expr){
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(expr)
}

# one real mother-daughter pair in the `raw` format, from the data files
# (skips when they are not reachable)
.real_pair_raw <- function(i = 1){
  pairs_path <- test_path("..", "..", "data", "all_MD_pairs.Rdata")
  data_path <- test_path("..", "..", "data", "LIC_data_201013.Rdata")
  skip_if_not(file.exists(pairs_path) && file.exists(data_path), "the real-data files are not reachable")
  env <- new.env()
  load(pairs_path, envir = env)
  load(data_path, envir = env)
  pairs <- env$all_MD_pairs
  mats <- env$mat_list_raw
  obs_time_full <- seq(0, 1, length.out = ncol(mats[[1]]))
  cell_mat <- function(idx) sapply(mats, function(m) m[idx, ])
  m_mat <- cell_mat(pairs$cell_idx_mother[i])
  d_mat <- cell_mat(pairs$cell_idx_daughter[i])
  m_rows <- which(apply(m_mat, 1, function(r) !anyNA(r)))
  d_rows <- which(apply(d_mat, 1, function(r) !anyNA(r)))
  protein_names <- c("Cdc10", "Stb3", "CLB5", "Whi5", "Xbp1", "Tup1")
  Y_M <- m_mat[m_rows, , drop = FALSE]
  Y_D <- d_mat[d_rows, , drop = FALSE]
  colnames(Y_M) <- colnames(Y_D) <- protein_names
  list(m_obs_time = obs_time_full[m_rows], Y_M = Y_M, d_obs_time = obs_time_full[d_rows], Y_D = Y_D,
       birth_time = obs_time_full[d_rows[1]], times_e = NULL, protein_names = protein_names)
}

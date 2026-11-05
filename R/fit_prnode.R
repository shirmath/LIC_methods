#' Fit the PRNODE networks of one mother-daughter pair
#'
#' The user-facing wrapper of the recommended estimator:
#' builds the two-block design with `make_prnode_design()` when given the raw
#' observations, fits `fit_prnode_path()` along `lambda`, and reads the networks
#' and coefficients back in the response's original units.
#'
#' @details
#' Units. The design divides protein `j` by its pooled
#' standard deviation `s_j` as response and as regulator (see
#' `make_prnode_design()` for the three steps of scaling and weighting), so
#' the solver's coefficient of regulator `k`'s `m`-th power in equation `j` is
#' in units of `s_j / s_k^m`. This function multiplies it back by `s_j / s_k^m`,
#' and the nuisance coefficients (initial state and slope) by `s_j`, so every
#' element named `*_est`, `intercepts` and `*_init_conditions` is in the units
#' of `Y_M` and `Y_D`; `path` keeps the solver's scaled units.
#'
#' Shape. For a single `lambda`, networks are `p x p`
#' logical matrices, coefficients `p x pq` matrices, `intercepts` a `2 x p`
#' matrix with rows `mother` and `daughter`, and the initial conditions
#' length-`p` vectors. For a path of `L` values every object gains a third
#' dimension of length `L` named by the `lambda` values (`intercepts` becomes
#' `2 x p x L`, the initial conditions `p x L`), and `fit[l]` returns the
#' single-`lambda` shape at `lambda[l]`.
#'
#' With a ready `public`, `bool_scale_response` (when given) must agree with
#' how the design was built, and a user `omega` reweights the design's rows
#' from the weights it was built with.
#'
#' @param public the design list of `make_prnode_design()` (or a dataset's
#' `public`); `NULL` builds it from `m_obs_time`, `Y_M`, `d_obs_time`, `Y_D`,
#' `times_e` and `q`.
#' @param lambda the penalty path, decreasing (one value allowed); or `NULL`
#' (the default), which chooses one value by `tune_prnode()` at `mu` with its
#' default settings (so the daughter must have at least 4 frames, and the
#' tuning's entry checks and warnings apply); or a `prnode_tuning` from
#' `tune_prnode()`, whose `lambda` and `mu` are used (a `mu` given here must
#' agree with the tuning's).
#' @param bool_scale_response passed to `make_prnode_design()`.
#' @param d_obs_time,Y_D the daughter's frames and observations (when `public`
#' is `NULL`).
#' @param edge_tol,max_iter,tol passed to `fit_prnode_path()`.
#' @param m_obs_time,Y_M the mother's frames and observations (when `public`
#' is `NULL`).
#' @param mu multipliers `c(S, M, D)` of the scaffold and deviation penalties.
#' @param noise_method passed to `make_prnode_design()`.
#' @param omega `NULL` (estimated) or a length-2 vector `c(M, D)` of loss
#' weights.
#' @param q degree of the polynomial basis (when `public` is `NULL`).
#' @param spar,spar_grid the smoothing level of the design when `public` is
#' `NULL`. With `spar = NULL` and `lambda = NULL` (the default) the level is
#' chosen from `spar_grid` together with the penalty by `tune_prnode_spar()`;
#' with a numeric `spar`, or a given `lambda`, the design is built at `spar`
#' (0.5 when `NULL`). Ignored when `public` is given: its own smoothing is
#' used.
#' @param times_e the fine grid (when `public` is `NULL`); `NULL` derives it
#' from the mother's frames.
#' @param verbose a non-negative number; `0` prints nothing. `1` reports this
#' function's phases (building the design, fitting, reading out). Each
#' further level is passed to `fit_prnode_path()` one lower, so `2` names each
#' protein, `3` reports each `lambda` and `4` each sweep (see `fit_prnode_path()`).
#' @returns an object of class `prnode_fit`: `converged` (`p x L` logical),
#' `design` (the design without `x` and `y`), `d_init_conditions`,
#' `d_network_est`, `d_ode_est`, `discount` (list `M`, `D`), `intercepts`,
#' `lambda`, `m_init_conditions`, `m_network_est`, `m_ode_est`, `mu`, `omega`
#' (named `M`, `D`), `path` (the `prnode_path`), `response_scale`,
#' `s_network_est`, `s_ode_est`; shapes as in Details; and, only when the
#' penalty was chosen by `tune_prnode()` (a `NULL` or `prnode_tuning`
#' `lambda`), `tuning`, that `prnode_tuning` (not carried through `fit[l]`).
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' mu <- c(S = 1, M = 0.9, D = 0.9)
#' lambda_max <- lambda_max_prnode(ds$public, mu = mu, bool_pooled = TRUE)
#' lambda_path <- lambda_max*10^seq(0, -2, length.out = 10)
#' fit <- fit_prnode(ds$public, lambda = lambda_path, mu = mu)
#' print(fit)
#' # the networks at the fifth lambda: TRUE where edge j <- k is present
#' fit$m_network_est[, , 5]
#' fit$d_network_est[, , 5]
#' # one point of the path, in the single-lambda shape
#' fit_5 <- fit[5]
#' round(fit_5$d_ode_est, 2)
#' dim(predict(fit_5, cell = "D"))
#' # the penalty chosen by tune_prnode() (a short path keeps the example quick)
#' tuning <- tune_prnode(ds$public, num_lambda = 20, num_decades = 2)
#' fit_tuned <- fit_prnode(ds$public, lambda = tuning)
#' fit_tuned$tuning$pick_idx
#' @export
fit_prnode <- function(public = NULL,
                      lambda = NULL,
                      bool_scale_response = TRUE,
                      d_obs_time = NULL,
                      edge_tol = 1e-8,
                      max_iter = 1e4,
                      m_obs_time = NULL,
                      mu = c(S = 1, M = 1, D = 1),
                      noise_method = c("diff2", "spline"),
                      omega = NULL,
                      q = 1,
                      spar = NULL,
                      spar_grid = c(0, 0.25, 0.5, 0.75),
                      times_e = NULL,
                      tol = 1e-10,
                      verbose = 0,
                      Y_D = NULL,
                      Y_M = NULL){
  noise_method <- match.arg(noise_method)
  if(length(verbose) != 1 || !is.numeric(verbose) || is.na(verbose) || verbose < 0){
    stop("`verbose` must be a single non-negative number (0 = silent, 1 = this function's phases, 2 and up = the solver's levels)")
  }
  if(is.null(public)){
    if(is.null(m_obs_time) || is.null(Y_M) || is.null(d_obs_time) || is.null(Y_D)){
      stop("without `public`, `m_obs_time`, `Y_M`, `d_obs_time` and `Y_D` are all needed to build the design")
    }
    if(is.null(lambda) && is.null(spar)){
      # raw observations and nothing fixed: the smoothing level and the
      # penalty are chosen together by the known-variance BIC
      if(verbose > 0) print("fit_prnode: choosing the smoothing level and the penalty with tune_prnode_spar()")
      chosen <- tune_prnode_spar(m_obs_time = m_obs_time, Y_M = Y_M, d_obs_time = d_obs_time, Y_D = Y_D,
                                 bool_scale_response = bool_scale_response, edge_tol = edge_tol,
                                 max_iter = max_iter, mu = mu, noise_method = noise_method, omega = omega,
                                 q = q, spar_grid = spar_grid, times_e = times_e, tol = tol,
                                 verbose = max(verbose - 1, 0))
      public <- chosen$public
      lambda <- chosen$tuning
      lambda$spar_criterion <- chosen$spar_criterion
    } else {
      if(verbose > 0) print("fit_prnode: building the two-block design from the raw observations")
      public <- make_prnode_design(m_obs_time = m_obs_time, Y_M = Y_M, d_obs_time = d_obs_time, Y_D = Y_D,
                                   bool_scale_response = bool_scale_response, noise_method = noise_method,
                                   omega = omega, q = q, spar = if(is.null(spar)) 0.5 else spar,
                                   times_e = times_e)
    }
  } else {
    if(!is.list(public) || is.null(public$x) || is.null(public$y)){
      stop("`public` must be the design list of make_prnode_design()")
    }
    if(!missing(bool_scale_response) && !is.null(public$response_scale)){
      bool_built_scaled <- any(public$response_scale != 1)
      if(bool_built_scaled != bool_scale_response){
        stop("`bool_scale_response = ", bool_scale_response, "` disagrees with the design, which was built ",
             if(bool_built_scaled) "with" else "without", " response scaling; rebuild the design or drop the argument")
      }
    }
    if(!is.null(omega)){
      if(length(omega) != 2 || any(omega <= 0)) stop("`omega` must be a length-2 positive vector c(M, D)")
      if(is.null(public$omega)) stop("`public` carries no `omega` to reweight from")
      m_rows <- seq_len(public$n_M)
      d_rows <- public$n_M + seq_len(public$n_D)
      factor_M <- sqrt(omega[1]/public$omega[[1]])
      factor_D <- sqrt(omega[2]/public$omega[[2]])
      public$x[m_rows, ] <- factor_M*public$x[m_rows, ]
      public$x[d_rows, ] <- factor_D*public$x[d_rows, ]
      public$y[m_rows, ] <- factor_M*public$y[m_rows, ]
      public$y[d_rows, ] <- factor_D*public$y[d_rows, ]
      public$omega <- c(M = unname(omega[1]), D = unname(omega[2]))
    }
  }
  # the penalty: given, chosen here by tune_prnode(), or read from a tuning
  # (after the omega reweighting above, so a penalty chosen here sees the
  # design that is fitted)
  tuning <- NULL
  if(is.null(lambda)){
    if(verbose > 0) print("fit_prnode: no `lambda` given; choosing it with tune_prnode()")
    tuning <- tune_prnode(public, edge_tol = edge_tol, max_iter = max_iter, mu = mu, tol = tol,
                          verbose = max(verbose - 1, 0))
    lambda <- tuning$lambda
  } else if(inherits(lambda, "prnode_tuning")){
    tuning <- lambda
    if(!missing(mu) && any(abs(.mu_vec(mu) - tuning$mu) > 1e-12)){
      stop("`mu` = (", paste(signif(mu, 3), collapse = ", "), ") disagrees with the tuning's mu = (",
           paste(signif(tuning$mu, 3), collapse = ", "), "); drop `mu` or pass the tuning's multipliers")
    }
    mu <- tuning$mu
    lambda <- tuning$lambda
  }
  if(verbose > 0){
    print(paste0("fit_prnode: fitting ", length(lambda), " lambda value", if(length(lambda) != 1) "s", " for ",
                 ncol(as.matrix(public$y)), " protein", if(ncol(as.matrix(public$y)) != 1) "s"))
  }
  path <- fit_prnode_path(public, lambda = lambda, edge_tol = edge_tol, max_iter = max_iter, mu = mu, tol = tol,
                     verbose = max(verbose - 1, 0))
  if(verbose > 0) print("fit_prnode: reading the networks and coefficients back in the response's original units")
  fit <- .prnode_fit_from_path(path)
  if(!is.null(tuning)) fit$tuning <- tuning
  fit
}

#' Build a prnode_fit from a prnode_path (the read-out in original units)
#' @noRd
.prnode_fit_from_path <- function(path){
  settings <- path$design
  p <- dim(path$theta_M)[1]
  pq <- dim(path$theta_M)[2]
  L <- dim(path$theta_M)[3]
  q <- if(!is.null(settings$q)) settings$q else pq/p
  s_vec <- settings$response_scale
  if(is.null(s_vec)) s_vec <- rep(1, p)
  reg_scale <- settings$profile$M$regulator_scale
  p_reg <- pq/q
  if(is.null(reg_scale)) reg_scale <- if(length(s_vec) == p_reg) s_vec else rep(1, p_reg)
  # coefficient (j, (k - 1) q + m) is multiplied by s_j / s_k^m, with s_j the
  # response scale of row j and s_k the regulator scale of column block k
  unscale_mat <- matrix(1, p, pq)
  for(k in seq_len(p_reg)) for(m in seq_len(q)){
    unscale_mat[, (k - 1)*q + m] <- s_vec/reg_scale[k]^m
  }
  unscale_coef <- function(arr){
    out <- arr*as.vector(unscale_mat)
    if(L == 1) matrix(out[, , 1], p, pq) else out
  }
  as_network <- function(arr){
    if(L == 1) matrix(arr[, , 1], p, p_reg) else arr
  }
  lambda_names <- format(path$lambda, digits = 15)
  intercepts <- array(NA_real_, c(2, p, L), dimnames = list(c("mother", "daughter"), NULL, lambda_names))
  m_init <- matrix(NA_real_, p, L, dimnames = list(NULL, lambda_names))
  d_init <- m_init
  for(l in seq_len(L)){
    intercepts["mother", , l] <- path$gamma_M[2, , l]*s_vec
    intercepts["daughter", , l] <- path$gamma_D[2, , l]*s_vec
    m_init[, l] <- path$gamma_M[1, , l]*s_vec
    d_init[, l] <- path$gamma_D[1, , l]*s_vec
  }
  if(L == 1){
    intercepts <- matrix(intercepts[, , 1], 2, p, dimnames = list(c("mother", "daughter"), NULL))
    m_init <- m_init[, 1]
    d_init <- d_init[, 1]
  }
  structure(list(converged = path$converged,
                 design = settings,
                 d_init_conditions = d_init,
                 d_network_est = as_network(path$networks$D),
                 d_ode_est = unscale_coef(path$theta_D),
                 discount = path$discount,
                 intercepts = intercepts,
                 lambda = path$lambda,
                 m_init_conditions = m_init,
                 m_network_est = as_network(path$networks$M),
                 m_ode_est = unscale_coef(path$theta_M),
                 mu = path$mu,
                 omega = settings$omega,
                 path = path,
                 response_scale = s_vec,
                 s_network_est = as_network(path$networks$S),
                 s_ode_est = unscale_coef(path$s)),
            class = "prnode_fit")
}

#' Subset a prnode_path to some lambda indices
#' @noRd
.subset_path <- function(path,
                         i){
  out <- path
  for(nm in c("s", "d_M", "d_D", "theta_M", "theta_D", "gamma_M", "gamma_D")){
    out[[nm]] <- path[[nm]][, , i, drop = FALSE]
  }
  for(nm in c("S", "M", "D")) out$networks[[nm]] <- path$networks[[nm]][, , i, drop = FALSE]
  for(nm in c("objective", "gap", "kkt_max", "iterations", "converged")){
    out[[nm]] <- path[[nm]][, i, drop = FALSE]
  }
  out$lambda <- path$lambda[i]
  out
}

#' Extract one or more lambda values of a prnode_fit
#'
#' `fit[l]` returns the fit at `lambda[l]` in the single-`lambda` shape;
#' several indices return a path-shaped fit.
#'
#' @param x a `prnode_fit`.
#' @param i whole-number indices into `x$lambda`, or a logical vector of length
#' `length(x$lambda)`.
#' @returns a `prnode_fit`.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' mu <- c(S = 1, M = 0.9, D = 0.9)
#' lambda_max <- lambda_max_prnode(ds$public, mu = mu, bool_pooled = TRUE)
#' lambda_path <- lambda_max*10^seq(0, -2, length.out = 10)
#' fit <- fit_prnode(ds$public, lambda = lambda_path, mu = mu)
#' fit[5]
#' fit[c(2, 5)]$lambda
#' @export
`[.prnode_fit` <- function(x, i){
  L <- length(x$lambda)
  if(is.logical(i)){
    if(length(i) != L) stop("a logical `i` must have one entry per lambda value (", L, ")")
    i <- which(i)
  }
  if(!is.numeric(i) || length(i) < 1 || any(is.na(i)) || any(i != round(i)) || any(i < 1) || any(i > L)){
    stop("`i` must be whole-number indices into the ", L, " lambda value(s) of the fit, or a logical of that length")
  }
  i <- as.integer(i)
  .prnode_fit_from_path(.subset_path(x$path, i))
}

#' Print a prnode_fit
#'
#' Reports the dimensions, `lambda`, `mu`, whether the response was scaled, the
#' edge counts of the scaffold (`s_`), mother and daughter networks, and the
#' maximum discount per cell.
#'
#' @param x a `prnode_fit`.
#' @param ... ignored.
#' @returns `x`, invisibly.
#' @export
print.prnode_fit <- function(x, ...){
  settings <- x$design
  L <- length(x$lambda)
  p <- length(x$response_scale)
  q <- if(!is.null(settings$q)) settings$q else dim(x$path$theta_M)[2]/p
  cat("PRNODE fit (class prnode_fit)\n")
  cat("  proteins p =", p, " basis degree q =", q, "\n")
  cat("  frames: mother n_M =", settings$n_M, " daughter n_D =", settings$n_D, "\n")
  if(L == 1){
    cat("  lambda =", signif(x$lambda, 3), "\n")
  } else {
    cat("  lambda path of", L, "values from", signif(x$lambda[1], 3), "to", signif(x$lambda[L], 3), "\n")
  }
  cat("  mu (S, M, D) =", paste(signif(x$mu, 3), collapse = ", "), "\n")
  if(!is.null(x$tuning)){
    cat("  lambda chosen by tune_prnode(): pick at path index", x$tuning$pick_idx, "of",
        length(x$tuning$lambda_path), if(!is.na(x$tuning$spar)) paste0(" (spar ", x$tuning$spar, ")"), "\n")
  }
  cat("  response scaled by each protein's pooled sd:", if(any(x$response_scale != 1)) "yes" else "no", "\n")
  if(!is.null(x$omega)) cat("  loss weights omega (M, D) =", paste(signif(x$omega, 3), collapse = ", "), "\n")
  count <- function(net) if(L == 1) sum(net) else paste(sum(net[, , 1]), "to", sum(net[, , L]))
  cat("  edges (fitted-contribution rule): s_network_est", count(x$s_network_est),
      " m_network_est", count(x$m_network_est), " d_network_est", count(x$d_network_est), "\n")
  cat("  maximum discount: mother", signif(max(x$discount$M, na.rm = TRUE), 3),
      " daughter", signif(max(x$discount$D, na.rm = TRUE), 3), "\n")
  cat("  converged:", sum(x$converged), "of", length(x$converged), "(protein, lambda) fits\n")
  invisible(x)
}

#' Fitted trajectories of a prnode_fit at the frames
#'
#' `Z_c gamma_c + Psi_c theta_c` for the requested cell, in the response's
#' original units.
#'
#' @param object a `prnode_fit`.
#' @param cell `"M"` (mother) or `"D"` (daughter).
#' @param ... ignored.
#' @returns an `n_c x p` matrix for one `lambda`, or `n_c x p x L` for a path.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' mu <- c(S = 1, M = 0.9, D = 0.9)
#' lambda_max <- lambda_max_prnode(ds$public, mu = mu, bool_pooled = TRUE)
#' lambda_path <- lambda_max*10^seq(0, -2, length.out = 10)
#' fit <- fit_prnode(ds$public, lambda = lambda_path, mu = mu)
#' Y_D_hat <- predict(fit[5], cell = "D")
#' dim(Y_D_hat)
#' dim(ds$raw$Y_D)
#' @export
predict.prnode_fit <- function(object,
                               cell = c("M", "D"),
                               ...){
  cell <- match.arg(cell)
  prof <- object$design$profile[[cell]]
  if(is.null(prof)) stop("the fit carries no profile for cell ", cell, "; it was fitted on a design without one")
  path <- object$path
  gamma_arr <- if(cell == "M") path$gamma_M else path$gamma_D
  theta_arr <- if(cell == "M") path$theta_M else path$theta_D
  p <- dim(theta_arr)[1]
  L <- dim(theta_arr)[3]
  s_vec <- object$response_scale
  out <- array(NA_real_, c(prof$n, p, L))
  for(l in seq_len(L)){
    fitted_scaled <- prof$Z %*% matrix(gamma_arr[, , l], 2, p) + prof$Psi %*% t(matrix(theta_arr[, , l], p))
    out[, , l] <- sweep(fitted_scaled, 2, s_vec, "*")
  }
  if(L == 1) matrix(out[, , 1], prof$n, p) else out
}

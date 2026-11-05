# The inheritance measure of a mother-daughter pair: the one-frame flow-map
# similarity of the two regulatory systems (the definition, on true or fitted
# coefficients), its plug-in estimator at the pick of a tuned PRNODE fit, and the
# wrong-mother permutation test. The estimator has no fitted constants. The simulator stores the
# same measure on the two true systems as `truth$inheritance`
# (`build_md_dataset()`), so an estimate is read against the truth's number.

#' The Jacobian of the additive polynomial field at a state
#'
#' For the field `f_j(x) = sum_k sum_m pars[j, (k - 1) q + m] x_k^m` of
#' `sim_md_network()` (column `(k - 1) q + m` = regulator `k`'s `m`-th power),
#' the `p x p` matrix `df_j / dx_k` at `x`. At `q = 1` it is `pars` itself,
#' whatever `x`.
#'
#' @param pars `p x pq` coefficient matrix, row = equation.
#' @param q the polynomial degree.
#' @param x length-`p` state (ignored at `q = 1`).
#' @returns a `p x p` matrix.
#' @noRd
.jacobian_at <- function(pars,
                         q,
                         x){
  p <- nrow(pars)
  if(q == 1) return(matrix(pars, p, p))
  jac <- matrix(0, p, p)
  for(k in seq_len(p)){
    for(m in seq_len(q)){
      jac[, k] <- jac[, k] + m*pars[, (k - 1)*q + m]*x[k]^(m - 1)
    }
  }
  jac
}

#' The displacement a linear system produces in one step: `e^{hA} - I`
#'
#' `NA` (a matrix of `NA`) when the exponential overflows, which happens only
#' for very large coefficients (the dense end of a `q = 3` path, where the raw
#' cubic coefficients reach `1e5` and more).
#'
#' @noRd
.flow_displacement <- function(A,
                               h){
  p <- nrow(A)
  out <- tryCatch(as.matrix(Matrix::expm(Matrix::Matrix(h*A))) - diag(p),
                  error = function(e) matrix(NA_real_, p, p))
  if(!all(is.finite(out))) out[] <- NA_real_
  out
}

#' The per-power unscaling matrix of a design's coefficients
#'
#' The design divides protein `j` by its pooled sd `s_j` as a response and as
#' a regulator, so the solver's coefficient of regulator `k`'s `m`-th power in
#' equation `j` is in units of `s_j / s_k^m`; multiplying by this `p x pq`
#' matrix returns it to the units of the observations (what `fit_prnode()`
#' does for its `*_ode_est`).
#'
#' @noRd
.unscale_mat <- function(p,
                         q,
                         response_scale,
                         regulator_scale){
  out <- matrix(1, p, p*q)
  for(k in seq_len(p)){
    for(m in seq_len(q)){
      out[, (k - 1)*q + m] <- response_scale/regulator_scale[k]^m
    }
  }
  out
}

#' The one-frame flow-map similarity of two regulatory systems
#'
#' The package's definition of the inheritance of a mother-daughter pair, a
#' number in `[0, 1]`: with `A_M` and `A_D` the two cells' system matrices
#' and `h` one sampling interval of the movie,
#' `tau = 1 - ||e^{h A_M} - e^{h A_D}||_F / (||e^{h A_M} - I||_F + ||e^{h A_D} - I||_F)`,
#' the similarity of the displacements the two regulatory systems produce in
#' one frame, self-regulation (the diagonal) included. It is 1 if and only if
#' the two systems coincide (and 1 by convention when both are empty), 0 when
#' their one-frame displacements point in opposite directions, and to first
#' order in `h` it is `1 - ||A_M - A_D||_F / (||A_M||_F + ||A_D||_F)`, so
#' every coefficient is weighted by its size and the sampling interval cancels.
#' Two unrelated systems do not read 0 but a chance level near
#' `1 - 1 / sqrt(2)` (about 0.3), which is why the absence of inheritance is
#' judged by `inheritance_test()` against wrong-mother pairings rather than
#' against zero.
#'
#' At `q = 1` the system matrix of a cell is its coefficient matrix. At
#' `q > 1` it is the Jacobian of that cell's polynomial field at **its own
#' mean state**, `x_ref$M` for the mother and `x_ref$D` for the daughter (the
#' column means of each cell's frames, `colMeans(Y_M)` and `colMeans(Y_D)`),
#' so each system is read where its own cell lives. A single common state
#' would extrapolate a short daughter's quadratic field to states she never
#' visits; each cell's linearized system is smallest at its own mean, so no
#' single state is neutral for both.
#'
#' @param m_ode_pars the mother's `p x pq` coefficient matrix (row `j` = the
#' equation of protein `j`, column `(k - 1) q + m` = regulator `k`'s `m`-th
#' power): a truth's `m_ode_pars`, or a fit's `m_ode_est`.
#' @param d_ode_pars the daughter's, same shape.
#' @param h the horizon, one sampling interval of the movie in the units of
#' the frame times; must be positive.
#' @param q the polynomial degree of the two matrices (`pq / p`).
#' @param x_ref at `q > 1` a list with `M` and `D`, the length-`p` states the
#' mother's and the daughter's fields are linearized at (each cell's own
#' mean state); ignored at `q = 1`.
#' @returns one number in `[0, 1]`, or `NA` when a matrix exponential
#' overflows.
#' @examples
#' truth <- sim_md_network(p = 6, q = 1, num_shared = 6, num_mother_only = 2,
#'                         num_daughter_only = 2, seed_number = 1)
#' flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars, h = 1/47)
#' flow_map_similarity(truth$m_ode_pars, truth$m_ode_pars, h = 1/47)   # a clone: 1
#' # the first-order form, close to the exact value at a short horizon
#' 1 - norm(truth$m_ode_pars - truth$d_ode_pars, "F")/
#'   (norm(truth$m_ode_pars, "F") + norm(truth$d_ode_pars, "F"))
#' # a quadratic truth: each field linearized at its own cell's mean state
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4, q = 2), q = 2,
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' flow_map_similarity(ds$truth$m_ode_pars, ds$truth$d_ode_pars, h = 1/35, q = 2,
#'                     x_ref = list(M = colMeans(ds$raw$Y_M), D = colMeans(ds$raw$Y_D)))
#' @export
flow_map_similarity <- function(m_ode_pars,
                                d_ode_pars,
                                h,
                                q = 1,
                                x_ref = NULL){
  m_ode_pars <- as.matrix(m_ode_pars)
  d_ode_pars <- as.matrix(d_ode_pars)
  stopifnot(is.numeric(m_ode_pars), is.numeric(d_ode_pars),
            length(h) == 1, is.numeric(h), is.finite(h), h > 0,
            length(q) == 1, is.numeric(q), q >= 1, q == round(q))
  p <- nrow(m_ode_pars)
  if(any(dim(d_ode_pars) != dim(m_ode_pars))){
    stop("`m_ode_pars` and `d_ode_pars` must have the same dimensions (p x pq)")
  }
  if(ncol(m_ode_pars) != p*q) stop("`m_ode_pars` must be p x pq with q = ", q, ", not ", p, " x ", ncol(m_ode_pars))
  if(!all(is.finite(m_ode_pars)) || !all(is.finite(d_ode_pars))) stop("the coefficient matrices must be finite")
  x_ref <- .check_x_ref(x_ref, p, q, remedy = ", list(M = colMeans(Y_M), D = colMeans(Y_D))")
  .flow_map_similarity(.jacobian_at(m_ode_pars, q, x_ref$M), .jacobian_at(d_ode_pars, q, x_ref$D), h)
}

#' Check the two cells' reference states
#'
#' At `q = 1` no state is needed: `x_ref` is ignored and `list(M = NULL,
#' D = NULL)` is returned (`.jacobian_at()` ignores the state there). At
#' `q > 1` `x_ref` must be a list with `M` and `D`, each a finite numeric
#' state of length `p`; the two are returned as plain numeric vectors.
#'
#' @param remedy text appended to the error message, naming what to pass.
#' @returns a list `M`, `D`, each `NULL` (at `q = 1`) or a length-`p` numeric.
#' @noRd
.check_x_ref <- function(x_ref,
                         p,
                         q,
                         remedy = ""){
  if(q == 1) return(list(M = NULL, D = NULL))
  bool_state <- function(x) is.numeric(x) && length(x) == p && all(is.finite(x))
  bool_ok <- is.list(x_ref) && !is.data.frame(x_ref) && all(c("M", "D") %in% names(x_ref)) &&
    bool_state(x_ref$M) && bool_state(x_ref$D)
  if(!bool_ok){
    stop("at q = ", q, " `x_ref` must be a list with `M` and `D`, the length-", p,
         " states the mother's and the daughter's fields are linearized at (each cell's own mean state",
         remedy, ")")
  }
  list(M = as.numeric(x_ref$M), D = as.numeric(x_ref$D))
}

#' The similarity of two system matrices, no checks
#' @noRd
.flow_map_similarity <- function(A_M,
                                 A_D,
                                 h){
  e_M <- .flow_displacement(A_M, h)
  e_D <- .flow_displacement(A_D, h)
  if(anyNA(e_M) || anyNA(e_D)) return(NA_real_)
  denom <- sqrt(sum(e_M^2)) + sqrt(sum(e_D^2))
  if(denom == 0) return(1)
  min(1, max(0, 1 - sqrt(sum((e_M - e_D)^2))/denom))
}

#' The flow-map similarity of two true systems in standardized units
#'
#' The target of `inheritance_flow()`: the two systems (at `q > 1`, each
#' cell's Jacobian at the mean of its own noiseless states) expressed in the
#' units the solver sees, `S^{-1} A S` with `S` the diagonal of the pooled sd
#' of the two cells' noiseless states, then compared by the one-frame flow-map
#' similarity.
#'
#' @param y_m_noiseless,y_d_noiseless the two cells' noiseless frames
#' (`n x p`).
#' @returns one number in `[0, 1]`.
#' @noRd
.flow_map_similarity_std <- function(m_ode_pars,
                                     d_ode_pars,
                                     y_m_noiseless,
                                     y_d_noiseless,
                                     h,
                                     q = 1){
  p <- nrow(m_ode_pars)
  x_ref <- if(q > 1) list(M = colMeans(y_m_noiseless), D = colMeans(y_d_noiseless)) else NULL
  s_vec <- vapply(seq_len(p), function(j) stats::sd(c(y_m_noiseless[, j], y_d_noiseless[, j])), numeric(1))
  to_std <- function(A) sweep(sweep(A, 1, s_vec, "/"), 2, s_vec, "*")
  .flow_map_similarity(to_std(.jacobian_at(m_ode_pars, q, x_ref$M)), to_std(.jacobian_at(d_ode_pars, q, x_ref$D)), h)
}

#' Read the tuned path and the design out of a prnode_fit or a prnode_tuning
#' @noRd
.tuned_path_of <- function(fit){
  if(inherits(fit, "prnode_fit")){
    if(is.null(fit$tuning)){
      stop("`fit` carries no tuning: the estimator reads the tuned path, so fit with fit_prnode(public) ",
           "(no `lambda`) or with `lambda` = a prnode_tuning from tune_prnode()")
    }
    tuning <- fit$tuning
  } else if(inherits(fit, "prnode_tuning")){
    tuning <- fit
  } else {
    stop("`fit` must be a prnode_fit fitted with the tuning rule (from fit_prnode()), or a prnode_tuning ",
         "(from tune_prnode())")
  }
  if(!inherits(tuning, "prnode_tuning") || !inherits(tuning$path, "prnode_path")){
    stop("the tuning must be a prnode_tuning whose `path` is a prnode_path")
  }
  tuning
}

#' Estimate the inheritance of a fitted pair: the flow-map plug-in at the pick
#'
#' The estimator of `flow_map_similarity()` from a PRNODE fit whose penalty
#' was chosen by `tune_prnode()`: the two fitted systems at the tuned `lambda` (scaffold plus
#' deviation) are compared by the flow-map similarity, read in the
#' **standardized** units of the design (each protein divided by its pooled
#' sd, as the solver sees it), so that a protein's scale does not weight its
#' row of the system. It has no fitted constant and no window along the path:
#' the known-variance BIC pick is the selector's own point, and averaging the
#' plug-in over a window past it did not improve the estimate in simulations.
#'
#' The estimate is compressed toward the middle of `[0, 1]` by estimation
#' noise and biased low at high noise (the pick's bias grows with the
#' daughter's length), so it ranks pairs and feeds the permutation test
#' `inheritance_test()`; it is not calibrated to the scale of the truth. Its
#' target is the truth of the standardized systems at one sampling interval.
#'
#' At `q > 1` the fitted fields are polynomial and each cell's system is the
#' Jacobian of its field at **that cell's own mean observed state**: the mother's at
#' `colMeans(Y_M)` and the daughter's at `colMeans(Y_D)`, which
#' `make_prnode_design()` stores as the design's `mean_state`, each divided
#' by the design's `response_scale` because the solver's coefficients are
#' standardized. Reading both cells at one common state would make a short
#' daughter's estimate an artefact of extrapolating her quadratic field to
#' states she never visits. The target `truth$inheritance_std` of `build_md_dataset()` linearizes the two true
#' fields the same way, each at the mean of its own cell's noiseless states.
#'
#' @param fit a `prnode_fit` from `fit_prnode()` with its penalty chosen by
#' the tuning rule (`lambda = NULL` or a `prnode_tuning`), or the
#' `prnode_tuning` itself. The tuned path is read from `fit$tuning$path`.
#' @param bool_path return the plug-in along the whole path beside the
#' estimate.
#' @param h the horizon; `NULL` (the default) is one sampling interval of the
#' mother, the median spacing of her frames in the design.
#' @param x_ref `NULL` (the default: at `q > 1` each fitted field is
#' linearized at its own cell's mean observed state, the design's
#' `mean_state`) or, at `q > 1`, a list with `M` and `D`, the two length-`p`
#' states in the units of the observations, for a sensitivity reading; each
#' is divided by the design's `response_scale` to the standardized state.
#' Ignored at `q = 1`.
#' @returns with `bool_path = FALSE`, one number in `[0, 1]`. With
#' `bool_path = TRUE`, a list: `estimate` (the number), `h`, `path_df` (data
#' frame, one row per path point in the fit's order: `index`, `lambda`,
#' `similarity` (the plug-in; `NA` where a matrix exponential overflowed),
#' `bool_pick`) and `pick_idx`.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' fit <- fit_prnode(ds$public)          # lambda chosen by tune_prnode()
#' inheritance_flow(fit)
#' out <- inheritance_flow(fit, bool_path = TRUE)
#' round(out$path_df$similarity[seq_len(out$pick_idx)], 3)
#' @export
inheritance_flow <- function(fit,
                             bool_path = FALSE,
                             h = NULL,
                             x_ref = NULL){
  tuning <- .tuned_path_of(fit)
  stopifnot(length(bool_path) == 1, is.logical(bool_path), !is.na(bool_path))
  path <- tuning$path
  design <- path$design
  p <- dim(path$theta_M)[1]
  pq <- dim(path$theta_M)[2]
  L <- dim(path$theta_M)[3]
  q <- if(!is.null(design$q)) design$q else pq/p
  if(pq != p*q){
    stop("the flow map needs every protein as a response (", p, " responses, ", pq/q,
         " regulators): fit the full system, not a subset of responses")
  }
  if(is.null(h)){
    if(is.null(design$m_obs_time) || length(design$m_obs_time) < 2){
      stop("`h` cannot be read from the design (no mother frames); pass the sampling interval")
    }
    h <- stats::median(diff(design$m_obs_time))
  }
  stopifnot(length(h) == 1, is.numeric(h), is.finite(h), h > 0)
  if(q > 1){
    if(is.null(x_ref)){
      x_ref <- design$mean_state
      if(is.null(x_ref)){
        stop("the design carries no `mean_state` (a design built before package 0.13.0, or by hand): at q = ", q,
             " pass `x_ref = list(M = colMeans(Y_M), D = colMeans(Y_D))`")
      }
    }
    .check_x_ref(x_ref, p, q, remedy = "; the design's `mean_state` is the default")
  }
  response_scale <- design$response_scale
  if(is.null(response_scale)) response_scale <- rep(1, p)
  # the solver's coefficients are in standardized units, so each standardized
  # system is linearized at its own cell's standardized state (NULL stays
  # NULL at q = 1, where .jacobian_at() ignores the state)
  z_ref <- lapply(.check_x_ref(x_ref, p, q), function(x) if(is.null(x)) NULL else x/response_scale)
  pick_idx <- tuning$pick_idx

  # Only the pick enters the estimate; the rest of the path is evaluated only
  # when the path is asked for (the dense end can overflow at q > 1).
  idx_eval <- if(bool_path) seq_len(L) else pick_idx
  similarity_vec <- rep(NA_real_, L)
  similarity_vec[idx_eval] <- vapply(idx_eval, function(l){
    theta_M <- matrix(path$theta_M[, , l], p, pq)
    theta_D <- matrix(path$theta_D[, , l], p, pq)
    .flow_map_similarity(.jacobian_at(theta_M, q, z_ref$M), .jacobian_at(theta_D, q, z_ref$D), h)
  }, numeric(1))
  estimate <- similarity_vec[pick_idx]
  if(!is.finite(estimate)){
    stop("the flow map is not finite at the pick (path index ", pick_idx,
         "): a matrix exponential overflowed; shorten `h`")
  }
  stopifnot(estimate >= 0, estimate <= 1)
  if(!bool_path) return(estimate)
  list(estimate = estimate,
       h = h,
       path_df = data.frame(index = seq_len(L),
                            lambda = path$lambda,
                            similarity = similarity_vec,
                            bool_pick = seq_len(L) == pick_idx),
       pick_idx = pick_idx)
}

#' The permutation p-value of an inheritance estimate against its null
#'
#' `(1 + #{null >= estimate}) / (1 + length(null))`: the share of the
#' wrong-pairing estimates at or above the pair's own, with the pair itself
#' counted once, so the smallest attainable value is `1 / (1 + length(null))`
#' and the p-value is exact under exchangeability of the pair's estimate with
#' its wrong pairings' (Phipson and Smyth, 2010).
#'
#' @param estimate the pair's inheritance estimate (one number).
#' @param null_estimates the estimates of the wrong pairings of the same
#' daughter (a numeric vector, at least one value, no `NA`).
#' @returns one number in `(0, 1]`.
#' @examples
#' inheritance_pvalue(0.8, c(0.5, 0.6, 0.9, 0.4))   # 2 / 5
#' @export
inheritance_pvalue <- function(estimate,
                               null_estimates){
  stopifnot(length(estimate) == 1, is.numeric(estimate), is.finite(estimate),
            is.numeric(null_estimates), length(null_estimates) >= 1, !anyNA(null_estimates))
  (1 + sum(null_estimates >= estimate))/(1 + length(null_estimates))
}

#' The wrong-mother permutation test of inheritance
#'
#' Is the daughter's network inherited from her mother rather than from any
#' cell? The daughter is paired with each of `null_mothers`, every wrong
#' pairing is fitted jointly by PRNODE at the tuning rule exactly as the true
#' pair is (`fit_prnode()` on the raw observations with `design_args` and
#' `fit_args` and no `lambda`, so every pairing chooses its own smoothing
#' level and penalty by `tune_prnode_spar()` unless `design_args$spar` fixes
#' the level; the same settings for every pairing, which the exchangeability
#' of the p-value needs), and each fit is read by `inheritance_flow()` at its
#' own pick. Each pairing therefore costs about four path fits (one per level of
#' the default `spar_grid`); fix `design_args$spar` to cut that to one. The p-value is
#' `inheritance_pvalue()` of the true pair's estimate against the wrong
#' pairings'. Under the null that the daughter's network is unrelated to
#' which mother she is paired with, her true estimate is exchangeable with
#' the wrong ones, so the p-value is exact at the resolution
#' `1 / (1 + length(null_mothers))`. Two features of the statistic: the test
#' is invariant to any monotone recalibration of the estimate, and the null
#' must come from joint refits, because the scaffold of the joint fit is what
#' the estimate reads (separately fitted systems compare at a lower level).
#'
#' @param m_obs_time,Y_M the true mother's frames and `n_M x p` observations.
#' @param d_obs_time,Y_D the daughter's frames and `n_D x p` observations.
#' @param null_mothers a list of wrong mothers, each an `n x p` matrix of
#' observations on the frames `m_obs_time`, or a list with `Y_M` and its own
#' `m_obs_time`; the names of the list name the null estimates.
#' @param bool_keep_fits keep every wrong pairing's `prnode_fit` in the
#' return (`fits`).
#' @param design_args a named list of further arguments for the design,
#' passed to `fit_prnode()` (`noise_method`, `omega`, `spar`, `spar_grid`,
#' `bool_scale_response`), applied to every pairing; a given `spar` fixes the
#' smoothing level and skips the search.
#' @param fit `NULL` (the true pair is fitted here) or the true pair's
#' `prnode_fit` from `fit_prnode()` at the tuning rule on the design built
#' from the same arguments, `design_args`, `fit_args` and `mu`, which saves
#' one fit. With `fit` given and `times_e` `NULL`, the wrong pairings are
#' built on the grid of `fit$design$times_e`, so their integrals match the
#' true pair's to the last digit.
#' @param fit_args a named list of further arguments for `fit_prnode()`
#' (`edge_tol`, `tol`, `max_iter`), applied to every pairing.
#' @param h the horizon of `inheritance_flow()`; `NULL` is one sampling
#' interval of the true mother.
#' @param mu the multipliers of `fit_prnode()`, the same for every pairing.
#' @param q the polynomial degree of the design.
#' @param times_e the fine grid of `make_prnode_design()`; `NULL` takes
#' `fit$design$times_e` when `fit` is given, else derives it from each
#' pairing's mother's frames.
#' @param x_ref `NULL` (at `q > 1`, every pairing's two fields are linearized
#' at their own cells' mean observed states, the design's `mean_state`) or a
#' list with `M` and `D`, two fixed length-`p` states applied to every
#' pairing (a sensitivity reading).
#' @param verbose a non-negative number; `1` reports each pairing fitted.
#' @returns an object of class `prnode_inheritance_test`: `estimate` (the
#' true pair's), `fit` (its `prnode_fit`), `fits` (`NULL`, or with
#' `bool_keep_fits` the list of wrong-pairing fits), `h`, `null` (the
#' wrong pairings' estimates, named as `null_mothers`), `num_null`,
#' `p_value`.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' # two wrong mothers: the mothers of two other simulated pairs
#' others <- lapply(2:3, function(s){
#'   build_md_dataset(seed = s, num_d = 18, noise_level_ratio = 0.1,
#'                    truth_args = list(p = 4),
#'                    obs_time = seq(0, 1, length.out = 36),
#'                    steps_per_frame = 100)$raw$Y_M
#' })
#' test <- inheritance_test(ds$raw$m_obs_time, ds$raw$Y_M, ds$raw$d_obs_time,
#'                          ds$raw$Y_D, null_mothers = others)
#' print(test)
#' @export
inheritance_test <- function(m_obs_time,
                             Y_M,
                             d_obs_time,
                             Y_D,
                             null_mothers,
                             bool_keep_fits = FALSE,
                             design_args = list(),
                             fit = NULL,
                             fit_args = list(),
                             h = NULL,
                             mu = c(S = 1, M = 1, D = 1),
                             q = 1,
                             times_e = NULL,
                             x_ref = NULL,
                             verbose = 0){
  Y_M <- as.matrix(Y_M)
  Y_D <- as.matrix(Y_D)
  stopifnot(is.numeric(Y_M), is.numeric(Y_D), ncol(Y_M) == ncol(Y_D),
            length(m_obs_time) == nrow(Y_M), length(d_obs_time) == nrow(Y_D),
            length(bool_keep_fits) == 1, is.logical(bool_keep_fits), !is.na(bool_keep_fits),
            length(q) == 1, is.numeric(q), q >= 1, q == round(q),
            is.list(design_args), is.list(fit_args))
  if(length(design_args) > 0 && (is.null(names(design_args)) || any(names(design_args) == ""))){
    stop("`design_args` must be a named list")
  }
  if(length(fit_args) > 0 && (is.null(names(fit_args)) || any(names(fit_args) == ""))){
    stop("`fit_args` must be a named list")
  }
  reserved <- c("m_obs_time", "Y_M", "d_obs_time", "Y_D", "q", "times_e")
  if(any(names(design_args) %in% reserved)){
    stop("`design_args` may not set ", paste(intersect(names(design_args), reserved), collapse = ", "),
         "; those are arguments of inheritance_test()")
  }
  if(any(names(fit_args) %in% c("public", "lambda", "mu", "verbose"))){
    stop("`fit_args` may not set `public`, `lambda`, `mu` or `verbose`")
  }
  # the pairing is fitted from the raw observations, so the smoothing level is
  # chosen per pairing by the rule (design_args$spar fixes it)
  fit_pairing <- function(m_time, Y){
    do.call(fit_prnode, c(list(m_obs_time = m_time, Y_M = Y, d_obs_time = d_obs_time, Y_D = Y_D, q = q,
                               times_e = times_e, mu = mu, verbose = max(verbose - 1, 0)),
                          design_args, fit_args))
  }
  if(length(verbose) != 1 || !is.numeric(verbose) || is.na(verbose) || verbose < 0){
    stop("`verbose` must be a single non-negative number")
  }
  if(!is.list(null_mothers) || length(null_mothers) < 1){
    stop("`null_mothers` must be a non-empty list of wrong mothers (matrices, or lists with `Y_M` and `m_obs_time`)")
  }
  p <- ncol(Y_M)
  null_list <- lapply(seq_along(null_mothers), function(i){
    wm <- null_mothers[[i]]
    if(is.list(wm) && !is.data.frame(wm)){
      if(is.null(wm$Y_M)) stop("`null_mothers[[", i, "]]` is a list without `Y_M`")
      wm_Y <- as.matrix(wm$Y_M)
      wm_time <- if(is.null(wm$m_obs_time)) m_obs_time else wm$m_obs_time
    } else {
      wm_Y <- as.matrix(wm)
      wm_time <- m_obs_time
    }
    if(!is.numeric(wm_Y) || ncol(wm_Y) != p || length(wm_time) != nrow(wm_Y)){
      stop("`null_mothers[[", i, "]]` must be an n x ", p, " matrix on ", length(m_obs_time),
           " frames (or a list with `Y_M` and a matching `m_obs_time`)")
    }
    list(m_obs_time = wm_time, Y_M = wm_Y)
  })
  null_names <- names(null_mothers)
  if(is.null(null_names)) null_names <- paste0("null_", seq_along(null_mothers))
  if(is.null(fit)){
    if(verbose > 0) print("inheritance_test: fitting the true pair")
    fit <- fit_pairing(m_obs_time, Y_M)
  } else {
    if(!inherits(fit, "prnode_fit") || is.null(fit$tuning)) stop("`fit` must be a prnode_fit carrying its tuning")
    if(any(abs(fit$mu - .mu_vec(mu)) > 1e-12)) stop("`fit` was fitted at another `mu`")
    searched <- !is.null(fit$tuning$spar_criterion)
    if(searched == !is.null(design_args$spar)){
      stop("`fit` was ", if(searched) "tuned with the spar search" else paste0("fitted at spar = ", fit$tuning$spar),
           " but `design_args$spar` is ", if(is.null(design_args$spar)) "NULL (the wrong pairings would search)" else
           "given (the wrong pairings would be fixed)", "; the wrong pairings must be fitted by the rule the true pair was")
    }
    if(!searched && !isTRUE(all.equal(fit$tuning$spar, design_args$spar))){
      stop("`fit` was fitted at spar = ", fit$tuning$spar, ", not `design_args$spar` = ", design_args$spar)
    }
    if(fit$design$q != q || fit$design$n_M != nrow(Y_M) || fit$design$n_D != nrow(Y_D) || fit$design$p != p){
      stop("`fit` does not match the data's dimensions (q, n_M, n_D, p)")
    }
    # at q > 1 the true pair is read at the states stored in `fit`'s design
    # and every wrong pairing at the states of the data passed here, so the
    # two must be the same data (the dimensions alone would not catch a
    # rescaled or reordered Y)
    if(q > 1 && is.null(x_ref)){
      stored <- fit$tuning$path$design$mean_state
      if(is.null(stored) || !isTRUE(all.equal(unname(stored$M), unname(colMeans(Y_M)), tolerance = 1e-10)) ||
         !isTRUE(all.equal(unname(stored$D), unname(colMeans(Y_D)), tolerance = 1e-10))){
        stop("`fit`'s design stores other mean states (`mean_state`) than the column means of `Y_M` and `Y_D`: ",
             "at q = ", q, " the true pair and the wrong pairings would be linearized at different states; ",
             "pass the data `fit` was built on")
      }
    }
    # The wrong pairings must integrate on the grid the true pair did, or the
    # tie of a pairing with the true mother is lost in the last digits.
    if(is.null(times_e) && !is.null(fit$design$times_e)) times_e <- fit$design$times_e
  }
  if(is.null(h)) h <- stats::median(diff(m_obs_time))
  estimate <- inheritance_flow(fit, h = h, x_ref = x_ref)

  null_estimates <- numeric(length(null_list))
  fits <- if(bool_keep_fits) vector("list", length(null_list)) else NULL
  for(i in seq_along(null_list)){
    if(verbose > 0) print(paste0("inheritance_test: fitting wrong pairing ", i, " of ", length(null_list)))
    wm <- null_list[[i]]
    fit_i <- fit_pairing(wm$m_obs_time, wm$Y_M)
    null_estimates[i] <- inheritance_flow(fit_i, h = h, x_ref = x_ref)
    if(bool_keep_fits) fits[[i]] <- fit_i
  }
  names(null_estimates) <- null_names
  if(bool_keep_fits) names(fits) <- null_names
  structure(list(estimate = estimate,
                 fit = fit,
                 fits = fits,
                 h = h,
                 null = null_estimates,
                 num_null = length(null_estimates),
                 p_value = inheritance_pvalue(estimate, null_estimates)),
            class = "prnode_inheritance_test")
}

#' Print a prnode_inheritance_test
#'
#' @param x a `prnode_inheritance_test`.
#' @param ... ignored.
#' @returns `x`, invisibly.
#' @export
print.prnode_inheritance_test <- function(x, ...){
  cat("Wrong-mother permutation test of inheritance (class prnode_inheritance_test)\n")
  cat("  estimate (flow-map similarity at the pick):", signif(x$estimate, 3), "\n")
  cat("  wrong pairings:", x$num_null, " null mean", signif(mean(x$null), 3), " null sd",
      if(x$num_null > 1) signif(stats::sd(x$null), 3) else NA, " null max", signif(max(x$null), 3), "\n")
  cat("  p-value:", signif(x$p_value, 3), " (smallest attainable", signif(1/(1 + x$num_null), 3), ")\n")
  cat("  horizon h =", signif(x$h, 3), "\n")
  invisible(x)
}

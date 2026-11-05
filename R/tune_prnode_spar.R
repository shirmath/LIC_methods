#' Choose the smoothing level and the penalty of a PRNODE fit together
#'
#' The smoothing level `spar` of `smooth_md_trajectory()` changes the design
#' (the integrated basis is built on the smoothed trajectories), so it cannot
#' be chosen by `tune_prnode()`, which sees one built design. This function
#' builds one design per level of `spar_grid` from the raw observations, tunes
#' each by `tune_prnode()` and keeps the level whose known-variance BIC is
#' smallest, the minimum over `spar` and over the path jointly. The criteria are comparable across `spar` because the profiled
#' response and the weights do not depend on it (with the default
#' `noise_method = "diff2"`; the `"spline"` noise estimate does depend on
#' `spar`, and then the criteria are comparable only approximately).
#'
#' @param m_obs_time,Y_M the mother's frames and `n_M x p` observations.
#' @param d_obs_time,Y_D the daughter's frames and observations.
#' @param bool_scale_response,noise_method,omega,q,times_e passed to
#' `make_prnode_design()`.
#' @param edge_tol,max_iter,mu,num_decades,num_lambda,tol,verbose passed to
#' `tune_prnode()` (`verbose` one lower).
#' @param spar_grid the smoothing levels to try, each in `[0, 1]`. The
#' default is `c(0, 0.25, 0.5, 0.75)`: level 1 is left out because its
#' near-collinear regressors make the solver slow.
#' Ties in the criterion go to the level that comes first in `spar_grid`.
#'
#' @returns a list: `public` (the chosen design), `spar` (the chosen level),
#' `spar_criterion` (the minimum criterion at each level, named by level) and
#' `tuning` (the chosen design's `prnode_tuning`).
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' out <- tune_prnode_spar(ds$raw$m_obs_time, ds$raw$Y_M, ds$raw$d_obs_time,
#'                         ds$raw$Y_D, spar_grid = c(0.25, 0.5),
#'                         times_e = ds$public$times_e,
#'                         num_lambda = 20, num_decades = 2)
#' out$spar
#' out$spar_criterion
#' @export
tune_prnode_spar <- function(m_obs_time,
                             Y_M,
                             d_obs_time,
                             Y_D,
                             bool_scale_response = TRUE,
                             edge_tol = 1e-8,
                             max_iter = 1e4,
                             mu = c(S = 1, M = 1, D = 1),
                             noise_method = c("diff2", "spline"),
                             num_decades = 4,
                             num_lambda = 80,
                             omega = NULL,
                             q = 1,
                             spar_grid = c(0, 0.25, 0.5, 0.75),
                             times_e = NULL,
                             tol = 1e-10,
                             verbose = 0){
  noise_method <- match.arg(noise_method)
  stopifnot(is.numeric(spar_grid), length(spar_grid) >= 1, !anyNA(spar_grid), all(spar_grid >= 0),
            all(spar_grid <= 1), !anyDuplicated(spar_grid))
  if(length(verbose) != 1 || !is.numeric(verbose) || is.na(verbose) || verbose < 0){
    stop("`verbose` must be a single non-negative number")
  }
  tuned <- lapply(spar_grid, function(spar){
    if(verbose > 0) print(paste0("tune_prnode_spar: spar = ", spar))
    public <- make_prnode_design(m_obs_time = m_obs_time, Y_M = Y_M, d_obs_time = d_obs_time, Y_D = Y_D,
                                 bool_scale_response = bool_scale_response, noise_method = noise_method,
                                 omega = omega, q = q, spar = spar, times_e = times_e)
    tuning <- tune_prnode(public, edge_tol = edge_tol, max_iter = max_iter, mu = mu, num_decades = num_decades,
                          num_lambda = num_lambda, tol = tol, verbose = max(verbose - 1, 0))
    list(public = public, tuning = tuning)
  })
  spar_criterion <- vapply(tuned, function(z) min(z$tuning$criterion), numeric(1))
  names(spar_criterion) <- as.character(spar_grid)
  best <- which.min(spar_criterion)
  list(public = tuned[[best]]$public,
       spar = spar_grid[best],
       spar_criterion = spar_criterion,
       tuning = tuned[[best]]$tuning)
}

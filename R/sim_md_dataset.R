#' Noise statistic of one observed series
#'
#' The standard deviation of the series' second differences divided by
#' `sqrt(6)` (the standard deviation of a second difference of iid noise with
#' unit variance), relative to the series' own standard deviation. A rough
#' noise-to-signal ratio that needs no smoother.  A constant series gives
#' `NaN` (0/0), which is left as is: it cannot happen for a
#' stable ODE with a nonzero intercept or initial state.
#'
#' @param y_vec numeric vector, one series in time order; at least 3 values.
#' @returns a non-negative scalar, or `NaN` for a constant series.
#' @export
noise_stat <- function(y_vec){
  stopifnot(is.numeric(y_vec), length(y_vec) >= 3)
  stats::sd(diff(y_vec, differences = 2))/sqrt(6)/stats::sd(y_vec)
}

#' Draw one accepted mother-daughter truth
#'
#' Draws a network pair with `sim_md_network()` and solves it with
#' `sim_md_trajectory()` without observation noise, with a fixed birth time (the
#' last `num_d` mother frames belong to the daughter), redrawing until the
#' daughter's trajectory stays bounded: the largest absolute value of her
#' noiseless observations is at most `max_growth` times the mother's. This
#' rejects a daughter whose ODE blows up within the movie (a deviation can
#' make her system unstable, and an exploding trajectory resembles no protein
#' intensity and is no test of an estimator). With the antisymmetric pairs of
#' `sim_md_network()` at `q = 1` the daughter's system is always stable
#' (`A + A^T` is twice the negative diagonal), so the check does not reject
#' default draws; it is kept for `q > 1`, where stability is not guaranteed.
#' With `deviation_scale` other than 1 the daughter's drawn deviation `d_delta`
#' is multiplied by it before her trajectory is solved (her coefficients become
#' `s_pars + deviation_scale * d_delta`; her adjacency is unchanged for a
#' nonzero scale, and at 0 her private edges are gone and the edge types are
#' relabelled), so the growth bound judges the daughter that is actually
#' simulated. The seed is set
#' once before the loop and the attempts consume the stream in sequence, so the
#' accepted draw is a function of `seed_number` alone and the stream is left
#' positioned after it for whatever the caller draws next.
#'
#' @param num_d number of daughter frames; the daughter is born at
#' `obs_time[length(obs_time) + 1 - num_d]`, so it must be between 1 and
#' `length(obs_time)`.
#' @param obs_time the mother's observation times, increasing.
#' @param truth_args a list of arguments for `sim_md_network()`; `p` is
#' required, and `q = 1`, `num_shared = p`, `num_mother_only = num_daughter_only
#' = min(2, (p(p - 1) - p) %/% 2)` (2 for `p >= 4`) and `frac_unequal = 0.5`
#' are the defaults for the others.
#' @param deviation_scale non-negative multiplier on the daughter's drawn
#' deviation; `1` keeps the draw as it is.
#' @param max_attempts give up after this many rejected draws.
#' @param max_growth the bound on the daughter's size relative to the
#' mother's (ratio of their largest absolute noiseless values); `Inf` accepts
#' the first draw.
#' @param steps_per_frame passed to `sim_md_trajectory()` when `times_e` is
#' `NULL`.
#' @param times_e the fine grid of `sim_md_trajectory()`; `NULL` derives it
#' from the frames.
#' @param seed_number seed set before the first attempt; `NULL` leaves the
#' stream untouched.
#' @returns a list with `num_attempts` (the accepted attempt, counted from 1)
#' and `sim`, the `sim_md_trajectory()` output of the accepted draw (`sigma =
#' 0`, so `Y == Y_noiseless`; `sim$truth` is the network draw).
#' @export
draw_md_truth <- function(num_d,
                          obs_time,
                          truth_args,
                          deviation_scale = 1,
                          max_attempts = 500,
                          max_growth = 3,
                          steps_per_frame = 2000,
                          times_e = NULL,
                          seed_number = 10){
  stopifnot(length(num_d) == 1, num_d >= 1, num_d <= length(obs_time),
            max_attempts >= 1, max_growth >= 0, is.list(truth_args))
  .check_deviation_scale(deviation_scale)
  truth_args <- .fill_truth_args(truth_args)
  if(is.null(times_e)) times_e <- .derive_grid(obs_time, steps_per_frame)
  if(!is.null(seed_number)) set.seed(seed_number)

  birth_time <- obs_time[length(obs_time) + 1 - num_d]
  for(attempt in seq_len(max_attempts)){
    truth <- do.call(sim_md_network, c(truth_args, list(seed_number = NULL)))
    # the scaling draws nothing, so the stream contract above is unchanged
    if(deviation_scale != 1) truth <- .scale_deviation(truth, deviation_scale)
    sim <- sim_md_trajectory(truth = truth, obs_time = obs_time, birth_time = birth_time,
                             sigma = rep(0, truth$p), times_e = times_e, seed_number = NULL)
    bool_accept <- max(abs(sim$daughter$Y_noiseless)) <=
      max_growth*max(abs(sim$mother$Y_noiseless))
    if(isTRUE(bool_accept)) break
    if(attempt == max_attempts){
      stop("no draw kept the daughter within `max_growth` = ", max_growth, " times the mother's largest ",
           "absolute value in `max_attempts` = ", max_attempts, " attempts")
    }
  }

  stopifnot(nrow(sim$daughter$Y_noiseless) == num_d)
  list(num_attempts = attempt,
       sim = sim)
}

#' Multiply a truth's daughter deviation by a scale
#'
#' The daughter's coefficients become `s_pars + scale * d_delta`; her adjacency
#' and the scaffold's follow, and at `scale = 0` the edge types are relabelled
#' (her private edges are absent; a shared-unequal edge whose mother deviation
#' is also zero is shared-equal). The one implementation behind
#' `deviation_scale` in `draw_md_truth()` and `build_md_dataset()`.
#' @param truth a `sim_md_network()` list whose `d_delta` is the unscaled draw.
#' @returns the truth with `d_delta`, `d_ode_pars`, `d_network`, `s_network`
#' and (at scale 0) `edge_type` replaced.
#' @noRd
.scale_deviation <- function(truth,
                             scale){
  truth$d_delta <- scale*truth$d_delta
  truth$d_ode_pars <- truth$s_pars + truth$d_delta
  truth$d_network <- .network_of(truth$d_ode_pars, truth$q)
  truth$s_network <- truth$m_network & truth$d_network
  if(scale == 0){
    truth$edge_type[truth$edge_type == "daughter_only"] <- "absent"
    first_cols <- seq_len(ncol(truth$edge_type))*truth$q - truth$q + 1
    truth$edge_type[truth$edge_type == "shared_unequal" & truth$m_delta[, first_cols] == 0] <- "shared_equal"
  }
  truth
}

#' Guard of the `deviation_scale` argument
#' @noRd
.check_deviation_scale <- function(deviation_scale){
  if(length(deviation_scale) != 1 || !is.finite(deviation_scale) || deviation_scale < 0){
    stop("`deviation_scale` must be one finite non-negative number, not ", paste(deviation_scale, collapse = ", "))
  }
  invisible(TRUE)
}

#' Defaults of the truth arguments of `draw_md_truth()`
#' @noRd
.fill_truth_args <- function(truth_args){
  if(is.null(truth_args$p)) stop("`truth_args` must give the number of proteins `p`")
  p <- truth_args$p
  # two private edges per cell when the p(p - 1) positions allow it after p shared ones
  num_private <- max(0, min(2, (p*(p - 1) - p) %/% 2))
  defaults <- list(q = 1, num_shared = p, num_mother_only = num_private, num_daughter_only = num_private,
                   frac_unequal = 0.5)
  for(nm in names(defaults)) if(is.null(truth_args[[nm]])) truth_args[[nm]] <- defaults[[nm]]
  truth_args
}

#' Draw per-protein noise ratios around a level
#'
#' Each protein's ratio is `noise_level_ratio` times a log-uniform jitter within
#' a factor of `noise_jitter`. One `stats::runif()` call of length `p`.
#'
#' @param p number of proteins.
#' @param noise_level_ratio the level (noise sd over the series' own sd).
#' @param noise_jitter jitter factor; `1` gives every protein the level exactly.
#' @param seed_number seed set before the draw; `NULL` leaves the stream untouched.
#' @returns a numeric vector of length `p`.
#' @export
draw_noise_ratio <- function(p,
                             noise_level_ratio,
                             noise_jitter = 1.5,
                             seed_number = 10){
  stopifnot(p >= 1, noise_level_ratio >= 0, noise_jitter >= 1)
  if(!is.null(seed_number)) set.seed(seed_number)
  noise_level_ratio*exp(stats::runif(p, -log(noise_jitter), log(noise_jitter)))
}

#' Add observation noise to a noiseless mother-daughter pair
#'
#' Protein `j` of each cell gets Gaussian noise with standard deviation
#' `sigma_j = noise_ratio_vec[j]` times the **mother's** noiseless standard
#' deviation of protein `j`, the same absolute noise in both cells (one
#' microscope). With `bool_per_cell_sd = TRUE` each cell is scaled
#' by its own noiseless sd instead, the convention of the triplet banks. The
#' mother's noise matrix is drawn first (one `stats::rnorm()` call filling the
#' matrix column-major), then the daughter's. With `noise = "ar1"` each
#' protein's errors are an AR(1) series with coefficient `rho` and stationary
#' variance `sigma_j^2`.
#'
#' @param m_noiseless_mat,d_noiseless_mat noiseless observations, rows are
#' frames and columns are proteins; the same number of columns.
#' @param noise_ratio_vec one ratio per protein (recycled if a scalar).
#' @param bool_per_cell_sd scale each cell's noise by its own sd (`TRUE`) or
#' both by the mother's (`FALSE`).
#' @param noise `"iid"` or `"ar1"`.
#' @param rho the AR(1) coefficient.
#' @param seed_number seed set before the draws; `NULL` leaves the stream untouched.
#' @returns a list with `sigma_D`, `sigma_M` (the absolute sds used) and `Y_D`,
#' `Y_M` (the noisy matrices).
#' @export
add_md_noise <- function(m_noiseless_mat,
                         d_noiseless_mat,
                         noise_ratio_vec,
                         bool_per_cell_sd = FALSE,
                         noise = c("iid", "ar1"),
                         rho = 0,
                         seed_number = 10){
  noise <- match.arg(noise)
  p <- ncol(m_noiseless_mat)
  stopifnot(is.matrix(m_noiseless_mat), is.matrix(d_noiseless_mat),
            ncol(d_noiseless_mat) == p, all(noise_ratio_vec >= 0),
            length(noise_ratio_vec) %in% c(1, p), abs(rho) < 1)
  noise_ratio_vec <- rep_len(noise_ratio_vec, p)
  if(!is.null(seed_number)) set.seed(seed_number)

  sigma_M <- noise_ratio_vec*apply(m_noiseless_mat, 2, stats::sd)
  sigma_D <- if(bool_per_cell_sd) noise_ratio_vec*apply(d_noiseless_mat, 2, stats::sd) else sigma_M
  Y_M <- m_noiseless_mat + .draw_errors(nrow(m_noiseless_mat), p, sigma_M, noise, rho)
  Y_D <- d_noiseless_mat + .draw_errors(nrow(d_noiseless_mat), p, sigma_D, noise, rho)
  list(sigma_D = sigma_D,
       sigma_M = sigma_M,
       Y_D = Y_D,
       Y_M = Y_M)
}

#' Smooth observed trajectories onto a fine grid
#'
#' A smoothing spline with a fixed `spar`: one `stats::smooth.spline()` per protein with a
#' knot at every observation and the fixed smoothing parameter `spar`, predicted
#' on `times_e`. `times_e` may extend beyond the observed times; the spline is
#' then extrapolated linearly by `stats::predict()`.
#'
#' @param obs_time observation times, increasing, length `n`.
#' @param obs_mat `n x p` matrix of observations, columns are proteins.
#' @param times_e the grid to predict on.
#' @param spar the smoothing parameter of `stats::smooth.spline()`; `NULL` lets
#' the spline choose it by its own generalized cross-validation (not what the
#' tuning rule does: `tune_prnode_spar()` chooses one fixed level per pair
#' from a grid, jointly with the penalty).
#' @returns a `length(times_e) x p` matrix of smoothed trajectories.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' X_smooth <- smooth_md_trajectory(ds$raw$m_obs_time, ds$raw$Y_M,
#'                                  ds$raw$times_e)
#' dim(X_smooth)
#' # the residual at the frames relative to each protein's spread; the fixed
#' # spar = 0.5 is a strong smoother on series this fast (frames sit on the
#' # grid to 1e-10, not bit for bit, so look them up by nearest point)
#' idx <- sapply(ds$raw$m_obs_time,
#'               function(tt) which.min(abs(ds$raw$times_e - tt)))
#' round(apply(X_smooth[idx, ] - ds$raw$Y_M, 2, stats::sd)/
#'         apply(ds$raw$Y_M, 2, stats::sd), 2)
#' @export
smooth_md_trajectory <- function(obs_time,
                                 obs_mat,
                                 times_e,
                                 spar = 0.5){
  stopifnot(is.matrix(obs_mat), nrow(obs_mat) == length(obs_time),
            !is.unsorted(obs_time), length(obs_time) >= 4)
  p <- ncol(obs_mat)

  smooth_mat <- matrix(NA_real_, nrow = length(times_e), ncol = p)
  for(j in seq_len(p)){
    fit <- stats::smooth.spline(obs_time, obs_mat[, j], all.knots = TRUE, spar = spar)
    smooth_mat[, j] <- stats::predict(fit, times_e)$y
  }
  smooth_mat
}

#' The daughter's fine time grid
#'
#' The part of the mother's fine grid from the daughter's birth time on, built as
#' `seq(times_e[idx], max, by = step)` from the grid point nearest the birth time,
#' as every pipeline script did (`seq(times_e[map_to(birth_time, times_e)], 1,
#' by = 1e-5)`). The integration step is read from the first two entries of
#' this grid, so it must be evenly spaced.
#'
#' @param birth_time the daughter's birth time.
#' @param times_e the mother's fine grid, evenly spaced and increasing.
#' @returns a numeric vector starting at the grid point nearest `birth_time`
#' and ending at `times_e[length(times_e)]`.
#' @export
make_daughter_grid <- function(birth_time,
                               times_e){
  stopifnot(length(birth_time) == 1, length(times_e) >= 2, !is.unsorted(times_e))
  idx <- map_to(birth_time, times_e)
  seq(times_e[idx], times_e[length(times_e)], by = times_e[2] - times_e[1])
}

#' Simulate one dataset: truth, noisy observations and the estimator's design
#'
#' The dataset pipeline of the PRNODE simulations: one
#' truth from `draw_md_truth()` (seeded with `seed`), per-protein noise ratios
#' from `draw_noise_ratio()`, noise from `add_md_noise()` (one absolute sd per
#' protein for both cells), all on one stream, then the two-block design from
#' `make_prnode_design()`. The `level` decides what the design is built from:
#' `"ode"` smooths the noisy observations (the full pipeline), `"regression"`
#' hands `make_prnode_design()` the true trajectories, so the linear model
#' holds up to quadrature error and the gap between the levels is the effect
#' of smoothing. `meta`, `raw` and `truth` do not depend
#' on the level.
#'
#' Two settings for the daughter. With `deviation_scale = 1` the edge counts
#' of `truth_args` are used as drawn. With another `deviation_scale` the
#' daughter's drawn deviation is multiplied by it before her trajectory is
#' solved, through `draw_md_truth()` (so the growth bound judges the scaled
#' daughter; the adjacency is unchanged for a nonzero scale, and `0` gives the
#' daughter the scaffold's coefficients with the edge types relabelled); this is one way to set the pair's inheritance level. The other is
#' the editor `rewire_network()`, whose result is handed in as `truth`.
#'
#' A truth drawn by hand with `sim_md_network()` (or assembled by
#' `md_truth_from_pars()`) can be passed as `truth`
#' instead of being drawn here (so that one truth runs through a whole
#' document or study). Then `truth_args` must be empty, `seed` governs the
#' noise alone (the stream is seeded once before the noise ratios are drawn),
#' the growth bound is checked once and errors instead of redrawing, and
#' `meta$n_attempts` is `NA`. `deviation_scale` still applies: it rescales the
#' given daughter deviation.
#'
#' `truth$inheritance_std` is the same similarity on the two systems in the
#' standardized units the estimator reads (`S^{-1} A S` with `S` the diagonal
#' of the pooled sd of the two cells' noiseless states), the target of
#' `inheritance_flow()`.
#'
#' `truth$inheritance` is the pair's true inheritance, the package's measure
#' (`flow_map_similarity()`) on the two true systems: the one-frame flow-map
#' similarity of `truth$m_ode_pars` and `truth$d_ode_pars` at the horizon of
#' one sampling interval (the median spacing of `obs_time`); at a truth degree
#' above 1 each field is linearized at the mean of its own cell's noiseless
#' frames (the estimator `inheritance_flow()` linearizes each fitted field at
#' the mean of that cell's observed frames, so at that degree the reference
#' states differ by the noise in those means; the truth's number is
#' noise-free by design). At degree 1 it is a function of the two coefficient
#' matrices alone, so it does not depend on the noise, on `num_d` or on the
#' trajectories, and the same truth reads the same number in every dataset
#' built from it (above degree 1 the daughter's mean state moves with
#' `num_d`). It is 1 if and
#' only if the daughter's system is the mother's (`d_delta` equal to
#' `m_delta`); `deviation_scale = 0` gives the daughter the scaffold, which is
#' a clone of the mother only when the mother has no deviation of her own.
#' Two unrelated systems read a chance level near 0.3, not 0 (see
#' `flow_map_similarity()`), and scaling the daughter's deviation up moves the
#' similarity toward that level, not monotonically, so the similarity is read
#' off a dataset rather than targeted by a search.
#'
#' @param seed the dataset's seed.
#' @param num_d number of daughter frames, between 4 (the noise-sd estimate of
#' the loss weights needs four frames) and `length(obs_time)`; the daughter is
#' born at `obs_time[length(obs_time) + 1 - num_d]`.
#' @param noise_level_ratio the noise level (noise sd over the mother's own sd).
#' @param truth_args a list of arguments for `sim_md_network()`; `p` is
#' required (see `draw_md_truth()` for the defaults).
#' @param deviation_scale non-negative multiplier on the daughter's deviation;
#' `1` leaves it as drawn or given.
#' @param level `"ode"` or `"regression"`.
#' @param max_attempts,max_growth passed to `draw_md_truth()`.
#' @param noise,rho the noise model of `add_md_noise()`.
#' @param noise_jitter passed to `draw_noise_ratio()`.
#' @param noise_level an optional label for the level, stored in `meta`.
#' @param obs_time the mother's observation times.
#' @param q degree of the polynomial basis of the **design** (the truth's
#' degree is `truth_args$q`, 1 by default).
#' @param spar passed to `make_prnode_design()`.
#' @param steps_per_frame grid steps per frame interval when `times_e` is `NULL`.
#' @param times_e the fine grid; `NULL` derives it from `obs_time`.
#' @param truth `NULL` (the truth is drawn by `draw_md_truth()`), or a
#' `sim_md_network()` list to use as the truth.
#' @returns a list with `meta` (`dataset_id`, `seed`, `num_d`, `n_M`, `n_D`,
#' `q`, `level`, `noise_level`, `noise_level_ratio`, `noise_ratio`,
#' `noise_stat_m`, `noise_stat_d`, `n_attempts` (`NA` for a given `truth`),
#' `deviation_scale`, `birth_time`, `length_ratio`), `raw` (`m_obs_time`, `Y_M`, `d_obs_time`,
#' `Y_D`, `birth_time`, `times_e`, `protein_names`: what an estimator reads),
#' `truth` (`s_pars`, `m_delta`, `d_delta`, `m_ode_pars`, `d_ode_pars`,
#' `m_intercept`, `d_intercept`, `s_network`, `m_network`, `d_network`,
#' `edge_type`, `X_M`, `X_D`, `Y_M_noiseless`, `Y_D_noiseless`, `sigma_M`,
#' `sigma_D`, `gamma_M`, `gamma_D`, `inheritance`, `inheritance_std`, `deviation_scale`) and
#' `public` (the `make_prnode_design()` list). `public` is a reference design
#' at the default weighting and `spar = 0.5`; the default fit path does not use
#' it, because `fit_prnode()` rebuilds the design from `raw` and chooses the
#' smoothing level with `tune_prnode_spar()`. Fitting `public` directly tunes
#' only `lambda` at `spar = 0.5`, which is usually a worse fit (use
#' `fit_prnode(m_obs_time = , Y_M = , d_obs_time = , Y_D = , times_e = )` on
#' `raw` instead). It stays for the existing callers and is a candidate for
#' removal.
#' @examples
#' # four proteins, 36 mother frames, the daughter born at frame 19; a coarse
#' # grid (100 steps per frame interval) keeps the example fast
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' str(ds$meta)
#' ds$truth$edge_type
#' dim(ds$raw$Y_M)
#' dim(ds$raw$Y_D)
#' names(ds$public)
#' # the pair's true inheritance, the flow-map similarity of the two systems
#' ds$truth$inheritance
#' # the same draw with the daughter's deviation doubled: a less inherited daughter
#' ds_scaled <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                               truth_args = list(p = 4),
#'                               obs_time = seq(0, 1, length.out = 36),
#'                               steps_per_frame = 100, deviation_scale = 2)
#' c(as_drawn = ds$truth$inheritance, doubled = ds_scaled$truth$inheritance)
#' # a truth drawn by hand, carried through unchanged
#' truth <- sim_md_network(p = 4, num_shared = 4, num_mother_only = 1,
#'                         num_daughter_only = 1, seed_number = 3)
#' ds_given <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                              obs_time = seq(0, 1, length.out = 36),
#'                              steps_per_frame = 100, truth = truth)
#' identical(ds_given$truth$m_ode_pars, truth$m_ode_pars)
#' @export
build_md_dataset <- function(seed,
                             num_d,
                             noise_level_ratio,
                             truth_args = list(),
                             deviation_scale = 1,
                             level = c("ode", "regression"),
                             max_attempts = 500,
                             max_growth = 3,
                             noise = c("iid", "ar1"),
                             noise_jitter = 1.5,
                             noise_level = NA,
                             obs_time = seq(0, 1, length.out = 48),
                             q = 1,
                             rho = 0,
                             spar = 0.5,
                             steps_per_frame = 2000,
                             times_e = NULL,
                             truth = NULL){
  level <- match.arg(level)
  noise <- match.arg(noise)
  n_M <- length(obs_time)
  if(!is.null(truth)){
    if(!is.list(truth) || is.null(truth$p) || !is.matrix(truth$m_ode_pars) || is.null(truth$edge_type)){
      stop("`truth` must be a `sim_md_network()` list (with `p`, `m_ode_pars`, `edge_type`, ...)")
    }
    if(length(truth_args) > 0){
      stop("`truth` is given, so `truth_args` must be empty: the networks are not drawn here")
    }
  }
  # 4, not the design's 3: the loss weights are estimated from the data here
  # (`omega = NULL`), which needs four frames per cell (test plan H8)
  if(length(num_d) != 1 || num_d < 4 || num_d > n_M){
    stop("`num_d` = ", num_d, " must lie between 4 and the number of mother frames (", n_M, ")")
  }
  stopifnot(length(seed) == 1, noise_level_ratio >= 0, is.list(truth_args))
  .check_deviation_scale(deviation_scale)
  if(is.null(times_e)) times_e <- .derive_grid(obs_time, steps_per_frame)

  if(is.null(truth)){
    truth_args <- .fill_truth_args(truth_args)
    truth_draw <- draw_md_truth(num_d = num_d, obs_time = obs_time, truth_args = truth_args,
                                deviation_scale = deviation_scale, max_attempts = max_attempts,
                                max_growth = max_growth, times_e = times_e, seed_number = seed)
    sim <- truth_draw$sim
    num_attempts <- truth_draw$num_attempts
  } else {
    # the given truth is solved once, with the stream seeded as draw_md_truth()
    # would have seeded it, so the initial state and the noise draws below
    # depend on `seed` alone
    set.seed(seed)
    if(deviation_scale != 1) truth <- .scale_deviation(truth, deviation_scale)
    sim <- sim_md_trajectory(truth = truth, obs_time = obs_time,
                             birth_time = obs_time[n_M + 1 - num_d], sigma = rep(0, truth$p),
                             times_e = times_e, seed_number = NULL)
    if(!(max(abs(sim$daughter$Y_noiseless)) <= max_growth*max(abs(sim$mother$Y_noiseless)))){
      stop("the given `truth`'s daughter exceeds `max_growth` = ", max_growth, " times the mother's ",
           "largest absolute value, and a given truth is not redrawn")
    }
    num_attempts <- NA_integer_
  }
  truth <- sim$truth
  p <- truth$p
  birth_time <- sim$daughter$birth_time
  d_obs_time <- sim$daughter$obs_time
  X_M <- sim$mother$X
  X_D <- sim$daughter$X
  Y_M_noiseless <- sim$mother$Y_noiseless
  Y_D_noiseless <- sim$daughter$Y_noiseless

  # the ratio continues the stream the truth left; one absolute sigma_j for
  # both cells, read from the mother's noiseless sd
  noise_ratio_vec <- draw_noise_ratio(p, noise_level_ratio, noise_jitter = noise_jitter,
                                      seed_number = NULL)
  sigma_vec <- noise_ratio_vec*apply(Y_M_noiseless, 2, stats::sd)

  noisy <- add_md_noise(Y_M_noiseless, Y_D_noiseless, noise_ratio_vec, bool_per_cell_sd = FALSE,
                        noise = noise, rho = rho, seed_number = NULL)
  stopifnot(isTRUE(all.equal(noisy$sigma_M, sigma_vec)))
  protein_names <- paste0("protein_", seq_len(p))
  Y_M <- noisy$Y_M
  Y_D <- noisy$Y_D
  colnames(Y_M) <- colnames(Y_D) <- protein_names

  # the pair's true inheritance: the package's measure on the two true
  # systems, at the horizon of one sampling interval (the convention of
  # inheritance_flow()); at degree 1 a function of the coefficient matrices
  # alone, above it each field linearized at its own cell's mean noiseless
  # state (so the daughter's length enters through her mean)
  inheritance <- flow_map_similarity(truth$m_ode_pars, truth$d_ode_pars,
                                     h = stats::median(diff(obs_time)), q = truth$q,
                                     x_ref = if(truth$q > 1) list(M = colMeans(Y_M_noiseless),
                                                                  D = colMeans(Y_D_noiseless)) else NULL)

  # the same similarity on the systems in the standardized units the
  # estimator reads (the pooled sd of the two cells' noiseless states)
  inheritance_std <- .flow_map_similarity_std(truth$m_ode_pars, truth$d_ode_pars, Y_M_noiseless,
                                              Y_D_noiseless, h = stats::median(diff(obs_time)), q = truth$q)

  x_list <- if(level == "regression") list(M = X_M, D = X_D) else NULL
  public <- make_prnode_design(m_obs_time = obs_time, Y_M = Y_M, d_obs_time = d_obs_time, Y_D = Y_D,
                               q = q, spar = spar, times_e = times_e, x_list = x_list)

  meta <- list(dataset_id = paste0("seed", seed, "_nD", num_d, "_noise", noise_level_ratio),
               seed = seed,
               num_d = num_d,
               n_M = n_M,
               n_D = num_d,
               q = q,
               level = level,
               noise_level = noise_level,
               noise_level_ratio = noise_level_ratio,
               noise_ratio = noise_ratio_vec,
               noise_stat_m = stats::median(apply(Y_M, 2, noise_stat)),
               noise_stat_d = stats::median(apply(Y_D, 2, noise_stat)),
               n_attempts = num_attempts,
               deviation_scale = deviation_scale,
               birth_time = birth_time,
               length_ratio = num_d/n_M)
  raw <- list(m_obs_time = obs_time,
              Y_M = Y_M,
              d_obs_time = d_obs_time,
              Y_D = Y_D,
              birth_time = birth_time,
              times_e = times_e,
              protein_names = protein_names)
  truth_out <- list(s_pars = truth$s_pars,
                    m_delta = truth$m_delta,
                    d_delta = truth$d_delta,
                    m_ode_pars = truth$m_ode_pars,
                    d_ode_pars = truth$d_ode_pars,
                    m_intercept = truth$m_intercept,
                    d_intercept = truth$d_intercept,
                    s_network = truth$s_network,
                    m_network = truth$m_network,
                    d_network = truth$d_network,
                    edge_type = truth$edge_type,
                    X_M = X_M,
                    X_D = X_D,
                    Y_M_noiseless = Y_M_noiseless,
                    Y_D_noiseless = Y_D_noiseless,
                    sigma_M = noisy$sigma_M,
                    sigma_D = noisy$sigma_D,
                    gamma_M = rbind(Y_M_noiseless[1, ], truth$m_intercept),
                    gamma_D = rbind(Y_D_noiseless[1, ], truth$d_intercept),
                    inheritance = inheritance,
                    inheritance_std = inheritance_std,
                    deviation_scale = deviation_scale)
  list(meta = meta,
       raw = raw,
       truth = truth_out,
       public = public)
}

#' Adjacency of a coefficient matrix with q columns per regulator
#' @noRd
.network_of <- function(pars,
                        q){
  p <- nrow(pars)
  mat <- matrix(FALSE, p, p)
  for(k in seq_len(p)) mat[, k] <- apply(pars[, (k - 1)*q + seq_len(q), drop = FALSE] != 0, 1, any)
  mat
}

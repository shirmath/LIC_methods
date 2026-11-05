#' Simulate a mother-daughter pair of trajectories from a drawn truth
#'
#' The ODE level of the PRNODE simulator. The mother's ODE (`truth$m_ode_pars`, `truth$m_intercept`) is
#' solved exactly by `sim_ode_trajectory()` on the fine grid `times_e` from
#' `init_vals` and observed at `obs_time`. The daughter is born at
#' `birth_time`, one of the mother's frames; she starts from the mother's
#' **noiseless** state at that grid point (never from the noisy observation), is solved with her own ODE on the part of the grid from her
#' birth on, and is observed at the mother's frames from her birth on. Frames
#' are read off the grid by exact index: an observation time more than `1e-10`
#' from every grid point is an error naming that time.
#'
#' Noise is Gaussian with a per-protein standard deviation. One `sigma` vector
#' gives both cells the same absolute noise (one microscope); a list
#' `list(M = , D = )` gives each cell its own. With `noise = "ar1"` each
#' protein's errors follow an AR(1) process with coefficient `rho` and
#' stationary variance `sigma^2`. The random stream is consumed in the order:
#' initial values (when `init_vals` is `NULL`), the mother's errors (one
#' `rnorm()` call, column-major), the daughter's errors. The errors are drawn
#' even when `sigma` is zero, so that a dataset built at two noise levels
#' shares every other draw.
#'
#' @param truth a `sim_md_network()` list (`m_ode_pars`, `d_ode_pars`,
#' `m_intercept`, `d_intercept`; `p` and `q` are read from it or derived from
#' the dimensions of `m_ode_pars`).
#' @param obs_time the mother's observation times, increasing, each on the
#' grid to `1e-10`.
#' @param birth_time the daughter's birth time, one of `obs_time`.
#' @param sigma noise standard deviations: a length-`p` vector used for both
#' cells, or a list with elements `M` and `D`, each length `p`.
#' @param init_vals the mother's state at `times_e[1]`; `NULL` draws
#' `rnorm(p)`.
#' @param noise `"iid"` or `"ar1"`.
#' @param rho the AR(1) coefficient when `noise = "ar1"`, in `(-1, 1)`.
#' @param steps_per_frame when `times_e` is `NULL`, the grid is derived from
#' the frames as `seq(obs_time[1], obs_time[n], length.out = (n - 1) *
#' steps_per_frame + 1)`, so that equally spaced frames lie on the grid to
#' within `1e-10` (not bit for bit: index the grid with a tolerance, never
#' with `match()`).
#' @param times_e the fine grid, increasing and uniform; `NULL` derives it
#' from the frames.
#' @param seed_number seed set before the first draw; `NULL` leaves the stream
#' to the caller.
#' @returns a list with `daughter` (a list with `birth_idx`, the index of the
#' birth frame in `obs_time`; `birth_time`; `obs_time`, her frames; `sigma`;
#' `X`, her trajectory on the grid from her birth on; `Y`, the noisy
#' observations; `Y_noiseless`), `mother` (a list with `obs_time`, `sigma`,
#' `X`, `Y`, `Y_noiseless`) and `truth`, the input truth with `p` and `q`
#' filled in.
#' @examples
#' truth <- sim_md_network(p = 4, num_shared = 4, num_mother_only = 1,
#'                         num_daughter_only = 1, seed_number = 7)
#' obs_time <- seq(0, 1, length.out = 36)
#' sim <- sim_md_trajectory(truth, obs_time = obs_time,
#'                          birth_time = obs_time[19], sigma = rep(0.05, 4),
#'                          steps_per_frame = 100, seed_number = 3)
#' dim(sim$mother$Y)
#' dim(sim$daughter$Y)
#' sim$daughter$birth_time
#' # the daughter starts from the mother's noiseless state at birth
#' all.equal(sim$daughter$Y_noiseless[1, ],
#'           sim$mother$Y_noiseless[sim$daughter$birth_idx, ])
#' @export
sim_md_trajectory <- function(truth,
                              obs_time,
                              birth_time,
                              sigma,
                              init_vals = NULL,
                              noise = c("iid", "ar1"),
                              rho = 0,
                              steps_per_frame = 2000,
                              times_e = NULL,
                              seed_number = NULL){
  noise <- match.arg(noise)
  stopifnot(is.list(truth), is.matrix(truth$m_ode_pars), is.matrix(truth$d_ode_pars),
            !is.unsorted(obs_time, strictly = TRUE), length(birth_time) == 1,
            abs(rho) < 1, steps_per_frame >= 1)
  p <- nrow(truth$m_ode_pars)
  if(is.null(truth$p)) truth$p <- p
  if(is.null(truth$q)) truth$q <- ncol(truth$m_ode_pars)/p
  q <- truth$q
  stopifnot(truth$p == p, ncol(truth$m_ode_pars) == p*q, length(truth$m_intercept) == p,
            length(truth$d_intercept) == p)
  n_M <- length(obs_time)
  if(is.null(times_e)){
    times_e <- .derive_grid(obs_time, steps_per_frame)
  }
  stopifnot(length(times_e) >= 2, !is.unsorted(times_e, strictly = TRUE))
  m_idx <- .frame_index(obs_time, times_e, "obs_time")
  birth_idx <- which(abs(obs_time - birth_time) < 1e-10)
  if(length(birth_idx) != 1){
    stop("`birth_time` = ", birth_time, " must be one of `obs_time`")
  }
  sigma_list <- .sigma_list(sigma, p)

  if(!is.null(seed_number)) set.seed(seed_number)
  if(is.null(init_vals)) init_vals <- stats::rnorm(p)
  stopifnot(length(init_vals) == p)

  X_M <- sim_ode_trajectory(truth$m_ode_pars, truth$m_intercept, init_vals, times_e, q = q)
  Y_M_noiseless <- X_M[m_idx, , drop = FALSE]
  d_grid <- make_daughter_grid(birth_time, times_e)
  d_obs_time <- obs_time[birth_idx:n_M]
  X_D <- .sim_daughter(truth$d_ode_pars, truth$d_intercept, X_M, birth_time, times_e)
  d_idx <- .frame_index(d_obs_time, d_grid, "obs_time")
  Y_D_noiseless <- X_D[d_idx, , drop = FALSE]

  # the mother's errors are one rnorm() call, then the daughter's; both are
  # drawn even at sigma = 0 so that the stream position does not depend on
  # the noise level
  Y_M <- Y_M_noiseless + .draw_errors(n_M, p, sigma_list$M, noise, rho)
  Y_D <- Y_D_noiseless + .draw_errors(length(d_obs_time), p, sigma_list$D, noise, rho)

  stopifnot(identical(X_D[1, ], X_M[m_idx[birth_idx], ]))
  list(daughter = list(birth_idx = birth_idx,
                       birth_time = obs_time[birth_idx],
                       obs_time = d_obs_time,
                       sigma = sigma_list$D,
                       X = X_D,
                       Y = Y_D,
                       Y_noiseless = Y_D_noiseless),
       mother = list(obs_time = obs_time,
                     sigma = sigma_list$M,
                     X = X_M,
                     Y = Y_M,
                     Y_noiseless = Y_M_noiseless),
       truth = truth)
}

#' The daughter's trajectory from the mother's state at birth
#'
#' The one solver call behind the daughter of `sim_md_trajectory()`, so that
#' a rescaled deviation (`deviation_scale`) is solved exactly as the generator
#' solves it.
#'
#' @param d_ode_pars,d_intercept the daughter's ODE.
#' @param X_M the mother's trajectory on `times_e`.
#' @param birth_time the birth time (a grid point).
#' @param times_e the mother's grid.
#' @returns the daughter's trajectory on `times_e[times_e >= birth_time]`.
#' @noRd
.sim_daughter <- function(d_ode_pars,
                          d_intercept,
                          X_M,
                          birth_time,
                          times_e){
  birth_grid_idx <- .frame_index(birth_time, times_e, "birth_time")
  # the one construction of the daughter's grid, shared with make_prnode_design()
  d_grid <- make_daughter_grid(birth_time, times_e)
  q <- ncol(d_ode_pars)/nrow(d_ode_pars)
  sim_ode_trajectory(d_ode_pars, d_intercept, X_M[birth_grid_idx, ], d_grid, q = q)
}

#' Exact-index lookup of frames on a grid
#'
#' @param frames the times to look up.
#' @param grid the grid.
#' @param arg_name the argument name used in the error message.
#' @param tol the tolerance.
#' @returns integer indices into `grid`; an error names the first frame that is
#' not on the grid.
#' @noRd
.frame_index <- function(frames,
                         grid,
                         arg_name = "obs_time",
                         tol = 1e-10){
  idx <- vapply(frames, function(tt){
    hits <- which(abs(grid - tt) < tol)
    if(length(hits) == 0) NA_integer_ else hits[1]
  }, integer(1))
  if(anyNA(idx)){
    bad <- frames[which(is.na(idx))[1]]
    stop("`", arg_name, "` contains the time ", format(bad, digits = 15),
         ", which is not on the fine grid (tolerance ", tol, "); frames are read off ",
         "the grid by exact index")
  }
  idx
}

#' The default fine grid derived from the frames
#'
#' @param obs_time the mother's frames.
#' @param steps_per_frame grid steps per frame interval.
#' @returns `seq(obs_time[1], obs_time[n], length.out = (n - 1) * steps_per_frame + 1)`.
#' @noRd
.derive_grid <- function(obs_time,
                         steps_per_frame){
  n <- length(obs_time)
  stopifnot(n >= 2, steps_per_frame == round(steps_per_frame))
  seq(obs_time[1], obs_time[n], length.out = (n - 1)*steps_per_frame + 1)
}

#' Normalize the `sigma` argument to a list of two length-p vectors
#' @noRd
.sigma_list <- function(sigma,
                        p){
  if(is.list(sigma)){
    stopifnot(all(c("M", "D") %in% names(sigma)))
    out <- list(M = as.numeric(sigma$M), D = as.numeric(sigma$D))
  } else {
    sigma <- as.numeric(sigma)
    if(length(sigma) == 1) sigma <- rep(sigma, p)
    out <- list(M = sigma, D = sigma)
  }
  stopifnot(length(out$M) == p, length(out$D) == p, all(out$M >= 0), all(out$D >= 0))
  out
}

#' Draw an n x p error matrix with per-column sd, iid or AR(1)
#'
#' One `rnorm()` call fills the matrix column-major; the AR(1) filter is
#' applied per column with `e_1 = z_1` and `e_t = rho e_{t-1} + sqrt(1 - rho^2)
#' z_t`, so the stationary variance is `sigma^2` from the first frame.
#' @noRd
.draw_errors <- function(n,
                         p,
                         sigma_vec,
                         noise,
                         rho){
  z_mat <- matrix(stats::rnorm(n*p), nrow = n, ncol = p)
  if(noise == "ar1" && n > 1){
    z_mat <- apply(z_mat, 2, function(z){
      c(stats::filter(c(z[1], sqrt(1 - rho^2)*z[-1]), filter = rho, method = "recursive"))
    })
    z_mat <- matrix(z_mat, nrow = n, ncol = p)
  }
  sweep(z_mat, 2, sigma_vec, "*")
}

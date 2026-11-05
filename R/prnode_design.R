#' Integrated polynomial basis of one cell's trajectory at its frames
#'
#' Evaluates `psi(x) = (x, ..., x^q)` on the cell's fine grid for every
#' protein, integrates each column cumulatively from the cell's first grid
#' point by the trapezoid rule with the grid's step, and reads the frames off
#' the grid by exact index. The first row is
#' exactly zero, so `Psi(t_1) = 0`; a frame more than `1e-10` from every grid
#' point is an error naming that frame. The grid must be
#' uniform and must start at the cell's first frame.
#'
#' @param x_mat `length(times_e) x p` trajectory matrix on the fine grid
#' (smoothed or true), columns are proteins.
#' @param times_e the cell's fine grid, increasing and uniform.
#' @param obs_time the cell's frames, each a grid point.
#' @param q degree of the polynomial basis.
#' @returns an `n x pq` matrix, `n = length(obs_time)`, columns regulator-major:
#' column `(k - 1) q + m` is the integral of `x_k^m`.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' grid <- ds$raw$times_e
#' X_smooth <- smooth_md_trajectory(ds$raw$m_obs_time, ds$raw$Y_M, grid)
#' Psi <- make_integrated_basis(X_smooth, grid, ds$raw$m_obs_time, q = 1)
#' dim(Psi)
#' Psi[1, ]            # zero at the first frame
#' round(Psi[2:4, ], 3)
#' @export
make_integrated_basis <- function(x_mat,
                                  times_e,
                                  obs_time,
                                  q = 1){
  x_mat <- as.matrix(x_mat)
  stopifnot(nrow(x_mat) == length(times_e), length(times_e) >= 1, q >= 1, q == round(q),
            !is.unsorted(times_e, strictly = TRUE))
  p <- ncol(x_mat)
  n_grid <- length(times_e)
  if(n_grid >= 2){
    step_vec <- diff(times_e)
    if(max(abs(step_vec - step_vec[1])) > 1e-8*abs(step_vec[1])){
      stop("`times_e` must be uniformly spaced; the integration step is read from its first two entries")
    }
    delta <- step_vec[1]
  } else {
    delta <- NA_real_
  }
  obs_idx <- .frame_index(obs_time, times_e, "obs_time")

  Psi_mat <- matrix(0, nrow = length(obs_time), ncol = p*q)
  for(k in seq_len(p)){
    for(m in seq_len(q)){
      psi_vec <- x_mat[, k]^m
      # cumulative trapezoid: I_i = Delta (sum_{l <= i} psi_l - (psi_1 + psi_i) / 2)
      integral_vec <- if(n_grid >= 2) delta*(cumsum(psi_vec) - (psi_vec[1] + psi_vec)/2) else 0
      integral_vec[1] <- 0
      Psi_mat[, (k - 1)*q + m] <- integral_vec[obs_idx]
    }
  }
  Psi_mat
}

#' Profile out one cell's intercept and baseline slope
#'
#' Forms the nuisance columns `Z = [1, t - t_1]`, their QR factors, and the
#' residuals of the response and of the integrated basis on them. The profiled quantities are what the
#' penalized problem sees; the unprofiled `y` and `Psi` are kept so that the
#' nuisance coefficients `gamma = qr.coef(qr_Z, y - Psi theta)` and the fitted
#' trajectory can be recovered later. Requires at least 3 frames at two
#' distinct times.
#'
#' @param y_mat `n x p` response matrix (or a vector for one protein).
#' @param Psi_mat `n x pq` integrated basis from `make_integrated_basis()`.
#' @param obs_time the cell's `n` frames.
#' @returns a list with `n`, `Psi` (the input basis), `Psi_tilde` (its
#' residual on `Z`), `qr_Z`, `y` (the input response as a matrix), `y_tilde`
#' (its residual on `Z`) and `Z` (`n x 2`).
#' @export
profile_cell <- function(y_mat,
                         Psi_mat,
                         obs_time){
  y_mat <- as.matrix(y_mat)
  Psi_mat <- as.matrix(Psi_mat)
  n <- length(obs_time)
  if(n < 3){
    stop("`obs_time` has ", n, " frames; profiling the intercept and slope needs at least 3")
  }
  if(nrow(y_mat) != n || nrow(Psi_mat) != n){
    stop("`y_mat` (", nrow(y_mat), " rows) and `Psi_mat` (", nrow(Psi_mat),
         " rows) must have one row per frame of `obs_time` (", n, ")")
  }
  Z <- cbind(1, obs_time - obs_time[1])
  qr_Z <- qr(Z)
  if(qr_Z$rank < 2){
    stop("the nuisance columns [1, t - t_1] are not of full rank: at least two distinct frame times are needed")
  }
  list(n = n,
       Psi = Psi_mat,
       Psi_tilde = qr.resid(qr_Z, Psi_mat),
       qr_Z = qr_Z,
       y = y_mat,
       y_tilde = qr.resid(qr_Z, y_mat),
       Z = Z)
}

#' Build the two-block PRNODE design for one mother-daughter pair
#'
#' The design of the PRNODE estimator. Each cell is smoothed onto its fine
#' grid with `smooth_md_trajectory()` (the daughter on
#' `make_daughter_grid()`), unless `x_list` supplies the trajectories (the
#' regression level: the true trajectories, no smoothing); each cell's
#' integrated basis is built and its intercept and slope profiled out; the two
#' cells are stacked, each scaled by the square root of its loss weight.
#'
#' @details
#' Response scaling and loss weights, in the order they are applied:
#' 1. With `bool_scale_response = TRUE`, protein `j` is divided by its pooled
#'    standard deviation `s_j = sd(c(Y_M[, j], Y_D[, j]))` everywhere it
#'    appears: as the response of its own equation and as a regulator of every
#'    equation (the smoothed or supplied trajectories are divided by the same
#'    `s_j` before the basis is built). The coefficient of regulator `k`'s
#'    `m`-th power in equation `j` is therefore fitted in units of
#'    `s_j / s_k^m`; `fit_prnode()` undoes this. `response_scale` records `s`.
#' 2. With `omega = NULL`, the noise standard deviation `sigma_cj` of each cell
#'    and protein is estimated on the **unscaled** series by
#'    `estimate_noise_sd(method = noise_method)` (at least 4 frames per cell)
#'    and divided by `s_j`.
#' 3. Each cell's proteins are pooled, `omega_c = 1 / mean_j (sigma_cj /
#'    s_j)^2`, and the two cells then share one weight, the frame-weighted
#'    pool `w = n / (n_M / omega_M + n_D / omega_D)` with `n = n_M + n_D`
#'    (so `omega_M = omega_D = w`). A user `omega = c(M, D)` replaces steps 2
#'    and 3 and may weight the cells differently.
#'
#' The stacked design is `x = rbind(sqrt(omega_M) [Psi~_M, Psi~_M, 0],
#' sqrt(omega_D) [Psi~_D, 0, Psi~_D])` with three column blocks of `pq`
#' (scaffold, mother deviation, daughter deviation), so the scaffold block is
#' exactly the sum of the two deviation blocks and there are no unpenalized
#' columns; `y = rbind(sqrt(omega_M) y~_M, sqrt(omega_D) y~_D)`.
#'
#' Frames: both cells must be imaged on one lattice. The daughter's frames
#' must be the mother's frames from the daughter's birth to the mother's last
#' frame, and every frame must lie on the fine grid to `1e-10`. With
#' `times_e = NULL` the grid is derived from the mother's frames,
#' `seq(t_1, t_n, length.out = (n - 1) steps_per_frame + 1)`, so equally
#' spaced frames (the simulations' `seq(0, 1, length.out = 48)`, or real
#' frames every 12 minutes) lie on the grid to within `1e-10` (not bit for
#' bit) and unequally spaced frames are refused naming the first offending
#' frame.
#'
#' @param m_obs_time,Y_M the mother's frames and `n_M x p` observations.
#' @param d_obs_time,Y_D the daughter's frames and `n_D x p` observations;
#' `d_obs_time[1]` is her birth time.
#' @param bool_scale_response divide each protein by its pooled sd (step 1).
#' @param noise_method the `method` of `estimate_noise_sd()` for step 2.
#' @param omega `NULL` (steps 2-3) or a length-2 vector `c(M, D)` of loss
#' weights.
#' @param q degree of the polynomial basis.
#' @param spar smoothing parameter of `smooth_md_trajectory()`.
#' @param steps_per_frame grid steps per frame interval when `times_e` is
#' `NULL`.
#' @param times_e the mother's fine grid, uniform, starting at her first
#' frame; `NULL` derives it from the frames.
#' @param x_list `NULL`, or `list(M = , D = )` of trajectories on the mother's
#' grid and on `make_daughter_grid(d_obs_time[1], times_e)` (the regression
#' level), in the unscaled units of `Y_M` and `Y_D`.
#' @returns a list (the `public` element of a dataset) with `D_groups`
#' (`2p + 1:p`), `d_obs_time`, `index` (`rep(1:(3p), each = q)`), `M_groups`
#' (`p + 1:p`), `m_obs_time`, `mean_state` (list `M`, `D`: the column means
#' of `Y_M` and of `Y_D` in their unscaled units, the states at which
#' `inheritance_flow()` linearizes each cell's fitted field when `q > 1`),
#' `n_D`, `n_M`, `omega` (named `M`, `D`), `p`,
#' `profile` (list `M`, `D` of `profile_cell()` outputs on the scaled series,
#' each with the added `regulator_scale`, the per-regulator `s_k`),
#' `protein_names`, `q`, `response_scale` (length `p`), `spar` (the smoothing level; `NA` when `x_list` supplies the trajectories), `times_e` (the grid
#' used), `x` (`(n_M + n_D) x 3pq`) and `y` (`(n_M + n_D) x p`).
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' # from the observed matrices alone, as for a real mother-daughter pair
#' design <- make_prnode_design(m_obs_time = ds$raw$m_obs_time,
#'                              Y_M = ds$raw$Y_M,
#'                              d_obs_time = ds$raw$d_obs_time,
#'                              Y_D = ds$raw$Y_D, steps_per_frame = 100)
#' dim(design$x)
#' dim(design$y)
#' design$index
#' design$M_groups
#' design$D_groups
#' design$omega
#' all.equal(design$x, ds$public$x)   # the dataset's public design
#' @export
make_prnode_design <- function(m_obs_time,
                               Y_M,
                               d_obs_time,
                               Y_D,
                               bool_scale_response = TRUE,
                               noise_method = c("diff2", "spline"),
                               omega = NULL,
                               q = 1,
                               spar = 0.5,
                               steps_per_frame = 2000,
                               times_e = NULL,
                               x_list = NULL){
  noise_method <- match.arg(noise_method)
  Y_M <- as.matrix(Y_M)
  Y_D <- as.matrix(Y_D)
  stopifnot(is.numeric(Y_M), is.numeric(Y_D), length(q) == 1, q >= 1, q == round(q),
            length(bool_scale_response) == 1, is.logical(bool_scale_response))
  if(ncol(Y_M) != ncol(Y_D)){
    stop("`ncol(Y_M)` = ", ncol(Y_M), " and `ncol(Y_D)` = ", ncol(Y_D),
         " differ; both cells must carry the same proteins")
  }
  if(nrow(Y_M) != length(m_obs_time)){
    stop("`Y_M` has ", nrow(Y_M), " rows but `m_obs_time` has ", length(m_obs_time), " frames")
  }
  if(nrow(Y_D) != length(d_obs_time)){
    stop("`Y_D` has ", nrow(Y_D), " rows but `d_obs_time` has ", length(d_obs_time), " frames")
  }
  if(is.unsorted(m_obs_time, strictly = TRUE) || is.unsorted(d_obs_time, strictly = TRUE)){
    stop("`m_obs_time` and `d_obs_time` must be strictly increasing")
  }
  p <- ncol(Y_M)
  n_M <- nrow(Y_M)
  n_D <- nrow(Y_D)
  for(cell in c("mother", "daughter")){
    n_cell <- if(cell == "mother") n_M else n_D
    if(n_cell < 3){
      stop("the ", cell, " has ", n_cell, " frames; the design needs at least 3 per cell")
    }
  }
  if(!is.null(omega)){
    if(length(omega) != 2 || !is.numeric(omega) || any(omega <= 0)){
      stop("`omega` must be a length-2 positive numeric vector c(M, D), not of length ", length(omega))
    }
  }
  if(is.null(times_e)){
    times_e <- .derive_grid(m_obs_time, steps_per_frame)
  }
  stopifnot(length(times_e) >= 2, !is.unsorted(times_e, strictly = TRUE))
  # every frame on the grid (named by the first offending frame), before the
  # lattice check, so that unequally spaced mother frames are reported as such
  .frame_index(m_obs_time, times_e, "m_obs_time")
  d_times_e <- make_daughter_grid(d_obs_time[1], times_e)
  .frame_index(d_obs_time, d_times_e, "d_obs_time")
  # one lattice: the daughter's frames are the mother's from her birth on
  m_suffix <- m_obs_time[m_obs_time >= d_obs_time[1] - 1e-10]
  bool_lattice <- length(m_suffix) == n_D && all(abs(m_suffix - d_obs_time) < 1e-10)
  if(!bool_lattice){
    stop("`d_obs_time` must be the mother's frames from the daughter's birth to the mother's ",
         "last frame (`m_obs_time[m_obs_time >= d_obs_time[1]]`); a daughter whose frames ",
         "precede the mother's first frame, stop before her last frame, or fall between her ",
         "frames cannot share the mother's lattice")
  }
  protein_names <- colnames(Y_M)
  if(is.null(protein_names)) protein_names <- paste0("protein_", seq_len(p))

  # trajectories on the grids: supplied (regression level) or smoothed
  if(is.null(x_list)){
    if(n_M < 4 || n_D < 4){
      stop("smoothing needs at least 4 frames per cell; supply `x_list` for shorter cells")
    }
    X_M <- smooth_md_trajectory(m_obs_time, Y_M, times_e, spar = spar)
    X_D <- smooth_md_trajectory(d_obs_time, Y_D, d_times_e, spar = spar)
  } else {
    if(!is.list(x_list) || !all(c("M", "D") %in% names(x_list))){
      stop("`x_list` must be a list with elements `M` and `D`")
    }
    X_M <- as.matrix(x_list$M)
    X_D <- as.matrix(x_list$D)
    if(nrow(X_M) != length(times_e) || ncol(X_M) != p){
      stop("`x_list$M` must be a ", length(times_e), " x ", p, " matrix on the mother's grid, not ",
           nrow(X_M), " x ", ncol(X_M))
    }
    if(nrow(X_D) != length(d_times_e) || ncol(X_D) != p){
      stop("`x_list$D` must be a ", length(d_times_e), " x ", p, " matrix on the daughter's grid ",
           "(the mother's grid from her birth on), not ", nrow(X_D), " x ", ncol(X_D))
    }
  }

  # step 1: per-protein response scaling, applied to the series as response
  # and as regulator
  response_scale <- if(bool_scale_response){
    vapply(seq_len(p), function(j) stats::sd(c(Y_M[, j], Y_D[, j])), numeric(1))
  } else {
    rep(1, p)
  }
  if(any(!is.finite(response_scale)) || any(response_scale <= 0)){
    stop("every protein needs a positive pooled standard deviation for `bool_scale_response = TRUE`")
  }
  Y_M_scaled <- sweep(Y_M, 2, response_scale, "/")
  Y_D_scaled <- sweep(Y_D, 2, response_scale, "/")
  X_M_scaled <- sweep(X_M, 2, response_scale, "/")
  X_D_scaled <- sweep(X_D, 2, response_scale, "/")

  # steps 2-3: noise sds on the unscaled series, divided by s_j, pooled per cell
  if(is.null(omega)){
    for(cell in c("mother", "daughter")){
      n_cell <- if(cell == "mother") n_M else n_D
      if(n_cell < 4){
        stop("the ", cell, " has ", n_cell, " frames; estimating the noise sd for the loss weights ",
             "needs at least 4 frames per cell. Pass `omega = c(M, D)` to fit a shorter cell")
      }
    }
    sigma_M <- estimate_noise_sd(m_obs_time, Y_M, method = noise_method, spar = spar)
    sigma_D <- estimate_noise_sd(d_obs_time, Y_D, method = noise_method, spar = spar)
    omega_cell <- c(M = 1/mean((sigma_M/response_scale)^2), D = 1/mean((sigma_D/response_scale)^2))
    # one weight for both cells: the frame-weighted pool of the per-cell
    # estimates, 1/w = (n_M/omega_M + n_D/omega_D)/n (per-cell weights over-trust the shorter cell's noise estimate)
    w <- (n_M + n_D)/(n_M/omega_cell[["M"]] + n_D/omega_cell[["D"]])
    omega <- c(M = w, D = w)
  } else {
    omega <- c(M = unname(omega[1]), D = unname(omega[2]))
  }

  Psi_M <- make_integrated_basis(X_M_scaled, times_e, m_obs_time, q = q)
  Psi_D <- make_integrated_basis(X_D_scaled, d_times_e, d_obs_time, q = q)
  prof_M <- profile_cell(Y_M_scaled, Psi_M, m_obs_time)
  prof_D <- profile_cell(Y_D_scaled, Psi_D, d_obs_time)
  # the regulator scales travel with the profile, so a design whose `y` and
  # `response_scale` are later subset to some proteins still unscales every
  # regulator column (test plan F2)
  prof_M$regulator_scale <- response_scale
  prof_D$regulator_scale <- response_scale

  pq <- p*q
  xm <- sqrt(omega[["M"]])*prof_M$Psi_tilde
  xd <- sqrt(omega[["D"]])*prof_D$Psi_tilde
  zero_M <- matrix(0, n_M, pq)
  zero_D <- matrix(0, n_D, pq)
  x <- rbind(cbind(xm, xm, zero_M),
             cbind(xd, zero_D, xd))
  y <- rbind(sqrt(omega[["M"]])*prof_M$y_tilde,
             sqrt(omega[["D"]])*prof_D$y_tilde)
  dimnames(x) <- NULL
  dimnames(y) <- NULL
  stopifnot(nrow(x) == n_M + n_D, ncol(x) == 3*pq, nrow(y) == n_M + n_D, ncol(y) == p)

  list(D_groups = 2*p + seq_len(p),
       d_obs_time = d_obs_time,
       index = rep(seq_len(3*p), each = q),
       M_groups = p + seq_len(p),
       m_obs_time = m_obs_time,
       mean_state = list(M = colMeans(Y_M), D = colMeans(Y_D)),
       n_D = n_D,
       n_M = n_M,
       omega = omega,
       p = p,
       profile = list(M = prof_M, D = prof_D),
       protein_names = protein_names,
       q = q,
       response_scale = response_scale,
       spar = if(is.null(x_list)) spar else NA_real_,
       times_e = times_e,
       x = x,
       y = y)
}

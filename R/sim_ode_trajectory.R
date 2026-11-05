#' Solve one cell's ODE on a time grid
#'
#' The noiseless trajectory of `dx_j/dt = intercept_j + sum_k psi(x_k)' theta_jk`
#' with `psi(x) = (x, ..., x^q)`, started from `init_vals` at `times_e[1]`.
#' For `q = 1` the system is
#' linear and is solved exactly by the matrix exponential of the augmented
#' matrix `G = [A, b; 0, 0]`: on a uniform grid `E = expm(Delta G)` is formed
#' once and its powers are applied by repeated squaring, so the round-off
#' grows with `log2(length(times_e))` rather than with the number of steps;
#' on a non-uniform grid one exponential is formed per step. For `q > 1`, or
#' when `method = "deSolve"`, the system is integrated by
#' `deSolve::ode(method = "lsoda", rtol = 1e-10, atol = 1e-12)`; a solver
#' warning or a non-finite state is an error naming `deSolve`.
#'
#' There is no noise and no observation step here: `sim_md_trajectory()` reads
#' the frames off the grid by exact index.
#'
#' @param ode_pars `p x pq` coefficient matrix, row = equation, column
#' `(k - 1) q + m` = regulator `k`'s `m`-th power (the layout of
#' `sim_md_network()`).
#' @param intercept length-`p` vector of intercepts.
#' @param init_vals length-`p` state at `times_e[1]`.
#' @param times_e increasing time grid; a one-point grid returns `init_vals`.
#' @param q degree of the polynomial basis; must equal `ncol(ode_pars) / p`.
#' @param method `"expm"` (exact, `q = 1` only) or `"deSolve"`.
#' @param ... further arguments to `deSolve::ode()`.
#' @returns a base `length(times_e) x p` matrix with `X[1, ] == init_vals`.
#' @examples
#' # a two-protein oscillator, x1' = -x2 and x2' = x1: row j of `ode_pars`
#' # holds the regulators of protein j, so one turn takes 2 pi
#' A <- matrix(c(0, 1, -1, 0), 2, 2)
#' times <- seq(0, 2*pi, length.out = 101)
#' X <- sim_ode_trajectory(A, intercept = c(0, 0), init_vals = c(1, 0),
#'                         times_e = times)
#' round(X[c(1, 26, 51, 101), ], 6)   # start, a quarter, a half and a full turn
#' if(requireNamespace("deSolve", quietly = TRUE)){
#'   X_ode <- sim_ode_trajectory(A, c(0, 0), c(1, 0), times, method = "deSolve")
#'   max(abs(X_ode - X))
#' }
#' @export
sim_ode_trajectory <- function(ode_pars,
                               intercept,
                               init_vals,
                               times_e,
                               q = 1,
                               method = c("expm", "deSolve"),
                               ...){
  method <- match.arg(method)
  ode_pars <- as.matrix(ode_pars)
  p <- nrow(ode_pars)
  intercept <- as.numeric(intercept)
  if(length(init_vals) != p){
    stop("`init_vals` has length ", length(init_vals), " but `ode_pars` has ", p, " rows")
  }
  if(ncol(ode_pars) != p*q){
    stop("`ode_pars` has ", ncol(ode_pars), " columns; p * q = ", p*q, " are needed")
  }
  if(length(intercept) != p){
    stop("`intercept` has length ", length(intercept), " but `ode_pars` has ", p, " rows")
  }
  if(length(times_e) < 1 || is.unsorted(times_e, strictly = TRUE)){
    stop("`times_e` must be a strictly increasing time grid")
  }
  if(q > 1 && method == "expm"){
    method <- "deSolve"
  }
  init_vals <- as.numeric(init_vals)
  n_grid <- length(times_e)
  if(n_grid == 1) return(matrix(init_vals, nrow = 1, ncol = p))

  if(method == "expm"){
    X <- .solve_linear_ode(ode_pars, intercept, init_vals, times_e)
  } else {
    X <- .solve_ode_desolve(ode_pars, intercept, init_vals, times_e, q, ...)
  }
  X[1, ] <- init_vals
  dimnames(X) <- NULL
  stopifnot(nrow(X) == n_grid, ncol(X) == p)
  X
}

#' Exact solution of a linear ODE by the matrix exponential
#'
#' Uses `(x(t), 1) = expm(t G) (x(0), 1)` with `G = [A, b; 0, 0]`. On a
#' uniform grid the one-step propagator is squared repeatedly and applied to
#' blocks of already-solved points (doubling), so `log2(n)` matrix products
#' replace `n` sequential steps.
#'
#' @param A `p x p` coefficient matrix.
#' @param b length-`p` intercept.
#' @param x0 initial state.
#' @param times_e the grid.
#' @returns `length(times_e) x p` matrix.
#' @noRd
.solve_linear_ode <- function(A,
                              b,
                              x0,
                              times_e){
  p <- length(x0)
  n_grid <- length(times_e)
  G <- rbind(cbind(A, b), 0)
  step_vec <- diff(times_e)
  # seq() grids carry step jitter of about 1e-16 absolute, i.e. 1e-11 relative
  # at a step of 1e-5, so the uniformity test is loose enough to take the fast
  # doubling route on them (one expm instead of one per step)
  bool_uniform <- max(abs(step_vec - step_vec[1])) <= 1e-8*max(abs(step_vec))
  # rows are grid points, columns the augmented state (x, 1)
  Z <- matrix(NA_real_, nrow = n_grid, ncol = p + 1)
  Z[1, ] <- c(x0, 1)
  if(bool_uniform){
    E <- as.matrix(Matrix::expm(step_vec[1]*G))
    num_done <- 1
    E_power <- E
    while(num_done < n_grid){
      num_new <- min(num_done, n_grid - num_done)
      # points num_done + 1 ... num_done + num_new are E^num_done times
      # points 1 ... num_new
      Z[num_done + seq_len(num_new), ] <- Z[seq_len(num_new), , drop = FALSE] %*% t(E_power)
      num_done <- num_done + num_new
      if(num_done < n_grid) E_power <- E_power %*% E_power
    }
  } else {
    for(i in seq_len(n_grid - 1)){
      E_i <- as.matrix(Matrix::expm(step_vec[i]*G))
      Z[i + 1, ] <- c(E_i %*% Z[i, ])
    }
  }
  Z[, seq_len(p), drop = FALSE]
}

#' Numerical solution by deSolve
#'
#' @param ode_pars,intercept,init_vals,times_e,q as in `sim_ode_trajectory()`.
#' @param ... passed to `deSolve::ode()`.
#' @returns `length(times_e) x p` matrix.
#' @noRd
.solve_ode_desolve <- function(ode_pars,
                               intercept,
                               init_vals,
                               times_e,
                               q,
                               ...){
  if(!requireNamespace("deSolve", quietly = TRUE)){
    stop("`deSolve` is needed for q > 1 or method = \"deSolve\"; please install it")
  }
  p <- length(init_vals)
  rhs <- function(t, x, parms){
    # psi in regulator-major order: regulator k's powers 1..q are adjacent
    psi_vec <- as.vector(t(outer(x, seq_len(q), "^")))
    list(intercept + c(ode_pars %*% psi_vec))
  }
  out <- withCallingHandlers(
    tryCatch(deSolve::ode(y = init_vals, times = times_e, func = rhs, parms = NULL,
                          method = "lsoda", rtol = 1e-10, atol = 1e-12, ...),
             error = function(e) stop("deSolve failed: ", conditionMessage(e))),
    warning = function(w) stop("deSolve reported a problem (the system may be unstable or ",
                               "non-finite): ", conditionMessage(w)))
  X <- unname(as.matrix(out)[, 1 + seq_len(p), drop = FALSE])
  if(nrow(X) != length(times_e) || !all(is.finite(X))){
    stop("deSolve returned a non-finite or truncated trajectory")
  }
  X
}

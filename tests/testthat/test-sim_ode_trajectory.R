# Plan section 2: the exact ODE solver sim_ode_trajectory() (R/sim_ode_trajectory.R, new).

# the triplet system of the retired generator (sim_md_triplets(), deleted
# 2026-10-22): A antisymmetric, intercepts
# (-2A, 2A, 0); its closed form is a rotation at num_cycles cycles per unit
# time about the conserved axis (1, 1, 1)
.triplet_system <- function(num_cycles){
  A0 <- num_cycles/sqrt(3)*2*pi
  A <- A0*matrix(c(0, 1, -1,
                   -1, 0, 1,
                   1, -1, 0), 3, 3, byrow = TRUE)
  list(A = A, b = A0*c(-2, 2, 0), omega = A0*sqrt(3))
}

.triplet_closed_form <- function(sys, x0, times){
  # x* solves A x* = -b within the plane orthogonal to the kernel (1, 1, 1):
  # derived, not typed, after the audit of 2026-10-21 found the constant
  # (0, 0, 2) off by a sign (the solution is (0, 0, -2)); d/dt (x - x*) =
  # A (x - x*); split x0 - x* into the kernel direction and the rotation plane
  x_star <- -solve(sys$A + tcrossprod(rep(1, 3))/3, sys$b)
  dev <- x0 - x_star
  kernel <- rep(mean(dev), 3)
  u0 <- dev - kernel
  Au0 <- c(sys$A %*% u0)
  t(sapply(times, function(tt) x_star + kernel + cos(sys$omega*tt)*u0 + sin(sys$omega*tt)*Au0/sys$omega))
}

## S1: two routes agree. The expm route and the deSolve route agree on a
## stable random system and on the triplet system, and both agree with the
## triplet's closed form. The closed form is the oracle that neither solver
## shares.
test_that("S1: expm and deSolve agree to 1e-8 relative, and both match the triplet closed form", {
  .skip_unless_built("sim_ode_trajectory")
  .skip_unless_built("sim_md_network")
  skip_if_not_installed("deSolve")
  times <- seq(0, 1, by = 1e-4)
  truth <- sim_md_network(6, 1, 6, 2, 2, seed_number = 21)
  set.seed(21)
  x0 <- stats::rnorm(6)
  X_expm <- sim_ode_trajectory(truth$m_ode_pars, truth$m_intercept, x0, times, q = 1, method = "expm")
  X_ode <- sim_ode_trajectory(truth$m_ode_pars, truth$m_intercept, x0, times, q = 1, method = "deSolve")
  expect_true(is.matrix(X_expm) && !inherits(X_expm, "Matrix"))
  expect_equal(dim(X_expm), c(length(times), 6))
  expect_equal(X_expm[1, ], x0)
  expect_equal(X_ode, X_expm, tolerance = 1e-8)

  for(num_cycles in c(3, 6)){
    sys <- .triplet_system(num_cycles)
    x0 <- c(-1, 0.5, 1.5)
    X_cf <- .triplet_closed_form(sys, x0, times)
    X_expm <- sim_ode_trajectory(sys$A, sys$b, x0, times, q = 1, method = "expm")
    X_ode <- sim_ode_trajectory(sys$A, sys$b, x0, times, q = 1, method = "deSolve")
    label <- as.character(num_cycles)
    expect_equal(X_expm, X_cf, tolerance = 1e-8, info = label)
    expect_equal(X_ode, X_cf, tolerance = 1e-8, info = label)
    # periodicity at num_cycles cycles per unit time, read at t = 1 (a whole
    # number of cycles; 1/num_cycles is not on the 1e-4 grid) and conservation
    # of x1 + x2 + x3
    expect_equal(X_expm[length(times), ], X_expm[1, ], tolerance = 1e-8, info = label)
    expect_equal(rowSums(X_expm), rep(sum(x0), nrow(X_expm)), tolerance = 1e-10, info = label)
  }
})

## S2: closed forms a solver cannot get wrong by accident.
test_that("S2: a diagonal A gives exponentials to 1e-12 and A = 0 gives a straight line", {
  .skip_unless_built("sim_ode_trajectory")
  times <- seq(0, 1, by = 1e-4)
  a_vec <- c(-0.5, 0.3, -2)
  x0 <- c(1, 2, 3)
  X <- sim_ode_trajectory(diag(a_vec), rep(0, 3), x0, times, q = 1)
  expect_equal(X, sapply(1:3, function(k) x0[k]*exp(a_vec[k]*times)), tolerance = 1e-12)
  b_vec <- c(1, -2, 0.5)
  X0 <- sim_ode_trajectory(matrix(0, 3, 3), b_vec, x0, times, q = 1)
  expect_equal(X0, sapply(1:3, function(k) x0[k] + b_vec[k]*times), tolerance = 1e-10)
  # a one-point grid returns the initial values
  expect_equal(sim_ode_trajectory(diag(a_vec), rep(0, 3), x0, 0.3, q = 1), matrix(x0, 1, 3), ignore_attr = TRUE)
})

## S3: a non-uniform grid (the daughter's shorter first step when a birth time
## is off the grid) gives the same values at common points as the uniform grid.
test_that("S3: a non-uniform grid agrees with the uniform grid at the common times to 1e-10", {
  .skip_unless_built("sim_ode_trajectory")
  .skip_unless_built("sim_md_network")
  truth <- sim_md_network(4, 1, 4, 1, 1, seed_number = 31)
  x0 <- c(0.5, -0.5, 1, 0)
  uniform <- seq(0, 1, by = 1e-4)
  nonuni <- c(0, 0.3, 0.3 + 1e-4, seq(0.3 + 2e-4, 1, by = 1e-4))
  X_u <- sim_ode_trajectory(truth$m_ode_pars, truth$m_intercept, x0, uniform, q = 1)
  X_n <- sim_ode_trajectory(truth$m_ode_pars, truth$m_intercept, x0, nonuni, q = 1)
  common_u <- which(uniform %in% nonuni | sapply(uniform, function(tt) any(abs(nonuni - tt) < 1e-12)))
  common_n <- sapply(uniform[common_u], function(tt) which(abs(nonuni - tt) < 1e-12))
  expect_gt(length(common_u), 100)
  expect_equal(X_n[common_n, ], X_u[common_u, ], tolerance = 1e-10)
  # the exact matrix exponential at the exact times, as a second oracle
  G <- rbind(cbind(truth$m_ode_pars, truth$m_intercept), 0)
  exact <- t(sapply(nonuni, function(tt) as.matrix(Matrix::expm(tt*G))[1:4, ] %*% c(x0, 1)))
  expect_equal(X_n, exact, tolerance = 1e-10, ignore_attr = TRUE)
})

## S4: guards, each by message; the q > 1 system x' = x^2 from x(0) = 10 blows
## up at t = 0.1, which deSolve must report as an error, not a warning
## (pitfalls row 46).
test_that("S4: unsorted times, wrong lengths and a blow-up each error by message", {
  .skip_unless_built("sim_ode_trajectory")
  skip_if_not_installed("deSolve")
  A <- diag(c(-1, -1))
  expect_error(sim_ode_trajectory(A, c(0, 0), c(1, 1), c(0, 0.5, 0.2), q = 1), "times_e")
  expect_error(sim_ode_trajectory(A, c(0, 0), c(1, 1, 1), seq(0, 1, by = 0.1), q = 1), "init_vals")
  expect_error(sim_ode_trajectory(matrix(0, 2, 3), c(0, 0), c(1, 1), seq(0, 1, by = 0.1), q = 1), "ode_pars")
  expect_error(sim_ode_trajectory(matrix(c(0, 1), 1, 2), 0, 10, seq(0, 1, by = 1e-3), q = 2, method = "deSolve"),
               "deSolve|finite")
})

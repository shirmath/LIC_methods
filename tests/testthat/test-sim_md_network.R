# Plan section 1: the truth generator sim_md_network() (R/sim_md_network.R, new).
# Every test skips until the function is built.

.edge_types <- c("absent", "shared_equal", "shared_unequal", "mother_only", "daughter_only", "self")

.draw_truth <- function(p, q = 1, num_shared, num_mother_only, num_daughter_only, seed = 1, ...){
  sim_md_network(p = p, q = q, num_shared = num_shared, num_mother_only = num_mother_only,
                 num_daughter_only = num_daughter_only, seed_number = seed, ...)
}

## N1: the counts of each edge type are exactly what was asked for, including
## zero counts of a type (pitfalls row 29: all() of an empty comparison is
## TRUE, so the counts are compared as integers, never as all(...)).
test_that("N1: edge counts, the self-edges, and edge_type agree position by position", {
  .skip_unless_built("sim_md_network")
  grid <- expand.grid(p = c(3, 6), q = c(1, 3), num_shared = c(0, 2), num_mother_only = c(0, 1),
                      num_daughter_only = c(0, 2))
  for(r in seq_len(nrow(grid))){
    g <- grid[r, ]
    label <- paste(names(g), unlist(g), collapse = " ")
    truth <- .draw_truth(g$p, g$q, g$num_shared, g$num_mother_only, g$num_daughter_only, seed = r)
    m_net <- .as_adj(truth$m_network)
    d_net <- .as_adj(truth$d_network)
    s_net <- .as_adj(truth$s_network)
    off <- row(m_net) != col(m_net)
    expect_equal(dim(truth$edge_type), c(g$p, g$p), info = label)
    expect_true(all(truth$edge_type %in% .edge_types), info = label)
    expect_true(all(diag(truth$edge_type) == "self"), info = label)
    expect_equal(sum(truth$edge_type[off] != "absent"), g$num_shared + g$num_mother_only + g$num_daughter_only,
                 info = label)
    expect_equal(sum(truth$edge_type == "mother_only"), g$num_mother_only, info = label)
    expect_equal(sum(truth$edge_type == "daughter_only"), g$num_daughter_only, info = label)
    expect_equal(sum(truth$edge_type %in% c("shared_equal", "shared_unequal")), g$num_shared, info = label)
    expect_equal(sum(m_net) - g$p, g$num_shared + g$num_mother_only, info = label)
    expect_equal(sum(d_net) - g$p, g$num_shared + g$num_daughter_only, info = label)
    expect_true(all(diag(m_net)) && all(diag(d_net)), info = label)
    m_from_type <- matrix(truth$edge_type %in% c("shared_equal", "shared_unequal", "mother_only", "self"), g$p, g$p)
    d_from_type <- matrix(truth$edge_type %in% c("shared_equal", "shared_unequal", "daughter_only", "self"), g$p, g$p)
    expect_identical(m_net, m_from_type, info = label)
    expect_identical(d_net, d_from_type, info = label)
    expect_identical(s_net, m_net & d_net, info = label)
    expect_equal(dim(truth$m_ode_pars), c(g$p, g$p*g$q), info = label)
    expect_equal(dim(truth$d_ode_pars), c(g$p, g$p*g$q), info = label)
    expect_equal(c(truth$p, truth$q), c(g$p, g$q), info = label)
  }
})

## N2: the decomposition theta = s + d holds exactly and each type has its
## sign structure; the regulator block of edge (j, k) is columns (k-1)q + 1:q.
test_that("N2: m_ode_pars = s_pars + m_delta, d_ode_pars = s_pars + d_delta, with the type's sign structure", {
  .skip_unless_built("sim_md_network")
  for(q in c(1, 3)){
    truth <- .draw_truth(6, q, num_shared = 8, num_mother_only = 3, num_daughter_only = 3, seed = 7, frac_unequal = 0.5)
    expect_identical(truth$m_ode_pars, truth$s_pars + truth$m_delta, info = q)
    expect_identical(truth$d_ode_pars, truth$s_pars + truth$d_delta, info = q)
    block <- function(mat, j, k) mat[j, (k - 1)*q + seq_len(q)]
    for(j in 1:6) for(k in 1:6){
      type <- truth$edge_type[j, k]
      s_blk <- block(truth$s_pars, j, k)
      m_blk <- block(truth$m_delta, j, k)
      d_blk <- block(truth$d_delta, j, k)
      label <- paste(q, j, k, type)
      if(type == "absent"){
        expect_true(all(c(s_blk, m_blk, d_blk) == 0), info = label)
      } else if(type %in% c("shared_equal", "self")){
        expect_true(all(s_blk != 0) && all(m_blk == 0) && all(d_blk == 0), info = label)
        if(type == "self") expect_true(all(s_blk < 0), info = label)
      } else if(type == "shared_unequal"){
        expect_true(all(s_blk != 0), info = label)
        expect_true(xor(any(m_blk != 0), any(d_blk != 0)), info = label)
        nz <- if(any(m_blk != 0)) m_blk else d_blk
        expect_true(all(sign(nz[nz != 0]) == sign(s_blk[nz != 0])), info = label)
      } else if(type == "mother_only"){
        expect_true(all(s_blk == 0) && all(m_blk != 0) && all(d_blk == 0), info = label)
      } else if(type == "daughter_only"){
        expect_true(all(s_blk == 0) && all(m_blk == 0) && all(d_blk != 0), info = label)
      }
    }
  }
})

## N3 resolved as (b): the diagonal is drawn nonzero in both cells (N4 depends
## on it) and scored on request only; here only its presence is pinned.
test_that("N3(b): the diagonal is negative in both cells and equal across cells", {
  .skip_unless_built("sim_md_network")
  truth <- .draw_truth(6, 1, 6, 2, 2, seed = 3, self_decay = c(0.5, 2))
  expect_true(all(diag(truth$m_ode_pars) < 0))
  expect_true(all(diag(truth$m_ode_pars) >= -2) && all(diag(truth$m_ode_pars) <= -0.5))
  expect_identical(diag(truth$m_ode_pars), diag(truth$d_ode_pars))
})

## N4: no conserved quantity (what Z1 was the failure of). The augmented matrix
## [A, b; 0, 0] has a one-dimensional left null space (only the trivial
## e_{p+1}) iff A is nonsingular, and the profiled design built from the true
## trajectory has full column rank with a well-separated smallest singular
## value, for both cells and every accepted draw.
test_that("N4: every draw has a nonsingular A in both cells and a full-rank profiled design", {
  .skip_unless_built("sim_md_network")
  .skip_unless_built("sim_ode_trajectory")
  .skip_unless_built("make_integrated_basis")
  .skip_unless_built("profile_cell")
  p <- 6
  grid <- .fx_times_coarse
  obs_time <- .fx_obs_time
  birth_time <- obs_time[48 - 12 + 1]
  for(seed in seq_len(20)){
    truth <- .draw_truth(p, 1, 6, 2, 2, seed = seed)
    for(cell in c("M", "D")){
      A <- if(cell == "M") truth$m_ode_pars else truth$d_ode_pars
      b <- if(cell == "M") truth$m_intercept else truth$d_intercept
      G <- rbind(cbind(A, b), 0)
      expect_equal(qr(A)$rank, p, info = paste(seed, cell))
      expect_equal(qr(t(G))$rank, p, info = paste(seed, cell))
      set.seed(seed)
      init_vals <- stats::rnorm(p)
      if(cell == "M"){
        cell_grid <- grid
        cell_obs <- obs_time
      } else {
        cell_grid <- grid[grid >= birth_time - 1e-12]
        cell_obs <- obs_time[obs_time >= birth_time]
      }
      X <- sim_ode_trajectory(A, b, init_vals, cell_grid, q = 1)
      Psi <- make_integrated_basis(X, cell_grid, cell_obs, q = 1)
      obs_idx <- sapply(cell_obs, function(tt) which(abs(cell_grid - tt) < 1e-10))
      prof <- profile_cell(X[obs_idx, , drop = FALSE], Psi, cell_obs)
      sv <- svd(prof$Psi_tilde, nu = 0, nv = 0)$d
      expect_equal(qr(prof$Psi_tilde)$rank, p, info = paste(seed, cell))
      expect_gt(min(sv)/max(sv), 1e-8, label = paste("singular-value ratio", seed, cell))
    }
  }
})

## N5: a design constant, not a property: the share of draws whose A_c is
## stable with the default self_decay. Pinned as a snapshot with the seed so
## that a change to the magnitude ranges is noticed.
test_that("N5: at least 90% of 200 default draws are stable in both cells [snapshot]", {
  .skip_unless_built("sim_md_network")
  stable <- sapply(seq_len(200), function(seed){
    truth <- .draw_truth(6, 1, 6, 2, 2, seed = 1000 + seed)
    all(Re(eigen(truth$m_ode_pars, only.values = TRUE)$values) < 0) &&
      all(Re(eigen(truth$d_ode_pars, only.values = TRUE)$values) < 0)
  })
  rate <- mean(stable)
  expect_gte(rate, 0.9)
  expect_snapshot(cat("stable share over seeds 1001-1200:", rate))
})

## N6: the RNG contract. The exact number of draws is an implementation detail
## (sample(), the unequal factors and the signs all consume the stream), so
## what is pinned is reproducibility: a seed gives identical truths, and with
## seed_number = NULL the caller's stream is consumed deterministically and
## left advanced.
test_that("N6: same seed_number gives identical truths; seed_number = NULL leaves the stream to the caller", {
  .skip_unless_built("sim_md_network")
  t1 <- .draw_truth(6, 1, 6, 2, 2, seed = 11)
  t2 <- .draw_truth(6, 1, 6, 2, 2, seed = 11)
  expect_identical(t1, t2)
  t3 <- .draw_truth(6, 1, 6, 2, 2, seed = 12)
  expect_false(identical(t1$m_ode_pars, t3$m_ode_pars))

  set.seed(5)
  seed_before <- .Random.seed
  u1 <- sim_md_network(6, 1, 6, 2, 2, seed_number = NULL)
  seed_after <- .Random.seed
  next1 <- stats::runif(1)
  set.seed(5)
  u2 <- sim_md_network(6, 1, 6, 2, 2, seed_number = NULL)
  next2 <- stats::runif(1)
  expect_identical(u1, u2)
  expect_identical(next1, next2)
  expect_false(identical(seed_before, seed_after))
})

## N7: argument guards, each by message.
test_that("N7: impossible counts, a bad frac_unequal, non-integer counts and q < 1 error by message", {
  .skip_unless_built("sim_md_network")
  expect_error(sim_md_network(3, 1, num_shared = 5, num_mother_only = 1, num_daughter_only = 1), "6")
  expect_error(sim_md_network(3, 1, 1, 1, 1, frac_unequal = 1.5), "frac_unequal")
  expect_error(sim_md_network(3, 1, 1, 1, 1, frac_unequal = -0.1), "frac_unequal")
  expect_error(sim_md_network(3, 1, 1.5, 1, 1), "num_shared")
  expect_error(sim_md_network(3, 0, 1, 1, 1), "q")
})

# Plan section 10, G1 only: the fork grplasso.default(), now a test helper
# (helper-grplasso_fork.R), reproduces the reference solver. G2-G8
# and test-residualize.R, which tested the fork itself, were removed on
# 2026-10-23. The fork's design and fit helpers are in
# helper-fixtures.R (.fork_design(), .fork_fit_protein()).

## G1 (T13): the fork, called as section 6.3 of the derivation notes prescribes, reproduces
## the reference solver; the only reason the fork stays. Compared on the
## objective, the fitted values and the supports rather than coefficient by
## coefficient (review pass 6): the fork's relative-change stopping rule at
## tol = 1e-12 (defect Z7) leaves it at an objective higher by up to 2e-7
## relative, and in the flat valleys of the q = 3 designs (coefficients of
## order 1e4) that gap moves coefficients by up to 4% while the fitted values
## move by at most 6e-5 of the response norm (an objective gap of 1e-7 is a
## fitted gap of about its square root). The reference's objective is never
## the larger one. The
## objective is written in the design's coordinates, 0.5 ||y - x beta||^2 +
## lambda sum_g sqrt(q) mu_g ||x_g beta_g||, and checked against the
## solver's own `objective` first.
test_that("G1: the fork at tol = 1e-12 reproduces fit_prnode_path() in objective, fitted values and supports", {
  .skip_unless_built("fit_prnode_path")
  cases <- list(reg_q1 = .fixture_dataset("regression", num_d = 12)$public,
                ode_q1 = .fixture_dataset("ode", num_d = 12)$public,
                reg_q3 = .fixture_design_q3(12),
                ode_q3 = .fixture_design_q3(12, level = "ode"))
  mu_list <- list(c(S = 1, M = 1, D = 1), c(S = 1, M = 0.5, D = 0.5), c(S = 1, M = 2, D = 1))
  edge_tol <- 1e-8
  num_masked <- 0
  num_compared <- 0
  for(case in names(cases)) for(mu in mu_list){
    design <- cases[[case]]
    label <- paste(case, paste(mu, collapse = ","))
    p <- design$p
    q <- design$q
    pq <- p*q
    lambda_vec <- .lambda_path(design, mu = mu, num_lambda = 3, num_decades = 1.5)
    ref <- fit_prnode_path(design, lambda = lambda_vec, mu = mu, tol = 1e-12, edge_tol = edge_tol)
    group_cols <- lapply(seq_len(3*p), function(g) (g - 1)*q + seq_len(q))
    group_mu <- rep(c(mu[["S"]], mu[["M"]], mu[["D"]]), each = p)
    objective_of <- function(beta, j, lambda){
      resid <- design$y[, j] - design$x %*% beta
      penalty <- sum(sapply(seq_along(group_cols), function(g){
        cols <- group_cols[[g]]
        sqrt(q)*group_mu[g]*sqrt(sum((design$x[, cols, drop = FALSE] %*% beta[cols])^2))
      }))
      0.5*sum(resid^2) + lambda*penalty
    }
    # the fitted-contribution norms of theta_c, regulator by regulator, on the
    # cell's rows and block, relative to the cell's response norm
    contribution_of <- function(theta, cell, j){
      rows <- if(cell == "M") seq_len(design$n_M) else design$n_M + seq_len(design$n_D)
      block <- (if(cell == "M") pq else 2*pq) + seq_len(pq)
      y_norm <- sqrt(sum(design$y[rows, j]^2))
      sapply(seq_len(p), function(k){
        cols <- block[(k - 1)*q + seq_len(q)]
        sqrt(sum((design$x[rows, cols, drop = FALSE] %*% theta[(k - 1)*q + seq_len(q)])^2))/y_norm
      })
    }
    for(j in seq_len(p)){
      fd <- .fork_design(design)
      lambda_fork <- 2*sqrt(fd$n)*mu[["S"]]*lambda_vec
      fork_fit <- suppressWarnings(
        grplasso.default(x = fd$x, y = design$y[, j], index = fd$index, M_groups = fd$M_groups,
                         D_groups = fd$D_groups, lambda = lambda_fork, model = LinReg(), center = FALSE,
                         standardize = TRUE, residualize = TRUE, mscale = mu[["M"]]/mu[["S"]],
                         dscale = mu[["D"]]/mu[["S"]],
                         control = grpl.control(tol = 1e-12, max.iter = 20000, trace = 0)))
      for(l in seq_along(lambda_vec)){
        info <- paste(label, j, l)
        beta_fork <- c(fork_fit$coefficients[-(1:4), l])
        beta_ref <- .path_beta(ref, j, l)
        obj_ref <- objective_of(beta_ref, j, lambda_vec[l])
        obj_fork <- objective_of(beta_fork, j, lambda_vec[l])
        expect_equal(obj_ref, ref$objective[j, l], tolerance = 1e-8, info = info)
        expect_true(obj_ref <= obj_fork*(1 + 1e-9), info = paste(info, "reference objective larger by",
                                                                 signif((obj_ref - obj_fork)/obj_fork, 3)))
        expect_true(abs(obj_fork - obj_ref) <= 1e-6*obj_ref, info = paste(info, "objective gap",
                                                                       signif((obj_fork - obj_ref)/obj_ref, 3)))
        fitted_gap <- sqrt(sum((design$x %*% (beta_ref - beta_fork))^2))/sqrt(sum(design$y[, j]^2))
        expect_true(fitted_gap <= 1e-4, info = paste(info, "fitted gap", signif(fitted_gap, 3)))
        for(cell in c("M", "D")){
          theta_fork <- beta_fork[seq_len(pq)] + beta_fork[(if(cell == "M") pq else 2*pq) + seq_len(pq)]
          contrib_ref <- contribution_of(.path_mat(ref, paste0("theta_", cell), l)[j, ], cell, j)
          contrib_fork <- contribution_of(theta_fork, cell, j)
          # a regulator within a factor of 10 of the edge threshold in either
          # fit is a coin flip between two solvers and is not compared
          clear <- pmin(contrib_ref, contrib_fork) > 10*edge_tol | pmax(contrib_ref, contrib_fork) < edge_tol/10
          num_masked <- num_masked + sum(!clear)
          num_compared <- num_compared + length(clear)
          expect_identical(.path_net(ref, cell, l)[j, clear], (contrib_fork > edge_tol)[clear],
                           info = paste(info, cell))
        }
      }
    }
  }
  # the mask may not swallow the comparison: at most 2% of the regulator
  # positions sit in the band around the threshold
  expect_lt(num_masked/num_compared, 0.02)
})

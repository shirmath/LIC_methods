# Plan section 12, L15: plot_md_network() (R/plot_md_network.R), the
# igraph drawing of a network with one signed weight per edge. Data-only
# (L10): the drawing runs on a null device and the returned edge table is
# checked. The weight of edge k -> j is the slope of regulator k's
# contribution g_jk(x) = sum_m theta_jkm x^m at the cell's mean state, which
# is the coefficient itself at q = 1.

.edge_colours <- function() c(activating = "#1B7837", inhibiting = "#762A83")

## L15(a): q = 1. Weights are the coefficients, the sign gives the colour, the
## width is linear in |w| with the largest edge at 4 (or on a supplied
## `width_max`), the diagonal is off by default, `network` restricts the edge
## set, and the edges are listed in equation-major order.
test_that("L15(a): plot_md_network() weights at q = 1 are the coefficients", {
  .skip_unless_built("plot_md_network")
  skip_if_not_installed("igraph")
  coef_mat <- matrix(c(-1, 2, 0, 0, -1, -3, 0.5, 0, -1), 3, 3, byrow = TRUE)
  out <- .on_null_device(plot_md_network(coef_mat))
  expect_equal(names(out), c("edges", "graph"))
  edges <- out$edges
  expect_equal(names(edges), c("from", "to", "weight", "sign", "monotone", "lwd", "colour"))
  expect_equal(edges$to, c("1", "2", "3"))
  expect_equal(edges$from, c("2", "3", "1"))
  expect_equal(edges$weight, c(2, -3, 0.5))
  expect_equal(edges$sign, c(1, -1, 1))
  expect_true(all(edges$monotone))
  expect_equal(edges$colour, unname(.edge_colours()[c("activating", "inhibiting", "activating")]))
  expect_equal(edges$lwd, 0.5 + 3.5*abs(edges$weight)/3)
  expect_true(inherits(out$graph, "igraph"))
  expect_equal(igraph::ecount(out$graph), 3)
  expect_equal(igraph::vcount(out$graph), 3)

  out_diag <- .on_null_device(plot_md_network(coef_mat, bool_diagonal = TRUE))
  expect_equal(nrow(out_diag$edges), 6)
  expect_equal(sum(out_diag$edges$from == out_diag$edges$to), 3)

  out_width <- .on_null_device(plot_md_network(coef_mat, width_max = 6))
  expect_equal(out_width$edges$lwd, 0.5 + 3.5*abs(edges$weight)/6)

  keep <- matrix(FALSE, 3, 3); keep[2, 3] <- TRUE
  out_net <- .on_null_device(plot_md_network(coef_mat, network = keep))
  expect_equal(nrow(out_net$edges), 1)
  expect_equal(out_net$edges$weight, -3)

  out_named <- .on_null_device(plot_md_network(coef_mat, protein_names = c("a", "b", "c")))
  expect_equal(out_named$edges$from, c("b", "c", "a"))
  expect_equal(igraph::V(out_named$graph)$name, c("a", "b", "c"))

  empty <- .on_null_device(plot_md_network(diag(-1, 3)))
  expect_equal(nrow(empty$edges), 0)

  # arguments this function sets itself can still be overridden through `...`
  out_dots <- .on_null_device(plot_md_network(coef_mat, vertex.size = 20, edge.curved = 0))
  expect_equal(out_dots$edges, edges)
  expect_error(.on_null_device(plot_md_network(coef_mat, protein_names = c("a", "a", "b"))), "protein_names")
  coef_na <- coef_mat; coef_na[1, 2] <- NA
  expect_error(.on_null_device(plot_md_network(coef_na)), "NA")
})

## L15(b): q = 2. The weight is the slope g'(xbar_k) at the column mean of
## `X` (the Jacobian entry at the mean state), or at `x_ref` when given; a
## contribution whose derivative changes sign on the column's range is flagged
## non-monotone (drawn dashed); `X` is required.
test_that("L15(b): plot_md_network() weights at q > 1 are slopes at the mean state", {
  .skip_unless_built("plot_md_network")
  skip_if_not_installed("igraph")
  # equation 1 <- regulator 2: g(x) = x + x^2, x = (0, 0, 0, 2): mean 0.5, slope
  #   1 + 2*0.5 = 2 (the net slope over [0, 2] would be 3), monotone
  # equation 2 <- regulator 1: g(x) = x - x^2, x = (-1, 0, 0, 1): mean 0, slope
  #   1, but g' = 1 - 2x changes sign on [-1, 1], so non-monotone
  coef_mat <- matrix(0, 2, 4)
  coef_mat[1, 3:4] <- c(1, 1)
  coef_mat[2, 1:2] <- c(1, -1)
  X <- cbind(c(-1, 0, 0, 1), c(0, 0, 0, 2))
  out <- .on_null_device(plot_md_network(coef_mat, X = X))
  edges <- out$edges
  expect_equal(edges$to, c("1", "2"))
  expect_equal(edges$from, c("2", "1"))
  expect_equal(edges$weight, c(2, 1))
  expect_equal(edges$monotone, c(TRUE, FALSE))
  expect_equal(edges$sign, c(1, 1))
  expect_equal(edges$lwd, c(4, 0.5 + 3.5/2))
  # at a supplied state: slopes 1 + 2*0 = 1 and 1 - 2*1 = -1; the dashing
  # still follows the range of X
  out_ref <- .on_null_device(plot_md_network(coef_mat, X = X, x_ref = c(1, 0)))
  expect_equal(out_ref$edges$weight, c(1, -1))
  expect_equal(out_ref$edges$sign, c(1, -1))
  expect_equal(out_ref$edges$monotone, c(TRUE, FALSE))
  # the weights are the off-diagonal Jacobian entries at the mean state
  jac <- LICmethods:::.jacobian_at(coef_mat, 2, colMeans(X))
  expect_equal(edges$weight, jac[cbind(c(1, 2), c(2, 1))])
  expect_error(.on_null_device(plot_md_network(coef_mat)), "X")
  expect_error(.on_null_device(plot_md_network(coef_mat, X = X[, 1, drop = FALSE])), "X")
  expect_error(.on_null_device(plot_md_network(coef_mat, X = X, x_ref = 1)), "x_ref")
  expect_error(.on_null_device(plot_md_network(coef_mat, q = 3)), "q")
  X_na <- X; X_na[, 2] <- NA
  expect_error(.on_null_device(plot_md_network(coef_mat, X = X_na)), "column 2 of `X`")
})

## L15(c): a fitted path point and its truth give edge tables with the same
## columns, and the truth's weights at q = 1 are its coefficients.
test_that("L15(c): plot_md_network() on a prnode_fit and on the truth", {
  .skip_unless_built("plot_md_network")
  .skip_unless_built("fit_prnode", "public")
  skip_if_not_installed("igraph")
  ds <- .fixture_dataset("regression", num_d = 12)
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  fit <- fit_prnode(ds$public, lambda = lambda)
  out_fit <- .on_null_device(plot_md_network(fit$d_ode_est, network = fit$d_network_est))
  expect_equal(nrow(out_fit$edges), sum(fit$d_network_est[row(fit$d_network_est) != col(fit$d_network_est)]))
  out_truth <- .on_null_device(plot_md_network(ds$truth$d_ode_pars))
  off <- ds$truth$d_ode_pars[row(ds$truth$d_ode_pars) != col(ds$truth$d_ode_pars)]
  expect_equal(sort(out_truth$edges$weight), sort(off[off != 0]))
})

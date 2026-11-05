# Plan section 12, L13: plot_md_adjacency() (R/plot_md_adjacency.R),
# the adjacency image the vignette used to draw with a local helper. Data-only
# (L10): the drawing runs on a null device and the returned matrix is checked.

## L13: the matrix drawn is `network != 0`, rows = equation j, columns =
## regulator k, named by `protein_names`; a numeric and a logical input give
## the same drawing; non-square or non-matrix input errors.
test_that("L13: plot_md_adjacency() returns the logical matrix it draws", {
  .skip_unless_built("plot_md_adjacency")
  coef_mat <- matrix(c(-1, 2, 0, 0, -1, -3, 0.5, 0, -1), 3, 3, byrow = TRUE)
  drawn <- .on_null_device(plot_md_adjacency(coef_mat))
  expect_true(is.logical(drawn))
  expect_equal(unname(drawn), coef_mat != 0)
  expect_equal(rownames(drawn), c("1", "2", "3"))
  expect_equal(colnames(drawn), c("1", "2", "3"))
  drawn_logical <- .on_null_device(plot_md_adjacency(coef_mat != 0, main = "x"))
  expect_equal(drawn_logical, drawn)
  named <- .on_null_device(plot_md_adjacency(coef_mat, protein_names = c("a", "b", "c")))
  expect_equal(rownames(named), c("a", "b", "c"))
  expect_equal(colnames(named), c("a", "b", "c"))
  expect_error(.on_null_device(plot_md_adjacency(coef_mat[1:2, ])), "square")
  expect_error(.on_null_device(plot_md_adjacency(c(1, 0, 1))), "matrix")
  expect_error(.on_null_device(plot_md_adjacency(coef_mat, protein_names = c("a", "b"))), "protein_names")
})

## L13(b): a fitted network (one lambda of a prnode_fit) draws as its own
## logical matrix, so the vignette's panels are the package function's output.
test_that("L13(b): plot_md_adjacency() accepts a prnode_fit's network", {
  .skip_unless_built("plot_md_adjacency")
  .skip_unless_built("fit_prnode", "public")
  ds <- .fixture_dataset("regression", num_d = 12)
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  fit <- fit_prnode(ds$public, lambda = lambda)
  drawn <- .on_null_device(plot_md_adjacency(fit$d_network_est))
  expect_equal(unname(drawn), unname(fit$d_network_est))
})

## L13(c): `bool_diagonal = FALSE`
## blanks the self edges, so the image shows exactly the edge set that
## plot_md_network() draws by default and score_md_networks() scores.
test_that("L13(c): bool_diagonal = FALSE draws the off-diagonal edges that plot_md_network() draws", {
  .skip_unless_built("plot_md_adjacency")
  coef_mat <- matrix(c(-1, 2, 0, 0, -1, -3, 0.5, 0, -1), 3, 3, byrow = TRUE)
  drawn <- .on_null_device(plot_md_adjacency(coef_mat, bool_diagonal = FALSE))
  expected <- coef_mat != 0
  diag(expected) <- FALSE
  expect_equal(unname(drawn), expected)
  expect_equal(unname(.on_null_device(plot_md_adjacency(coef_mat, bool_diagonal = TRUE))), coef_mat != 0)
  expect_error(.on_null_device(plot_md_adjacency(coef_mat, bool_diagonal = NA)), "bool_diagonal")
  skip_if_not_installed("igraph")
  out <- .on_null_device(plot_md_network(coef_mat))
  edge_mat <- matrix(FALSE, 3, 3)
  edge_mat[cbind(match(out$edges$to, c("1", "2", "3")), match(out$edges$from, c("1", "2", "3")))] <- TRUE
  expect_equal(unname(drawn), edge_mat)
})

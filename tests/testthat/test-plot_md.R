# Plan section 12: plot_md_pair() and plot_md_pca3d() (R/plot_md_pair.R
# and R/plot_md_pca3d.R),
# L1-L12, data-only (L10: no pixel snapshots). Every drawing call runs on a
# null device and the returned data are what is checked.

# hcl.colors() lists a sequential palette dark to light, so the darkest entry
# is the first (the last is near white: #F9F9F9 for Blues 3)
.mother_blue <- function() grDevices::hcl.colors(3, "Blues 3")[1]
.daughter_red <- function() grDevices::hcl.colors(3, "Reds 3")[1]

.pca_cells <- function(ds, extra = 0, seed = 1){
  Y_list <- list(mother = ds$raw$Y_M, daughter = ds$raw$Y_D)
  t_list <- list(mother = ds$raw$m_obs_time, daughter = ds$raw$d_obs_time)
  set.seed(seed)
  for(i in seq_len(extra)){
    n_i <- sample(10:40, 1)
    Y_list[[2 + i]] <- matrix(stats::rnorm(n_i*ncol(ds$raw$Y_M)), n_i) + 0.5*ds$raw$Y_M[seq_len(n_i), ]
    t_list[[2 + i]] <- ds$raw$m_obs_time[seq_len(n_i)]
  }
  list(Y = Y_list, t = t_list)
}

.bind_panels <- function(panel_list) do.call(rbind, panel_list)

## L1: the returned data, with and without truth and fit.
test_that("L1: plot_md_pair() returns one panel per protein with observed, true and fitted values as documented", {
  .skip_unless_built("plot_md_pair")
  .skip_unless_built("fit_prnode", "public")
  ds <- .fixture_dataset("ode", num_d = 12)
  lambda <- .lambda_path(ds$public, num_lambda = 4)[3]
  fit <- fit_prnode(ds$public, lambda = lambda)
  panels <- .on_null_device(plot_md_pair(ds$raw, truth = ds$truth, fit = fit, times_e = ds$raw$times_e))
  expect_length(panels, 6)
  df <- .bind_panels(panels)
  expect_true(all(c("protein", "cell", "time", "observed", "true", "fitted", "colour") %in% names(df)))
  expect_equal(nrow(df), 6*(48 + 12))
  m_idx <- sapply(ds$raw$m_obs_time, function(tt) which(abs(ds$raw$times_e - tt) < 1e-10))
  pred_M <- stats::predict(fit, cell = "M")
  pred_D <- stats::predict(fit, cell = "D")
  for(j in 1:6){
    panel <- panels[[j]]
    m_rows <- panel$cell == "mother"
    d_rows <- panel$cell == "daughter"
    expect_equal(sum(m_rows), 48, info = j)
    expect_equal(sum(d_rows), 12, info = j)
    expect_equal(panel$observed[m_rows], unname(ds$raw$Y_M[, j]), info = j)
    expect_equal(panel$observed[d_rows], unname(ds$raw$Y_D[, j]), info = j)
    expect_equal(panel$time[m_rows], ds$raw$m_obs_time, info = j)
    expect_equal(panel$true[m_rows], unname(ds$truth$X_M[m_idx, j]), tolerance = 1e-12, info = j)
    expect_equal(panel$true[d_rows], unname(ds$truth$Y_D_noiseless[, j]), tolerance = 1e-12, info = j)
    expect_equal(panel$fitted[m_rows], unname(pred_M[, j]), tolerance = 1e-12, info = j)
    expect_equal(panel$fitted[d_rows], unname(pred_D[, j]), tolerance = 1e-12, info = j)
  }
  bare <- .on_null_device(plot_md_pair(ds$raw))
  df_bare <- .bind_panels(bare)
  expect_equal(nrow(df_bare), 6*60)
  expect_true(all(is.na(df_bare$true)))
  expect_true(all(is.na(df_bare$fitted)))
})

## L2: real data, no truth (the serialization row of the regimes table).
test_that("L2: plot_md_pair() draws a real pair with six panels and no error", {
  .skip_unless_built("plot_md_pair")
  raw <- .real_pair_raw(1)
  expect_gte(nrow(raw$Y_D), 3)
  panels <- .on_null_device(plot_md_pair(raw))
  expect_length(panels, 6)
  expect_equal(nrow(.bind_panels(panels)), 6*(nrow(raw$Y_M) + nrow(raw$Y_D)))
})

## L3: guards.
test_that("L3: mismatched proteins, a truth or fit of the wrong p, and a missing times_e each error by message", {
  .skip_unless_built("plot_md_pair")
  .skip_unless_built("fit_prnode", "public")
  ds <- .fixture_dataset("ode", num_d = 12)
  raw <- ds$raw
  bad_raw <- raw
  bad_raw$Y_D <- raw$Y_D[, 1:5]
  expect_error(.on_null_device(plot_md_pair(bad_raw)), "ncol|protein")
  wrong_truth <- ds$truth
  wrong_truth$X_M <- wrong_truth$X_M[, 1:5]
  expect_error(.on_null_device(plot_md_pair(raw, truth = wrong_truth, times_e = raw$times_e)), "truth|p")
  expect_error(.on_null_device(plot_md_pair(raw, truth = ds$truth)), "times_e")
  small <- .fixture_dataset("regression", num_d = 12)
  small_raw <- small$raw
  small_raw$Y_M <- small_raw$Y_M[, 1:3]
  small_raw$Y_D <- small_raw$Y_D[, 1:3]
  small_raw$protein_names <- small_raw$protein_names[1:3]
  fit6 <- fit_prnode(ds$public, lambda = .lambda_path(ds$public, num_lambda = 4)[3])
  expect_error(.on_null_device(plot_md_pair(small_raw, fit = fit6)), "fit|p")
})

## L4: the PCA is prcomp() on the stacked rows of all cells.
test_that("L4: scores equal prcomp() of the stacked matrix and the rotation is orthonormal", {
  .skip_unless_built("plot_md_pca3d")
  ds <- .fixture_dataset("ode", num_d = 12)
  cells <- .pca_cells(ds, extra = 3)
  out <- .on_null_device(plot_md_pca3d(cells$Y, cells$t, pair = c(mother = 1, daughter = 2), birth_time = ds$raw$birth_time))
  stacked <- do.call(rbind, cells$Y)
  pc <- stats::prcomp(stacked, center = TRUE, scale. = FALSE)
  scores_all <- pc$x[, 1:3]
  n_M <- nrow(cells$Y[[1]])
  n_D <- nrow(cells$Y[[2]])
  expect_equal(unname(out$scores[[1]]), unname(scores_all[seq_len(n_M), ]), tolerance = 1e-10)
  expect_equal(unname(out$scores[[2]]), unname(scores_all[n_M + seq_len(n_D), ]), tolerance = 1e-10)
  expect_equal(unname(out$rotation), unname(pc$rotation[, 1:3]), tolerance = 1e-10)
  expect_equal(crossprod(out$rotation), diag(3), tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(unname(out$sdev[1:3]), unname(pc$sdev[1:3]), tolerance = 1e-10)
})

## L5: the reference set matters.
test_that("L5: adding cells changes the drawn pair's projection; with the pair alone the scores are prcomp() of the pair", {
  .skip_unless_built("plot_md_pca3d")
  ds <- .fixture_dataset("ode", num_d = 12)
  two <- .pca_cells(ds, extra = 0)
  ten <- .pca_cells(ds, extra = 8)
  out_two <- .on_null_device(plot_md_pca3d(two$Y, two$t))
  out_ten <- .on_null_device(plot_md_pca3d(ten$Y, ten$t))
  expect_false(isTRUE(all.equal(out_two$projection[[1]], out_ten$projection[[1]])))
  pc <- stats::prcomp(rbind(ds$raw$Y_M, ds$raw$Y_D), center = TRUE, scale. = FALSE)
  expect_equal(unname(out_two$scores[[1]]), unname(pc$x[seq_len(48), 1:3]), tolerance = 1e-10)
})

## L6: the fixed isometric map.
test_that("L6: the projection is scores %*% P for the stated map, whose three axes are 120 degrees apart and equal", {
  .skip_unless_built("plot_md_pca3d")
  ds <- .fixture_dataset("ode", num_d = 12)
  cells <- .pca_cells(ds)
  out <- .on_null_device(plot_md_pca3d(cells$Y, cells$t))
  P <- rbind(c(cos(pi/6), -sin(pi/6)), c(-cos(pi/6), -sin(pi/6)), c(0, 1))
  for(i in 1:2) expect_equal(unname(out$projection[[i]]), unname(out$scores[[i]] %*% P), tolerance = 1e-10, info = i)
  axes <- diag(3) %*% P
  lengths <- sqrt(rowSums(axes^2))
  expect_equal(lengths, rep(lengths[1], 3), tolerance = 1e-12)
  angle <- function(u, v) acos(sum(u*v)/(sqrt(sum(u^2))*sqrt(sum(v^2))))
  expect_equal(c(angle(axes[1, ], axes[2, ]), angle(axes[2, ], axes[3, ]), angle(axes[1, ], axes[3, ])),
               rep(2*pi/3, 3), tolerance = 1e-12)
})

## L7: birth continuity (M9's visible consequence).
test_that("L7: on noiseless data the daughter's first projected point is the mother's at birth; on noisy data it is not", {
  .skip_unless_built("plot_md_pca3d")
  ds <- .fixture_dataset("regression", num_d = 12)
  birth_idx <- which(abs(ds$raw$m_obs_time - ds$raw$birth_time) < 1e-10)
  noiseless <- list(ds$truth$Y_M_noiseless, ds$truth$Y_D_noiseless)
  t_list <- list(ds$raw$m_obs_time, ds$raw$d_obs_time)
  out0 <- .on_null_device(plot_md_pca3d(noiseless, t_list, birth_time = ds$raw$birth_time))
  expect_identical(unname(out0$projection[[2]][1, ]), unname(out0$projection[[1]][birth_idx, ]))
  out1 <- .on_null_device(plot_md_pca3d(list(ds$raw$Y_M, ds$raw$Y_D), t_list, birth_time = ds$raw$birth_time))
  expect_false(isTRUE(all.equal(unname(out1$projection[[2]][1, ]), unname(out1$projection[[1]][birth_idx, ]))))
})

## L8/L11: the shade runs from light at the first frame to dark at the last.
test_that("L8: each cell's colours are non-increasing in lightness along time", {
  .skip_unless_built("plot_md_pca3d")
  ds <- .fixture_dataset("ode", num_d = 12)
  cells <- .pca_cells(ds)
  out <- .on_null_device(plot_md_pca3d(cells$Y, cells$t))
  for(i in 1:2){
    cols <- out$colours[[i]]
    expect_length(cols, nrow(cells$Y[[i]]))
    lightness <- colMeans(grDevices::col2rgb(cols))
    expect_true(all(diff(lightness) <= 1e-9), info = i)
    expect_lt(lightness[length(lightness)], lightness[1])
  }
})

## L9: guards.
test_that("L9: out-of-range pair indices, k != 3, a one-row cell and mismatched list lengths each error", {
  .skip_unless_built("plot_md_pca3d")
  ds <- .fixture_dataset("ode", num_d = 12)
  cells <- .pca_cells(ds, extra = 1)
  expect_error(.on_null_device(plot_md_pca3d(cells$Y, cells$t, pair = c(mother = 1, daughter = 9))), "pair")
  expect_error(.on_null_device(plot_md_pca3d(cells$Y, cells$t, k = 2)), "k")
  short <- cells
  short$Y[[2]] <- short$Y[[2]][1, , drop = FALSE]
  short$t[[2]] <- short$t[[2]][1]
  expect_error(.on_null_device(plot_md_pca3d(short$Y, short$t)), "row|frame")
  expect_error(.on_null_device(plot_md_pca3d(cells$Y, cells$t[1:2])), "obs_time_list|length")
})

## L12 (L11): the colour convention, defined once in .md_palette().
test_that("L12: the mother is always the darkest Blues 3 and the daughter the darkest Reds 3; the ramps are hcl.colors()", {
  .skip_unless_built("plot_md_pair")
  .skip_unless_built("plot_md_pca3d")
  ds <- .fixture_dataset("ode", num_d = 12)
  panels <- .on_null_device(plot_md_pair(ds$raw))
  df <- .bind_panels(panels)
  expect_true(all(df$colour[df$cell == "mother"] == .mother_blue()))
  expect_true(all(df$colour[df$cell == "daughter"] == .daughter_red()))
  cells <- .pca_cells(ds)
  out <- .on_null_device(plot_md_pca3d(cells$Y, cells$t))
  expect_true(all(out$colours[[1]] %in% grDevices::hcl.colors(48, "Blues 3")))
  expect_true(all(out$colours[[2]] %in% grDevices::hcl.colors(12, "Reds 3")))
  # the light end of each ramp is near white (saturation 0, hue 0 by
  # convention), so the hue is read on the saturated entries only
  hue <- function(cols){
    hsv_mat <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
    hsv_mat["h", hsv_mat["s", ] > 0.05]
  }
  mother_hue <- hue(c(.mother_blue(), out$colours[[1]]))
  daughter_hue <- hue(c(.daughter_red(), out$colours[[2]]))
  expect_gt(length(mother_hue), 40)
  expect_gt(length(daughter_hue), 8)
  expect_true(all(mother_hue >= 0.5 & mother_hue <= 0.75))
  expect_true(all(daughter_hue < 0.1 | daughter_hue > 0.9))
})

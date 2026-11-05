#' Draw a mother-daughter pair as polylines in the PCA space of many cells
#'
#' The rows of every cell in `Y_list` (one `p`-valued time point each) are
#' stacked and reduced by `stats::prcomp(center = TRUE, scale. = FALSE)`; the
#' first three components are the coordinates, so the reference set of cells
#' fixes the axes. Only the two cells named by `pair` are
#' drawn, each as a polyline through its frames in time order, the mother
#' through the blue ramp and the daughter through the red ramp of
#' `.md_palette()`, light at the first frame and dark at the last; the
#' daughter's birth is a large filled dot. The 3-D coordinates are projected
#' to the page by the fixed isometric map `x' = (x - y) cos 30`, `y' = z - (x
#' + y) sin 30` (the third axis up, the first two down-right and down-left,
#' pairwise 120 degrees apart), with the three axes drawn from the origin.
#'
#' @param Y_list a list of `n_i x p` matrices, one per cell.
#' @param obs_time_list a list of the cells' frame times, same length as
#' `Y_list`, each of length `n_i`.
#' @param pair indices into `Y_list` of the mother and the daughter, in that
#' order.
#' @param birth_time optional; drawn in the title.
#' @param k number of components; only 3 is supported.
#' @param ... further graphical parameters for `graphics::plot()`.
#' @returns invisibly, a list with `colours` (list of two per-point colour
#' vectors), `pair`, `projection` (list of two `n_i x 2` matrices), `rotation`
#' (`p x 3`), `scores` (list of two `n_i x 3` matrices) and `sdev`.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' out <- plot_md_pca3d(list(ds$raw$Y_M, ds$raw$Y_D),
#'                      list(ds$raw$m_obs_time, ds$raw$d_obs_time),
#'                      birth_time = ds$raw$birth_time)
#' round(out$sdev[1:3], 3)
#' @export
plot_md_pca3d <- function(Y_list,
                          obs_time_list,
                          pair = c(mother = 1, daughter = 2),
                          birth_time = NULL,
                          k = 3,
                          ...){
  if(k != 3) stop("`k` must be 3: only a three-component projection is drawn")
  stopifnot(is.list(Y_list), length(Y_list) >= 1)
  if(!is.list(obs_time_list) || length(obs_time_list) != length(Y_list)){
    stop("`obs_time_list` must be a list of the same length as `Y_list`")
  }
  Y_list <- lapply(Y_list, as.matrix)
  p <- ncol(Y_list[[1]])
  for(i in seq_along(Y_list)){
    if(ncol(Y_list[[i]]) != p) stop("every cell in `Y_list` must have p = ", p, " columns")
    if(nrow(Y_list[[i]]) != length(obs_time_list[[i]])){
      stop("cell ", i, " has ", nrow(Y_list[[i]]), " rows but ", length(obs_time_list[[i]]), " frame times")
    }
  }
  pair <- as.integer(pair)
  if(length(pair) != 2 || any(is.na(pair)) || any(pair < 1) || any(pair > length(Y_list))){
    stop("`pair` must give two indices into `Y_list` (1 to ", length(Y_list), ")")
  }
  for(i in pair){
    if(nrow(Y_list[[i]]) < 2) stop("cell ", i, " has fewer than 2 rows (frames); a polyline needs at least 2")
  }
  if(p < 3) stop("at least 3 proteins are needed for three components")

  stacked <- do.call(rbind, Y_list)
  pc <- stats::prcomp(stacked, center = TRUE, scale. = FALSE)
  scores_all <- pc$x[, 1:3, drop = FALSE]
  ends <- cumsum(vapply(Y_list, nrow, integer(1)))
  starts <- c(1, ends[-length(ends)] + 1)
  cell_scores <- lapply(pair, function(i) unname(scores_all[starts[i]:ends[i], , drop = FALSE]))
  # the isometric map: axes at 330, 210 and 90 degrees, 120 apart (plan L6);
  # with + sin(pi/6) in the first two rows the third axis would sit 60 degrees
  # from each of the others
  P <- rbind(c(cos(pi/6), -sin(pi/6)), c(-cos(pi/6), -sin(pi/6)), c(0, 1))
  projection <- lapply(cell_scores, function(sc) sc %*% P)
  n_M <- nrow(cell_scores[[1]])
  n_D <- nrow(cell_scores[[2]])
  palette <- .md_palette(n_M, n_D)
  colours <- list(palette$M, palette$D)

  axis_len <- max(abs(unlist(projection)))*0.5
  axes <- axis_len*diag(3) %*% P
  all_xy <- rbind(do.call(rbind, projection), axes, c(0, 0))
  graphics::plot(all_xy[, 1], all_xy[, 2], type = "n", xlab = "", ylab = "", axes = FALSE,
                 main = if(is.null(birth_time)) "mother (blue) and daughter (red) in PCA space" else
                   paste0("mother (blue) and daughter (red) in PCA space; birth at t = ", signif(birth_time, 3)),
                 asp = 1, ...)
  for(a in 1:3){
    graphics::arrows(0, 0, axes[a, 1], axes[a, 2], length = 0.08, col = "grey50")
    graphics::text(axes[a, 1], axes[a, 2], paste0("PC", a), pos = 3, col = "grey30", cex = 0.8)
  }
  for(i in 1:2){
    proj <- projection[[i]]
    n_i <- nrow(proj)
    graphics::segments(proj[-n_i, 1], proj[-n_i, 2], proj[-1, 1], proj[-1, 2], col = colours[[i]][-1], lwd = 2)
    graphics::points(proj[, 1], proj[, 2], col = colours[[i]], pch = 16, cex = 0.8)
  }
  graphics::points(projection[[2]][1, 1], projection[[2]][1, 2], col = colours[[2]][1], bg = colours[[2]][n_D],
                   pch = 21, cex = 2, lwd = 2)
  invisible(list(colours = colours,
                 pair = pair,
                 projection = projection,
                 rotation = pc$rotation[, 1:3, drop = FALSE],
                 scores = cell_scores,
                 sdev = pc$sdev))
}

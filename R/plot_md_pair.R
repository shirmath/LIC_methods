#' Plot one mother-daughter pair, one panel per protein
#'
#' Base graphics, a `ceiling(p / ncol) x ncol` grid of panels. Each panel
#' shows the observed values as points (mother blue, daughter red), the true
#' trajectories as solid lines when `truth` supplies `X_M` and `X_D` on
#' `times_e`, the fitted trajectories (`stats::predict()` of a `prnode_fit`,
#' joined by dashed lines) when `fit` is given, and a vertical line at the
#' daughter's birth. Works on real data (`truth = NULL`). The drawn data are
#' returned so that tests check numbers, not pixels.
#'
#' @param raw a dataset's `raw` list: `m_obs_time`, `Y_M`, `d_obs_time`, `Y_D`,
#' `birth_time` (optional), `protein_names` (optional).
#' @param truth `NULL` or a dataset's `truth` with `X_M`, `X_D` on `times_e`.
#' @param fit `NULL` or a `prnode_fit` at one `lambda`, fitted on this pair.
#' @param times_e the fine grid of `truth$X_M`; required when `truth` is given.
#' @param ncol number of panel columns.
#' @param protein_names panel titles; defaults to `raw$protein_names` or the
#' column names of `Y_M`.
#' @param ... further graphical parameters for `graphics::plot()`.
#' @returns invisibly, a list of `p` data frames with columns `protein`,
#' `cell`, `time`, `observed`, `true`, `fitted`, `colour`, one row per
#' observed frame of each cell.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' mu <- c(S = 1, M = 0.9, D = 0.9)
#' lambda_max <- lambda_max_prnode(ds$public, mu = mu, bool_pooled = TRUE)
#' lambda_path <- lambda_max*10^seq(0, -2, length.out = 10)
#' fit <- fit_prnode(ds$public, lambda = lambda_path, mu = mu)
#' panels <- plot_md_pair(ds$raw, truth = ds$truth, fit = fit[5],
#'                        times_e = ds$raw$times_e)
#' utils::head(panels[[1]])
#' @export
plot_md_pair <- function(raw,
                         truth = NULL,
                         fit = NULL,
                         times_e = NULL,
                         ncol = 2,
                         protein_names = NULL,
                         ...){
  stopifnot(is.list(raw), !is.null(raw$Y_M), !is.null(raw$Y_D))
  Y_M <- as.matrix(raw$Y_M)
  Y_D <- as.matrix(raw$Y_D)
  if(ncol(Y_M) != ncol(Y_D)){
    stop("`ncol(Y_M)` = ", ncol(Y_M), " and `ncol(Y_D)` = ", ncol(Y_D), " differ: the two cells must carry the same proteins")
  }
  p <- ncol(Y_M)
  m_obs_time <- raw$m_obs_time
  d_obs_time <- raw$d_obs_time
  stopifnot(length(m_obs_time) == nrow(Y_M), length(d_obs_time) == nrow(Y_D))
  birth_time <- if(!is.null(raw$birth_time)) raw$birth_time else d_obs_time[1]
  if(is.null(protein_names)) protein_names <- raw$protein_names
  if(is.null(protein_names)) protein_names <- colnames(Y_M)
  if(is.null(protein_names)) protein_names <- paste0("protein_", seq_len(p))

  true_M <- true_D <- NULL
  if(!is.null(truth)){
    if(is.null(times_e)) stop("`times_e` (the grid of `truth$X_M`) is required when `truth` is given")
    X_M <- as.matrix(truth$X_M)
    X_D <- as.matrix(truth$X_D)
    if(ncol(X_M) != p || ncol(X_D) != p){
      stop("`truth` has ", ncol(X_M), " proteins but the pair has p = ", p)
    }
    if(nrow(X_M) != length(times_e)) stop("`truth$X_M` must have one row per entry of `times_e`")
    m_idx <- .frame_index(m_obs_time, times_e, "m_obs_time")
    d_grid <- times_e[times_e >= birth_time - 1e-12]
    if(nrow(X_D) != length(d_grid)) stop("`truth$X_D` must be on the grid from the daughter's birth on")
    d_idx <- .frame_index(d_obs_time, d_grid, "d_obs_time")
    true_M <- X_M[m_idx, , drop = FALSE]
    true_D <- X_D[d_idx, , drop = FALSE]
  }
  fitted_M <- fitted_D <- NULL
  if(!is.null(fit)){
    if(!inherits(fit, "prnode_fit")) stop("`fit` must be a prnode_fit")
    if(length(fit$response_scale) != p || length(fit$lambda) != 1){
      stop("`fit` has p = ", length(fit$response_scale), " proteins and ", length(fit$lambda),
           " lambda value(s); one lambda and the pair's p = ", p, " are needed")
    }
    fitted_M <- stats::predict(fit, cell = "M")
    fitted_D <- stats::predict(fit, cell = "D")
    if(nrow(fitted_M) != nrow(Y_M) || nrow(fitted_D) != nrow(Y_D)) stop("`fit` was not fitted on this pair's frames")
  }

  palette <- .md_palette(3, 3)
  col_M <- palette$M[3]
  col_D <- palette$D[3]
  panel_list <- vector("list", p)
  old_par <- graphics::par(mfrow = c(ceiling(p/ncol), ncol), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(old_par))
  for(j in seq_len(p)){
    panel_df <- data.frame(protein = protein_names[j],
                           cell = c(rep("mother", nrow(Y_M)), rep("daughter", nrow(Y_D))),
                           time = c(m_obs_time, d_obs_time),
                           observed = c(Y_M[, j], Y_D[, j]),
                           true = if(is.null(true_M)) NA_real_ else c(true_M[, j], true_D[, j]),
                           fitted = if(is.null(fitted_M)) NA_real_ else c(fitted_M[, j], fitted_D[, j]),
                           colour = c(rep(col_M, nrow(Y_M)), rep(col_D, nrow(Y_D))),
                           stringsAsFactors = FALSE)
    ylim <- range(c(panel_df$observed, panel_df$true, panel_df$fitted,
                    if(!is.null(truth)) c(truth$X_M[, j], truth$X_D[, j])), na.rm = TRUE)
    graphics::plot(panel_df$time, panel_df$observed, col = panel_df$colour, pch = 16, cex = 0.7,
                   xlab = "time", ylab = protein_names[j], main = protein_names[j], ylim = ylim, ...)
    if(!is.null(truth)){
      graphics::lines(times_e, truth$X_M[, j], col = col_M, lwd = 1.5)
      graphics::lines(times_e[times_e >= birth_time - 1e-12], truth$X_D[, j], col = col_D, lwd = 1.5)
    }
    if(!is.null(fit)){
      graphics::lines(m_obs_time, fitted_M[, j], col = col_M, lty = 2, lwd = 1.5)
      graphics::lines(d_obs_time, fitted_D[, j], col = col_D, lty = 2, lwd = 1.5)
    }
    graphics::abline(v = birth_time, col = "grey50", lty = 3)
    panel_list[[j]] <- panel_df
  }
  invisible(panel_list)
}


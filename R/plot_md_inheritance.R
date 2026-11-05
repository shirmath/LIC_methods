#' Plot the flow-map similarity along the tuned path, with the pick
#'
#' `inheritance_flow(fit, bool_path = TRUE)` drawn against `lambda` on a log
#' axis: the plug-in similarity of the two fitted systems (standardized
#' units) at every point of the tuned path (1 at `lambda_max`, where both
#' systems are empty), the estimate (the plug-in at the pick) as a solid
#' horizontal line, and the pick of the tuning rule as a dotted vertical
#' line. When `truth` is given (a dataset's `truth`, which carries the
#' noiseless states), the truth's own similarity on the standardized systems
#' at the same horizon (`truth$inheritance_std`'s definition) is drawn as a
#' dashed line, so the path can be read against the value the daughter was
#' simulated at. The drawn numbers are
#' returned so that tests check numbers, not pixels.
#'
#' @param fit a `prnode_fit` from `fit_prnode()` with its penalty chosen by
#' the tuning rule, or a `prnode_tuning` (see `inheritance_flow()`).
#' @param h,x_ref passed to `inheritance_flow()`.
#' @param truth `NULL` or a dataset's `truth` (with `m_ode_pars`, `d_ode_pars`,
#' `q`, `Y_M_noiseless` and `Y_D_noiseless`).
#' @param ... further graphical parameters for `graphics::plot()`.
#' @returns invisibly, a list: `estimate`, `path_df` (the data frame of
#' `inheritance_flow()`: `index`, `lambda`, `similarity`, `bool_pick`),
#' `reference` (the truth's standardized similarity, or `NA` when `truth` is
#' `NULL`), `pick_idx` and `ylim` (the y-range drawn).
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' fit <- fit_prnode(ds$public)
#' out <- plot_md_inheritance(fit, truth = ds$truth)
#' c(estimate = out$estimate, truth = out$reference)
#' @export
plot_md_inheritance <- function(fit,
                                h = NULL,
                                truth = NULL,
                                x_ref = NULL,
                                ...){
  flow <- inheritance_flow(fit, bool_path = TRUE, h = h, x_ref = x_ref)
  reference <- NA_real_
  if(!is.null(truth)){
    if(!is.list(truth) || is.null(truth$m_ode_pars) || is.null(truth$d_ode_pars)){
      stop("`truth` must be a dataset's truth, with `m_ode_pars` and `d_ode_pars`")
    }
    if(is.null(truth$Y_M_noiseless) || is.null(truth$Y_D_noiseless)){
      stop("`truth` must carry `Y_M_noiseless` and `Y_D_noiseless` (the standardization's pooled sds)")
    }
    q_truth <- if(!is.null(truth$q)) truth$q else ncol(as.matrix(truth$m_ode_pars))/nrow(as.matrix(truth$m_ode_pars))
    reference <- .flow_map_similarity_std(as.matrix(truth$m_ode_pars), as.matrix(truth$d_ode_pars),
                                          truth$Y_M_noiseless, truth$Y_D_noiseless, h = flow$h, q = q_truth)
  }
  path_df <- flow$path_df
  ylim <- range(c(path_df$similarity, flow$estimate, reference), na.rm = TRUE)
  graphics::plot(path_df$lambda, path_df$similarity, log = "x", type = "n", xlab = "lambda",
                 ylab = "flow-map similarity of the fitted systems", main = "inheritance along the tuned path",
                 ylim = ylim, ...)
  graphics::lines(path_df$lambda, path_df$similarity, type = "b", pch = 16)
  graphics::abline(h = flow$estimate, lty = 1)
  if(!is.na(reference)) graphics::abline(h = reference, lty = 2)
  graphics::abline(v = path_df$lambda[flow$pick_idx], lty = 3)
  invisible(list(estimate = flow$estimate,
                 path_df = path_df,
                 reference = reference,
                 pick_idx = flow$pick_idx,
                 ylim = ylim))
}

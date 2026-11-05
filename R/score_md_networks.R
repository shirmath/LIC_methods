#' Score estimated mother-daughter networks against a truth, by edge type
#'
#' Confusion counts, MCC and F1 for the scaffold, mother and daughter
#' networks, broken down by the truth's edge type, plus two rates: the **leak rate** (the share of the other cell's
#' private edges that appear in this cell: `daughter_only` edges in the
#' mother's network, `mother_only` edges in the daughter's) and the **borrow
#' rate** (the share of shared edges the cell recovers). The scalar MCC of the
#' daughter alone hides a leak. The diagonal (self-degradation,
#' `edge_type == "self"`) is excluded unless `bool_score_diagonal = TRUE`.
#' MCC and F1 are 0 when their denominator vanishes (the usual convention); a rate with no
#' qualifying edge is `NA`.
#'
#' @section The two rates:
#' Write `A_hat_c[j, k]` for the estimated adjacency of cell `c` (1 when the
#' estimate has the edge `j <- k`) and partition the scored (off-diagonal)
#' positions by the truth's `edge_type`: `E_shared` (`shared_equal` or
#' `shared_unequal`, in both true networks), `E_M` (`mother_only`) and `E_D`
#' (`daughter_only`). Then
#'
#' `borrow_rate[c] = mean over (j, k) in E_shared of A_hat_c[j, k]`,
#'
#' `leak_rate["D"] = mean over (j, k) in E_M of A_hat_D[j, k]`,
#' `leak_rate["M"] = mean over (j, k) in E_D of A_hat_M[j, k]`.
#'
#' The borrow rate is the cell's true-positive rate on the shared edges
#' (whatever the fit's degree of borrowing; a fit with no borrowing still has
#' one), and the leak rate is its false-positive rate on the edges that
#' belong to the other cell alone. Both partition the positions by the
#' **true** edge types, so they exist only on simulated data; on real data
#' there is no private edge to leak and no shared edge to recover, and the
#' nearest observable quantities are the discount `fit$discount` (a property
#' of the design and `mu`, which needs no truth) and the overlap of the two
#' estimated networks, which is not a leak rate.
#'
#' @param fit any list with `m_network_est`, `d_network_est`, `s_network_est`
#' (`p x p`, logical or 0/1; a `prnode_fit` at one `lambda`).
#' @param truth a `sim_md_network()` list or a dataset's `truth` (`m_network`,
#' `d_network`, `s_network`, `edge_type`).
#' @param bool_score_diagonal include the diagonal positions.
#' @returns a list with `borrow_rate` (named `M`, `D`), `confusion` (list `S`,
#' `M`, `D` of matrices with one row per edge type and columns `TP`, `FP`,
#' `FN`, `TN`), `f1` and `mcc` (named `S`, `M`, `D`) and `leak_rate` (named
#' `M`, `D`).
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' mu <- c(S = 1, M = 0.9, D = 0.9)
#' lambda_max <- lambda_max_prnode(ds$public, mu = mu, bool_pooled = TRUE)
#' lambda_path <- lambda_max*10^seq(0, -2, length.out = 10)
#' fit <- fit_prnode(ds$public, lambda = lambda_path, mu = mu)
#' scores <- score_md_networks(fit[5], ds$truth)
#' scores$mcc
#' scores$confusion$D
#' # the daughter's MCC along the path
#' sapply(seq_along(lambda_path),
#'        function(l) score_md_networks(fit[l], ds$truth)$mcc[["D"]])
#' @export
score_md_networks <- function(fit,
                              truth,
                              bool_score_diagonal = FALSE){
  for(nm in c("m_network_est", "d_network_est", "s_network_est")){
    if(is.null(fit[[nm]])) stop("`fit` must carry `", nm, "`")
    if(length(dim(fit[[nm]])) == 3) stop("`fit$", nm, "` has a lambda dimension; score one lambda at a time (fit[l])")
  }
  for(nm in c("m_network", "d_network", "s_network", "edge_type")){
    if(is.null(truth[[nm]])) stop("`truth` must carry `", nm, "`")
  }
  as_adj <- function(mat) matrix(as.logical(as.matrix(mat) != 0), nrow(mat), ncol(mat))
  est <- list(S = as_adj(fit$s_network_est), M = as_adj(fit$m_network_est), D = as_adj(fit$d_network_est))
  tru <- list(S = as_adj(truth$s_network), M = as_adj(truth$m_network), D = as_adj(truth$d_network))
  p <- nrow(tru$M)
  edge_type <- truth$edge_type
  stopifnot(all(dim(edge_type) == c(p, p)), all(vapply(est, function(m) all(dim(m) == c(p, p)), logical(1))))
  scored <- if(bool_score_diagonal) matrix(TRUE, p, p) else row(edge_type) != col(edge_type)
  type_levels <- c("absent", "shared_equal", "shared_unequal", "mother_only", "daughter_only", "self")
  type_levels <- intersect(type_levels, unique(c(edge_type[scored], type_levels[1:5])))

  confusion <- list()
  mcc <- f1 <- c(S = NA_real_, M = NA_real_, D = NA_real_)
  for(cell in c("S", "M", "D")){
    e <- est[[cell]][scored]
    t <- tru[[cell]][scored]
    types <- edge_type[scored]
    conf <- t(vapply(type_levels, function(ty){
      sel <- types == ty
      c(TP = sum(e[sel] & t[sel]), FP = sum(e[sel] & !t[sel]), FN = sum(!e[sel] & t[sel]),
        TN = sum(!e[sel] & !t[sel]))
    }, numeric(4)))
    rownames(conf) <- type_levels
    confusion[[cell]] <- conf
    tp <- sum(conf[, "TP"])
    fp <- sum(conf[, "FP"])
    fn <- sum(conf[, "FN"])
    tn <- sum(conf[, "TN"])
    denom <- sqrt(as.numeric(tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
    mcc[[cell]] <- if(denom == 0) 0 else (tp*tn - fp*fn)/denom
    f1[[cell]] <- if(2*tp + fp + fn == 0) 0 else 2*tp/(2*tp + fp + fn)
  }
  rate <- function(est_mat, types_wanted){
    sel <- scored & (edge_type %in% types_wanted)
    if(sum(sel) == 0) return(NA_real_)
    mean(est_mat[sel])
  }
  shared_types <- c("shared_equal", "shared_unequal")
  list(borrow_rate = c(M = rate(est$M, shared_types), D = rate(est$D, shared_types)),
       confusion = confusion,
       f1 = f1,
       leak_rate = c(M = rate(est$M, "daughter_only"), D = rate(est$D, "mother_only")),
       mcc = mcc)
}

#' Choose the penalty of a PRNODE fit
#'
#' The tuning rule: the
#' penalty `lambda` is chosen from the data along a path at a fixed multiplier
#' triplet `mu` by the known-variance BIC, the minimum of
#' `sum_j [ RSS_jl + log(n) k_jl ]` over the path. The multipliers are
#' not tuned: every per-dataset rule for `(mu_M, mu_D)` did worse in
#' simulations than the fixed pair, so they are a design constant, held at the
#' package default.
#'
#' @details
#' Path. Unless `path` or `lambda_path` is given, the path is `num_lambda`
#' values log-spaced over `num_decades` decades below the pooled
#' `lambda_max_prnode()`, fitted by `fit_prnode_path()` at `mu`. Every protein
#' shares the one path. Only the supports along the path are used below; every
#' coefficient is refit by least squares on its support, with a column aliased
#' within the support (the scaffold column when both deviation columns are in,
#' since `x_S = x_M + x_D`) set to zero.
#'
#' Criterion. With `n = n_M + n_D`, `RSS_jl` the refit's residual sum of
#' squares for protein `j` at path point `l` and `k_jl` the number of active
#' groups (a group is the `q` columns of one regulator in one block; it is
#' active when any of its coefficients is nonzero, and the refit uses all its
#' columns), the criterion is `sum_j [ RSS_jl + log(n) k_jl ]` and the pick is
#' its (first, sparsest) minimum. The design already weights each row by the
#' inverse noise variance (one common weight for the two cells, see
#' `make_prnode_design()`), so the noise variance is 1 and the Gaussian
#' criterion with known variance applies. At `q = 1` a group is one column and
#' `k_jl` is the number of nonzero coefficients.
#'
#' The smoothing level `spar` of the design is not chosen here: `public` is
#' one design. `fit_prnode()` on raw observations chooses it, by minimizing
#' this criterion jointly over `spar` and `lambda`.
#'
#' Checks and warnings. `public$x` must be a base numeric matrix (a sparse or
#' dense `Matrix` object is refused) and `public$x` and `public$y` must be
#' finite: a non-finite entry stops before any path is fitted, naming the
#' protein (in `y`) or the regulator and block (in `x`). Non-convergence of
#' `fit_prnode_path()` at some (protein, `lambda`) passes through as the
#' solver's own warning, with `path$converged` recording where; when the pick
#' itself sits on a non-converged fit, a second warning names the path index
#' and the protein(s), since the support there may be inflated.
#'
#' Cost. The path is the dominant cost: on a `p = 6` dataset with 48 mother
#' frames the reference solver takes 1 to 3 seconds for the 80-point path,
#' against a fraction of a second for the criterion. A path fitted elsewhere
#' (`fit_prnode_path()`, or an equivalent solver at `q = 1`) can be passed as
#' `path`.
#'
#' @param public the design list of `make_prnode_design()` (or a dataset's
#' `public`), with a base numeric `x` and finite `x` and `y`.
#' @param edge_tol,max_iter,tol passed to `fit_prnode_path()` when the path is
#' fitted here.
#' @param lambda_path `NULL` (the default grid), or a decreasing positive
#' vector to fit the path on; ignored when `path` is given.
#' @param mu multipliers `c(S, M, D)` of the scaffold and deviation penalties,
#' held fixed; with `path`, it must equal `path$mu` (or be left at its
#' default, in which case it is read from the path).
#' @param num_decades decades below the pooled `lambda_max` the default grid
#' spans.
#' @param num_lambda number of points of the default grid.
#' @param path `NULL`, or a `prnode_path` from `fit_prnode_path()` on this
#' `public`; then the path is not refitted and `lambda_path`, `num_decades`,
#' `num_lambda`, `edge_tol`, `max_iter` and `tol` are ignored.
#' @param verbose a non-negative number; `0` prints nothing, `1` reports the
#' phases (path, criterion), higher levels are passed to `fit_prnode_path()`
#' one lower.
#' @returns an object of class `prnode_tuning`: `criterion` (length `L`, the
#' known-variance BIC along the path), `k` (length `L`, active groups summed
#' over proteins), `lambda` (the chosen value, a scalar), `lambda_path`
#' (length `L`, decreasing), `mu` (named `S`, `M`, `D`), `path` (the
#' `prnode_path`), `pick_idx` (index of `lambda` on the path), `rss` (length
#' `L`, summed over proteins) and `spar` (the design's smoothing level, `NA`
#' when the design does not record one).
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' # a short path keeps the example quick; the default is 80 points over 4 decades
#' tuning <- tune_prnode(ds$public, num_lambda = 20, num_decades = 2)
#' print(tuning)
#' fit <- fit_prnode(ds$public, lambda = tuning)
#' fit$d_network_est
#' @export
tune_prnode <- function(public,
                        edge_tol = 1e-8,
                        lambda_path = NULL,
                        max_iter = 1e4,
                        mu = c(S = 1, M = 1, D = 1),
                        num_decades = 4,
                        num_lambda = 80,
                        path = NULL,
                        tol = 1e-10,
                        verbose = 0){
  if(!is.list(public) || is.null(public$x) || is.null(public$y) || is.null(public$n_M) || is.null(public$n_D)){
    stop("`public` must be the design list of make_prnode_design()")
  }
  if(length(verbose) != 1 || !is.numeric(verbose) || is.na(verbose) || verbose < 0){
    stop("`verbose` must be a single non-negative number (0 = silent, 1 = the phases, 2 and up = the solver's levels)")
  }
  n_M <- public$n_M
  n_D <- public$n_D
  n <- n_M + n_D
  x_mat <- public$x
  y_mat <- as.matrix(public$y)
  if(nrow(x_mat) != n || nrow(y_mat) != n){
    stop("`public$x` and `public$y` must have `n_M + n_D` = ", n, " rows (the mother's frames then the daughter's)")
  }
  p <- ncol(y_mat)
  # q from the design (its `q`, or its `index`), not from `y`, which may hold a subset of the proteins
  q <- .design_q(public, pq = ncol(x_mat)/3, p_y = p)
  .check_design_entries(x_mat = x_mat,
                        y_mat = y_mat,
                        q = q,
                        protein_names = public$protein_names,
                        arg = "public")

  # Step 1: the path (fitted here, or taken as given)
  if(is.null(path)){
    mu <- .mu_vec(mu)
    if(is.null(lambda_path)){
      if(length(num_lambda) != 1 || num_lambda < 2 || num_lambda != round(num_lambda)){
        stop("`num_lambda` must be a whole number of at least 2")
      }
      if(length(num_decades) != 1 || !is.finite(num_decades) || num_decades <= 0){
        stop("`num_decades` must be a positive number")
      }
      lambda_max <- lambda_max_prnode(public, mu = mu, bool_pooled = TRUE)
      lambda_path <- lambda_max*10^seq(0, -num_decades, length.out = num_lambda)
    }
    if(verbose > 0){
      print(paste0("tune_prnode: fitting the path of ", length(lambda_path), " lambda values for ", p, " protein",
                   if(p != 1) "s"))
    }
    path <- fit_prnode_path(public, lambda = lambda_path, edge_tol = edge_tol, max_iter = max_iter, mu = mu,
                            tol = tol, verbose = max(verbose - 1, 0))
  } else {
    if(!inherits(path, "prnode_path")) stop("`path` must be a prnode_path from fit_prnode_path()")
    if(!missing(mu) && any(abs(.mu_vec(mu) - path$mu) > 1e-12)){
      stop("`path` was fitted at mu = (", paste(signif(path$mu, 3), collapse = ", "), "), not the `mu` given; ",
           "pass the path's multipliers or refit the path")
    }
    mu <- .mu_vec(path$mu)
    if(dim(path$s)[1] != p || 3*dim(path$s)[2] != ncol(x_mat)){
      stop("`path` does not match `public`: it has ", dim(path$s)[1], " proteins and ", 3*dim(path$s)[2],
           " design columns against ", p, " and ", ncol(x_mat))
    }
  }
  lambda_path <- path$lambda
  L <- length(lambda_path)
  # the supports, each expanded to whole groups (the solver zeroes a group whole, so this changes nothing unless
  # a coordinate of an active group is exactly zero; it makes the refit "all columns of the active groups")
  support_list <- lapply(seq_len(p), function(j) .expand_groups(.path_coef_mat(path, j) != 0, q = q))

  # Step 2: the known-variance BIC and its minimum
  if(verbose > 0) print("tune_prnode: computing the known-variance BIC along the path")
  ic <- .known_variance_bic(x_mat = x_mat, y_mat = y_mat, support_list = support_list, q = q)
  pick_idx <- which.min(ic$criterion)
  if(verbose > 0){
    print(paste0("tune_prnode: pick at path index ", pick_idx, " of ", L,
                 " (lambda = ", signif(lambda_path[pick_idx], 3), ")"))
  }
  .warn_if_unconverged(converged_mat = path$converged, pick_idx = pick_idx)

  stopifnot(length(ic$criterion) == L, pick_idx >= 1, pick_idx <= L)
  structure(list(criterion = ic$criterion,
                 k = ic$k,
                 lambda = lambda_path[pick_idx],
                 lambda_path = lambda_path,
                 mu = mu,
                 path = path,
                 pick_idx = pick_idx,
                 rss = ic$rss,
                 spar = if(is.null(public$spar)) NA_real_ else public$spar),
            class = "prnode_tuning")
}

#' Print a prnode_tuning
#'
#' Reports the path, the pick and its number of active groups (summed over
#' proteins; at `q = 1` the nonzero coefficients), and the smoothing level of
#' the design when it records one.
#'
#' @param x a `prnode_tuning`.
#' @param ... ignored.
#' @returns `x`, invisibly.
#' @export
print.prnode_tuning <- function(x, ...){
  L <- length(x$lambda_path)
  cat("PRNODE tuning (class prnode_tuning)\n")
  cat("  path of", L, "lambda values from", signif(x$lambda_path[1], 3), "to", signif(x$lambda_path[L], 3),
      "at mu (S, M, D) =", paste(signif(x$mu, 3), collapse = ", "), "\n")
  cat("  pick (known-variance BIC): index", x$pick_idx, " lambda =", signif(x$lambda, 3),
      " active groups =", x$k[x$pick_idx], "\n")
  if(!is.na(x$spar)) cat("  smoothing level spar =", x$spar, "\n")
  invisible(x)
}

#' The coefficient matrix of one protein along a path
#'
#' Stacks the scaffold, mother-deviation and daughter-deviation coefficients
#' of protein `j` in the order of the design's columns (`index = rep(1:(3p),
#' each = q)`: the `p` scaffold groups, then the `p` mother groups, then the
#' `p` daughter groups).
#'
#' @param path a `prnode_path`.
#' @param j the protein (row of the path's arrays).
#' @returns a `3pq x L` matrix.
#' @noRd
.path_coef_mat <- function(path,
                           j){
  pq <- dim(path$s)[2]
  stack <- function(arr) matrix(arr[j, , , drop = FALSE], nrow = pq)
  rbind(stack(path$s), stack(path$d_M), stack(path$d_D))
}

#' The group size `q` of a design
#'
#' Read from the design's `q`, else from its `index` (`rep(seq_len(3p), each
#' = q)`, so `3p` distinct values), else assumed from the number of response
#' columns; never from `y` first, since `fit_prnode_path()` accepts a `y`
#' holding a subset of the proteins (test plan F2) while `x` keeps every
#' regulator.
#'
#' @param design the design list, or a path's `design` (the list without `x`
#' and `y`).
#' @param pq the number of columns per block (`ncol(x) / 3`).
#' @param p_y the number of response columns, the last resort.
#' @returns a whole number dividing `pq`.
#' @noRd
.design_q <- function(design,
                      pq,
                      p_y){
  q <- design$q
  if(is.null(q) && !is.null(design$index)) q <- pq/(length(unique(design$index))/3)
  if(is.null(q)) q <- pq/p_y
  if(length(q) != 1 || !is.numeric(q) || !is.finite(q) || q < 1 || q != round(q) || pq %% q != 0){
    stop("cannot read the group size `q` of the design: `q` = ", paste(q, collapse = ", "), " against ", pq,
         " columns per block")
  }
  as.integer(q)
}

#' Expand a support to whole groups
#'
#' @param support_mat `3pq x L` logical.
#' @param q rows per group.
#' @returns the `3pq x L` logical with every row of an active group `TRUE`.
#' @noRd
.expand_groups <- function(support_mat,
                           q){
  stopifnot(length(q) == 1, q >= 1, q == round(q), nrow(support_mat) %% q == 0)
  group_id <- rep(seq_len(nrow(support_mat)/q), each = q)
  active <- rowsum(support_mat + 0L, group_id, reorder = FALSE) > 0
  out <- active[group_id, , drop = FALSE]
  dimnames(out) <- NULL
  out
}

#' The number of active groups at each path point
#'
#' A group is `q` consecutive rows of the support (one regulator in one
#' block, in the order of `.path_coef_mat()`); it is active when any of its
#' rows is `TRUE`. The one definition of "support size" the criterion and
#' `print.prnode_tuning()` share.
#'
#' @param support_mat `3pq x L` logical.
#' @param q columns per group; `nrow(support_mat)` must be a multiple of it.
#' @returns a length-`L` numeric vector of counts in `0:(3p)`.
#' @noRd
.active_group_count <- function(support_mat,
                                q){
  stopifnot(length(q) == 1, q >= 1, q == round(q), nrow(support_mat) %% q == 0)
  group_id <- rep(seq_len(nrow(support_mat)/q), each = q)
  colSums(rowsum(support_mat + 0L, group_id, reorder = FALSE) > 0)
}

#' Warn when the pick sits on a non-converged fit
#'
#' The solver's own warning already says that some (protein, `lambda`) hit
#' `max_iter`; this second one fires only when the point the rule returns is
#' among them, since its support may then be inflated.
#'
#' @param converged_mat the path's `p x L` logical `converged` (or `NULL`,
#' which is silent).
#' @returns `NULL`, invisibly; called for its warning.
#' @noRd
.warn_if_unconverged <- function(converged_mat,
                                 pick_idx){
  if(is.null(converged_mat)) return(invisible(NULL))
  bad <- which(!converged_mat[, pick_idx])
  if(length(bad) > 0){
    warning("tune_prnode: the pick (path index ", pick_idx, ") did not converge for protein",
            if(length(bad) > 1) "s " else " ", paste(bad, collapse = ", "),
            "; the support there may be inflated (see `path$converged`; raise `max_iter` or loosen `tol`)")
  }
  invisible(NULL)
}

#' Least-squares refit of every support along a path
#'
#' Each distinct support is refit once; a column aliased within its support
#' (an `NA` coefficient from `lm.fit()`'s pivoting) is set to zero, and a
#' column outside the support is zero.
#'
#' @param x_mat `n x 3pq` design (or its training rows).
#' @param y_vec the response column (same rows).
#' @param support_mat `3pq x L` logical, `TRUE` where the path coefficient is nonzero.
#' @returns a `3pq x L` coefficient matrix.
#' @noRd
.refit_support <- function(x_mat,
                           y_vec,
                           support_mat){
  key_vec <- apply(support_mat, 2, paste, collapse = "")
  first_idx <- which(!duplicated(key_vec))
  beta_mat <- matrix(0, nrow(support_mat), length(first_idx))
  for(i in seq_along(first_idx)){
    active_idx <- which(support_mat[, first_idx[i]])
    if(length(active_idx) > 0){
      coef_vec <- stats::lm.fit(x_mat[, active_idx, drop = FALSE], y_vec)$coefficients
      coef_vec[is.na(coef_vec)] <- 0
      beta_mat[active_idx, i] <- coef_vec
    }
  }
  beta_mat[, match(key_vec, key_vec[first_idx]), drop = FALSE]
}

#' The known-variance BIC along a path
#'
#' `sum_j [ RSS_jl + log(n) k_jl ]`, with `RSS_jl` the least-squares refit on
#' the path's support (all columns of the active groups) and `k_jl` the number
#' of active groups (`.active_group_count()`; the number of nonzero columns at
#' `q = 1`). The design's weights put the noise at unit variance (up to the
#' estimate of the common weight), which is why the Gaussian criterion with a
#' known variance applies and the profile-likelihood BIC, whose variance
#' estimate every added column shrinks, does not.
#'
#' @param support_list per protein, the `3pq x L` expanded support.
#' @param q columns per group of each support (`nrow(support) = 3pq`).
#' @returns list `criterion`, `k` and `rss` (length `L`; `k` and `rss` summed
#' over proteins).
#' @noRd
.known_variance_bic <- function(x_mat,
                                y_mat,
                                support_list,
                                q = 1){
  n <- nrow(x_mat)
  L <- ncol(support_list[[1]])
  k_sum <- rss_sum <- numeric(L)
  for(j in seq_along(support_list)){
    beta_mat <- .refit_support(x_mat, y_mat[, j], support_list[[j]])
    # a zero RSS cannot happen on noisy data; guard all the same
    rss_sum <- rss_sum + pmax(colSums((y_mat[, j] - x_mat %*% beta_mat)^2), .Machine$double.eps)
    k_sum <- k_sum + .active_group_count(support_list[[j]], q = q)
  }
  list(criterion = unname(rss_sum + log(n)*k_sum),
       k = unname(k_sum),
       rss = unname(rss_sum))
}

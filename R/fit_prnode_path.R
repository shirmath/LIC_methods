#' Estimate each protein's observation-noise standard deviation
#'
#' Per protein, `"diff2"` is the standard deviation of the series' second
#' differences divided by `sqrt(6)` (the sd of a second difference of iid
#' noise with unit variance; the numerator of `noise_stat()`), and `"spline"`
#' is the residual standard deviation of `stats::smooth.spline(all.knots =
#' TRUE, spar)`, `sqrt(RSS / (n - df))`. Both are equivariant to rescaling
#' the series and invariant to shifting it or the time axis.
#' `"diff2"` needs no smoother but is biased upward when the signal curves
#' within two frames (on a six-cycle sinusoid over 48 frames the bias is
#' large); it is the provisional default of the loss weights.
#'
#' @param obs_time the `n` frames, increasing; at least 4.
#' @param Y `n x p` matrix of observations, columns are proteins.
#' @param method `"diff2"` or `"spline"`.
#' @param spar smoothing parameter for `method = "spline"`.
#' @returns an unnamed numeric vector of length `p`.
#' @examples
#' obs_time <- seq(0, 1, length.out = 48)
#' set.seed(1)
#' Y <- cbind(sin(2*pi*obs_time), cos(2*pi*obs_time)) +
#'   matrix(stats::rnorm(96, sd = 0.1), 48)
#' estimate_noise_sd(obs_time, Y)                       # near 0.1
#' estimate_noise_sd(obs_time, Y, method = "spline")
#' @export
estimate_noise_sd <- function(obs_time,
                              Y,
                              method = c("diff2", "spline"),
                              spar = 0.5){
  method <- match.arg(method)
  Y <- as.matrix(Y)
  n <- nrow(Y)
  if(n < 4){
    stop("`Y` has ", n, " frames; estimating the noise sd needs at least 4 (two second ",
         "differences, or four points for a smoothing spline)")
  }
  if(length(obs_time) != n) stop("`obs_time` must have one entry per row of `Y`")
  if(method == "diff2"){
    sigma_vec <- apply(Y, 2, function(y) stats::sd(diff(y, differences = 2))/sqrt(6))
  } else {
    sigma_vec <- apply(Y, 2, function(y){
      fit <- stats::smooth.spline(obs_time, y, all.knots = TRUE, spar = spar)
      resid_vec <- y - stats::fitted(fit)
      sqrt(sum(resid_vec^2)/max(n - fit$df, 1))
    })
  }
  unname(as.numeric(sigma_vec))
}

#' The smallest penalty at which every group is zero
#'
#' `max_g ||Q_g' y_j|| / (sqrt(q) mu_g)` per protein, over the orthonormalized
#' groups of the design: the KKT bound at
#' `b = 0`, including the multipliers.
#'
#' @param design the `public` list of `make_prnode_design()`.
#' @param mu multipliers `c(S, M, D)` of the scaffold and the two deviation
#' penalties; names optional, in that order.
#' @param bool_pooled if `TRUE`, the maximum over proteins (one scalar).
#' @returns an unnamed numeric vector of length `p`, or a scalar.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' lambda_max_prnode(ds$public)                       # one value per protein
#' lambda_max_prnode(ds$public, bool_pooled = TRUE)   # their maximum
#' @export
lambda_max_prnode <- function(design,
                              mu = c(S = 1, M = 1, D = 1),
                              bool_pooled = FALSE){
  mu <- .mu_vec(mu)
  blocks <- .design_blocks(design, mu)
  lambda_vec <- vapply(seq_len(ncol(design$y)), function(j){
    max(vapply(blocks$groups, function(gr){
      sqrt(sum(crossprod(gr$Q, design$y[, j])^2))/gr$w
    }, numeric(1)))
  }, numeric(1))
  if(bool_pooled) max(lambda_vec) else unname(lambda_vec)
}

#' Fit the PRNODE group lasso along a penalty path
#'
#' The reference solver of the PRNODE estimator. Per protein, each of the
#' `3p` column groups of the design is orthonormalized by a thin QR, and the
#' problem `min_b 1/2 ||y - sum_g Q_g b_g||^2 + lambda sum_g w_g ||b_g||` with
#' `w_g = sqrt(q) mu_g` is solved by cyclic group soft-thresholding: full sweeps
#' over the groups alternate with sweeps over the active set, with warm starts
#' along the decreasing `lambda` path. Once the active set has been stable for
#' three full sweeps and the gap criterion is not yet met, a Newton polish on
#' the active set's stationarity conditions accelerates the last digits. The
#' sequential strong rule (`screen = TRUE`) discards groups whose score at the
#' previous `lambda` is far below the threshold and re-admits any that violate
#' their KKT condition, checked whenever the screened subproblem is
#' stationary (not only once the full certificate below is met, which a
#' wrongly screened group would keep open), so the solution is the same as
#' without screening.
#'
#' Convergence is certified by the duality gap `F(b) - D(nu)` at the dual point
#' `nu = r min{1, min_g lambda w_g / ||Q_g' r||}`; the solver stops when
#' `gap <= tol F(b)` and the maximum KKT violation is at most `tol lambda`
#' (both relative; a zero response column converges at once with `gap = 0`)
#' and reports the gap, the maximum KKT violation and the number of full
#' sweeps. The second criterion pins the coefficients, which a gap of `tol F`
#' alone leaves only `sqrt(tol)`-accurate. Reaching `max_iter` sweeps sets `converged = FALSE`, keeps the
#' last iterate, and raises one warning naming the proteins and `lambda`.
#'
#' Read-out: `beta_g = R_g^{-1} b_g`; `theta_c = s + d_c`; `gamma_c =
#' qr.coef(qr_Z_c, y_c - Psi_c theta_c)` from the design's profile; an edge
#' `k -> j` in cell `c` when the fitted-contribution norm `||Psi~_ck theta_ck||`
#' exceeds `edge_tol ||y~_c||` (`edge_tol = 0` reduces to any nonzero
#' coefficient). The discount `eta_c[k]` is reported
#' per regulator from the design's group norms `a_g` (Frobenius norms of the
#' group's columns, the column norm at `q = 1`), clamped at 1: `min(1, 1 -
#' min(W_c, W_S - W_c') / W_c)` with `W_g = mu_g a_g`.
#'
#' @param design the `public` list of `make_prnode_design()`: `x`, `y`,
#' `index`, `M_groups`, `D_groups`, `n_M`, `n_D`, `p`, `q` and (for `gamma`)
#' `profile`. `x` must be a base numeric matrix (a `Matrix` object is
#' refused) and `x` and `y` finite; a non-finite entry stops naming the
#' protein (in `y`) or the regulator and block (in `x`).
#' @param lambda the penalty path, positive and decreasing (one value allowed);
#' an increasing path is sorted with a warning.
#' @param edge_tol relative threshold of the edge rule.
#' @param max_iter maximum number of full sweeps per `lambda` and protein.
#' @param mu multipliers `c(S, M, D)`, positive.
#' @param screen use the sequential strong rule.
#' @param tol relative duality-gap tolerance, positive.
#' @param verbose a non-negative number; `0` prints nothing. `1` names each
#' protein as its path starts and, when it ends, how many `lambda` values
#' converged and the most sweeps any of them took. `2` adds one line per
#' `lambda`: the number of sweeps, the relative duality gap `gap / objective`,
#' the relative KKT violation `kkt_max / lambda`, the number of nonzero groups
#' out of `3p`, and whether the fit converged or hit `max_iter`. `3` adds one
#' line per sweep (every sweep up to 10, every tenth to 100, then every
#' hundredth) with the same certificate, plus a line whenever the Newton
#' polish is tried, which is where to look when a fit stalls.
#' @returns an object of class `prnode_path`: `s`, `d_M`, `d_D`, `theta_M`,
#' `theta_D` (`p x pq x L` arrays, row = equation, column `(k - 1) q + m`; the
#' rows are the columns of `design$y`, which may be a subset of the `p`
#' regulators),
#' `gamma_M`, `gamma_D` (`2 x p x L`: initial state and slope, in the design's
#' units), `networks` (list `S`, `M`, `D` of `p x p x L` logical arrays),
#' `objective`, `gap`, `kkt_max`, `iterations`, `converged` (`p x L`),
#' `discount` (list `M`, `D`, length `p`), `design` (the design without `x` and
#' `y`), `edge_tol`, `lambda`, `mu`, `screen` and `tol`.
#' @examples
#' ds <- build_md_dataset(seed = 1, num_d = 18, noise_level_ratio = 0.1,
#'                        truth_args = list(p = 4),
#'                        obs_time = seq(0, 1, length.out = 36),
#'                        steps_per_frame = 100)
#' lambda_max <- lambda_max_prnode(ds$public, bool_pooled = TRUE)
#' lambda_path <- lambda_max*10^seq(0, -2, length.out = 5)
#' path <- fit_prnode_path(ds$public, lambda = lambda_path)
#' names(path)
#' path$converged
#' path$networks$D[, , 5]     # the daughter's network at the smallest lambda
#' round(path$theta_D[, , 5], 2)
#' @export
fit_prnode_path <- function(design,
                       lambda,
                       edge_tol = 1e-8,
                       max_iter = 1e4,
                       mu = c(S = 1, M = 1, D = 1),
                       screen = TRUE,
                       tol = 1e-10,
                       verbose = 0){
  mu <- .mu_vec(mu)
  if(length(tol) != 1 || !is.finite(tol) || tol <= 0) stop("`tol` must be a positive number, not ", tol)
  if(length(verbose) != 1 || !is.numeric(verbose) || is.na(verbose) || verbose < 0){
    stop("`verbose` must be a single non-negative number (0 = silent, 1 = per protein, 2 = per lambda, 3 = per sweep)")
  }
  if(length(max_iter) != 1 || max_iter < 1) stop("`max_iter` must be at least 1")
  if(length(edge_tol) != 1 || edge_tol < 0) stop("`edge_tol` must be non-negative")
  if(!is.numeric(lambda) || length(lambda) < 1 || any(!is.finite(lambda)) || any(lambda <= 0)){
    stop("`lambda` must be a vector of positive finite values (an unpenalized fit is least squares, not a lasso path)")
  }
  if(is.unsorted(rev(lambda))){
    warning("`lambda` should decrease along the path; sorting it in decreasing order")
    lambda <- sort(lambda, decreasing = TRUE)
  }
  blocks <- .design_blocks(design, mu)
  groups <- blocks$groups
  p <- blocks$p
  q <- blocks$q
  pq <- p*q
  num_resp <- ncol(as.matrix(design$y))
  L <- length(lambda)
  N <- nrow(design$x)
  Q_full <- do.call(cbind, lapply(groups, function(gr) gr$Q))
  g_cols <- blocks$q_cols
  w_vec <- vapply(groups, function(gr) gr$w, numeric(1))

  arr <- function() array(0, c(num_resp, pq, L), dimnames = list(NULL, NULL, format(lambda, digits = 15)))
  out <- list(s = arr(), d_M = arr(), d_D = arr(), theta_M = arr(), theta_D = arr())
  gamma_M <- array(NA_real_, c(2, num_resp, L), dimnames = list(c("init", "slope"), NULL, format(lambda, digits = 15)))
  gamma_D <- gamma_M
  net <- function() array(FALSE, c(num_resp, p, L), dimnames = list(NULL, NULL, format(lambda, digits = 15)))
  networks <- list(S = net(), M = net(), D = net())
  # unnamed, so that `objective[j, l]` is a bare scalar (the arrays carry the
  # lambda names on their third dimension instead)
  scalar_mat <- function(fill) matrix(fill, num_resp, L)
  objective <- scalar_mat(NA_real_)
  gap <- scalar_mat(NA_real_)
  kkt_max <- scalar_mat(NA_real_)
  iterations <- scalar_mat(NA_integer_)
  converged <- scalar_mat(NA)

  m_rows <- seq_len(design$n_M)
  d_rows <- design$n_M + seq_len(design$n_D)
  s_cols <- seq_len(pq)
  m_cols <- pq + seq_len(pq)
  d_cols <- 2*pq + seq_len(pq)
  # which column of each cell's profile is response j: the same position when
  # `y` holds every protein, otherwise matched by value (a design whose `y`
  # was subset to some proteins keeps the full profile; test plan F2)
  col_M <- .profile_columns(design, "M", m_rows)
  col_D <- .profile_columns(design, "D", d_rows)
  not_converged <- character(0)
  for(j in seq_len(num_resp)){
    y_j <- as.matrix(design$y)[, j]
    if(verbose > 0){
      print(paste0("fit_prnode_path: protein ", j, " of ", num_resp, " (", L, " lambda value", if(L != 1) "s", ")"))
    }
    path <- .group_lasso_path(Q_full, g_cols, w_vec, y_j, lambda, tol = tol, max_iter = max_iter,
                              screen = screen, verbose = verbose - 1)
    if(verbose > 0){
      print(paste0("fit_prnode_path: protein ", j, " done: ", sum(path$converged), " of ", L, " converged, most sweeps ",
                   "at one lambda = ", max(path$iterations)))
    }
    for(l in seq_len(L)){
      b_vec <- path$b[, l]
      beta <- numeric(ncol(design$x))
      for(g in seq_along(groups)){
        beta[groups[[g]]$cols] <- backsolve(groups[[g]]$R, b_vec[g_cols[[g]]])
      }
      out$s[j, , l] <- beta[s_cols]
      out$d_M[j, , l] <- beta[m_cols]
      out$d_D[j, , l] <- beta[d_cols]
      out$theta_M[j, , l] <- beta[s_cols] + beta[m_cols]
      out$theta_D[j, , l] <- beta[s_cols] + beta[d_cols]
      objective[j, l] <- path$objective[l]
      gap[j, l] <- path$gap[l]
      kkt_max[j, l] <- path$kkt_max[l]
      iterations[j, l] <- path$iterations[l]
      converged[j, l] <- path$converged[l]
      if(!path$converged[l]) not_converged <- c(not_converged, paste0("protein ", j, " at lambda = ",
                                                                       signif(lambda[l], 6)))
      if(!is.null(design$profile)){
        gamma_M[, j, l] <- .nuisance_coef(design$profile$M, col_M[j], out$theta_M[j, , l])
        gamma_D[, j, l] <- .nuisance_coef(design$profile$D, col_D[j], out$theta_D[j, , l])
      }
      for(k in seq_len(p)){
        k_cols <- (k - 1)*q + seq_len(q)
        networks$M[j, k, l] <- .edge_present(design$x[m_rows, m_cols[k_cols], drop = FALSE],
                                             out$theta_M[j, k_cols, l], y_j[m_rows], edge_tol)
        networks$D[j, k, l] <- .edge_present(design$x[d_rows, d_cols[k_cols], drop = FALSE],
                                             out$theta_D[j, k_cols, l], y_j[d_rows], edge_tol)
        networks$S[j, k, l] <- .edge_present(design$x[, s_cols[k_cols], drop = FALSE],
                                             out$s[j, k_cols, l], y_j, edge_tol)
      }
    }
  }
  if(length(not_converged) > 0){
    warning("fit_prnode_path(): max_iter = ", max_iter, " full sweeps reached before the duality gap fell ",
            "below tol for ", paste(not_converged, collapse = "; "),
            "; the last iterate is returned with converged = FALSE")
  }

  design_settings <- design[setdiff(names(design), c("x", "y"))]
  structure(list(converged = converged,
                 d_D = out$d_D,
                 d_M = out$d_M,
                 design = design_settings,
                 discount = .discount(design, mu),
                 edge_tol = edge_tol,
                 gamma_D = gamma_D,
                 gamma_M = gamma_M,
                 gap = gap,
                 iterations = iterations,
                 kkt_max = kkt_max,
                 lambda = lambda,
                 mu = mu,
                 networks = networks,
                 objective = objective,
                 s = out$s,
                 screen = screen,
                 theta_D = out$theta_D,
                 theta_M = out$theta_M,
                 tol = tol),
            class = "prnode_path")
}

#' Normalize the multiplier argument
#' @noRd
.mu_vec <- function(mu){
  if(length(mu) != 3 || !is.numeric(mu) || any(!is.finite(mu)) || any(mu <= 0)){
    stop("`mu` must be three positive multipliers c(S, M, D); a zero multiplier leaves a group unpenalized, ",
         "for which the duality-gap certificate is undefined")
  }
  mu <- as.numeric(mu)
  names(mu) <- c("S", "M", "D")
  mu
}

#' Check the entries of a design, naming the offending column
#'
#' The one entry validator of every route into the solver (`.design_blocks()`,
#' hence `fit_prnode_path()`, `lambda_max_prnode()` and a numeric-`lambda`
#' `fit_prnode()`) and of `tune_prnode()`, which calls it before any path is
#' fitted. `x` must be a base numeric matrix (a sparse or dense `Matrix`
#' object is refused), `y` numeric, and both finite. A column of `y` is named
#' by its column name, else by `protein_names` when that has one name per
#' response column, else by its index; a column of `x` is mapped through the
#' column order of `.path_coef_mat()` (the `S`, `M`, `D` blocks of `pq`
#' columns, `q` per regulator) to its regulator and block, since
#' `make_prnode_design()` keeps no column names on `x`, the regulator named by
#' `protein_names` when that has one name per regulator, else by its index.
#'
#' @param q the group size.
#' @param arg the caller's name for the design (`"public"` or `"design"`),
#' used in the messages.
#' @returns `NULL`, invisibly; called for its errors.
#' @noRd
.check_design_entries <- function(x_mat,
                                  y_mat,
                                  q,
                                  protein_names,
                                  arg){
  if(inherits(x_mat, "Matrix")){
    stop("`", arg, "$x` must be a base matrix; sparse (Matrix) formats are not supported")
  }
  if(!is.matrix(x_mat) || !is.numeric(x_mat)) stop("`", arg, "$x` must be a numeric matrix")
  if(!is.numeric(y_mat)) stop("`", arg, "$y` must be numeric")
  p_y <- ncol(y_mat)
  y_names <- colnames(y_mat)
  if(is.null(y_names)) y_names <- if(length(protein_names) == p_y) protein_names else as.character(seq_len(p_y))
  bad_y <- which(colSums(!is.finite(y_mat)) > 0)
  if(length(bad_y) > 0){
    stop("`", arg, "$y` has non-finite values (NA, NaN or Inf) for protein(s) ",
         paste(y_names[bad_y], collapse = ", "))
  }
  bad_x <- which(colSums(!is.finite(x_mat)) > 0)
  if(length(bad_x) > 0){
    pq <- ncol(x_mat)/3
    p_reg <- pq/q
    reg_names <- if(length(protein_names) == p_reg) protein_names else as.character(seq_len(p_reg))
    block <- c("S", "M", "D")[ceiling(bad_x/pq)]
    regulator <- ceiling(((bad_x - 1) %% pq + 1)/q)
    label <- unique(paste0(reg_names[regulator], " (", block, " block)"))
    stop("`", arg, "$x` has non-finite values (NA, NaN or Inf) for regulator(s) ", paste(label, collapse = ", "))
  }
  invisible(NULL)
}

#' Validate a design and orthonormalize its groups
#'
#' @returns list `groups` (per group, in `index` label order of first
#' appearance: `label`, `cols`, `Q`, `R`, `type`, `w`), `q_cols` (each group's
#' columns in the stacked orthonormal matrix), `p`, `q`.
#' @noRd
.design_blocks <- function(design,
                           mu){
  if(!is.list(design)) stop("`design` must be the list returned by make_prnode_design()")
  for(nm in c("x", "y", "index", "M_groups", "D_groups", "n_M", "n_D")){
    if(is.null(design[[nm]])) stop("`design` is missing the element `", nm, "`")
  }
  x <- design$x
  y <- as.matrix(design$y)
  if(nrow(y) != nrow(x)) stop("`design$y` has ", nrow(y), " rows but `design$x` has ", nrow(x))
  if(length(design$index) != ncol(x)){
    stop("`index` has length ", length(design$index), " but `x` has ", ncol(x),
         " columns; every column needs a group label")
  }
  if(anyNA(design$index)) stop("`index` has NA labels; the two-block design has no unpenalized columns")
  if(nrow(x) != design$n_M + design$n_D) stop("`x` must have n_M + n_D rows")
  labels <- unique(design$index)
  if(length(labels) %% 3 != 0){
    stop("`index` has ", length(labels), " groups; 3p groups (scaffold, mother, daughter per regulator) are expected")
  }
  # p regulators from the groups; `y` may hold a subset of the proteins as
  # responses (test plan F2), so its column count is not p
  p <- length(labels)/3
  sizes <- vapply(labels, function(g) sum(design$index == g), integer(1))
  if(length(unique(sizes)) != 1) stop("every group of `index` must have the same size q")
  q <- unname(sizes[1])
  if(!is.null(design$q) && design$q != q) stop("`design$q` disagrees with the group sizes of `index`")
  .check_design_entries(x_mat = x,
                        y_mat = y,
                        q = q,
                        protein_names = design$protein_names,
                        arg = "design")
  groups <- lapply(seq_along(labels), function(g){
    cols <- which(design$index == labels[g])
    decomp <- qr(x[, cols, drop = FALSE])
    R <- qr.R(decomp)
    diag_abs <- abs(diag(R))
    if(decomp$rank < length(cols) || min(diag_abs) <= 1e-10*max(diag_abs)){
      stop("the group with label ", labels[g], " (columns ", paste(cols, collapse = ", "),
           ") is numerically singular (condition number of R about ",
           format(max(diag_abs)/max(min(diag_abs), .Machine$double.eps), digits = 3), ")")
    }
    type <- if(labels[g] %in% design$M_groups) "M" else if(labels[g] %in% design$D_groups) "D" else "S"
    list(label = labels[g], cols = cols, Q = qr.Q(decomp), R = R, type = type, w = sqrt(q)*mu[[type]])
  })
  starts <- cumsum(c(0, rep(q, length(groups) - 1)))
  q_cols <- lapply(seq_along(groups), function(g) starts[g] + seq_len(q))
  list(groups = groups, q_cols = q_cols, p = p, q = q)
}

#' Cyclic group soft-thresholding along a lambda path, orthonormal groups
#'
#' @param Q matrix with orthonormal column groups.
#' @param g_cols list of each group's columns of `Q`.
#' @param w_vec per-group penalty weights.
#' @param y the response.
#' @param lambda_vec decreasing path.
#' @param verbose `fit_prnode_path()`'s level minus one: above 0 one line per
#' lambda, above 1 one line per (throttled) sweep and per polish attempt.
#' @returns list `b` (`ncol(Q) x L`), `objective`, `gap`, `kkt_max`,
#' `iterations`, `converged` (length `L`).
#' @noRd
.group_lasso_path <- function(Q,
                              g_cols,
                              w_vec,
                              y,
                              lambda_vec,
                              tol,
                              max_iter,
                              screen,
                              inner_max = 10,
                              verbose = 0){
  G <- length(g_cols)
  L <- length(lambda_vec)
  b <- numeric(ncol(Q))
  r <- y
  b_mat <- matrix(0, ncol(Q), L)
  objective <- gap <- kkt_max <- numeric(L)
  iterations <- integer(L)
  converged <- logical(L)
  group_norm <- function(vec) vapply(g_cols, function(cols) sqrt(sum(vec[cols]^2)), numeric(1))
  score_prev <- NULL

  for(l in seq_len(L)){
    lambda <- lambda_vec[l]
    candidates <- seq_len(G)
    if(screen && !is.null(score_prev) && l > 1){
      # sequential strong rule on ||z_g + b_g|| at the previous lambda
      keep <- score_prev >= w_vec*(2*lambda - lambda_vec[l - 1])
      candidates <- which(keep | group_norm(b) > 0)
      if(verbose > 1){
        print(paste0("  lambda ", l, " of ", L, ": the strong rule keeps ", length(candidates), " of ", G, " groups"))
      }
    }
    iter <- 0L
    stable_count <- 0L
    active_prev <- NULL
    bool_done <- FALSE
    repeat{
      b_before <- b
      sweep <- .cd_sweep(Q, g_cols, w_vec, lambda, candidates, b, r)
      b <- sweep$b
      r <- sweep$r
      iter <- iter + 1L
      active <- which(group_norm(b) > 0)
      for(it in seq_len(inner_max)){
        if(length(active) == 0) break
        sweep <- .cd_sweep(Q, g_cols, w_vec, lambda, active, b, r)
        b <- sweep$b
        r <- sweep$r
        if(sweep$max_change <= 1e-8*(1 + sqrt(sum(b^2)))) break
      }
      # extrapolate along the displacement of this outer iteration while the
      # objective keeps falling: cyclic descent crawls along flat valleys of
      # ill-conditioned designs, and this safeguarded over-relaxation cuts
      # the number of sweeps by orders of magnitude there
      note <- ""
      if(iter >= 2){
        acc <- .accelerate(Q, g_cols, w_vec, lambda, b, r, y, b - b_before)
        if(verbose > 1 && any(acc$b != b)) note <- ", extrapolated"
        b <- acc$b
        r <- acc$r
      }
      cert <- .gap_certificate(Q, g_cols, w_vec, lambda, b, r, y)
      if(verbose > 1 && .verbose_tick(iter)){
        print(paste0("    sweep ", iter, ": active groups = ", sum(group_norm(b) > 0), " of ", G,
                     ", gap/objective = ", signif(.ratio(cert$gap, cert$objective), 3),
                     ", kkt/lambda = ", signif(cert$kkt_max/lambda, 3), note))
      }
      # both criteria scale with tol: the gap certifies the objective, the
      # KKT residual the coefficients (a gap of tol F alone leaves the
      # coefficients only sqrt(tol)-accurate, which flips near-boundary
      # supports between sweep orders; test plan I1, R6, R7)
      # re-admission: a screened-out group violating its KKT condition (its
      # score is computed with the certificate, so the set is free every sweep)
      viol <- setdiff(which(cert$score > w_vec*lambda*(1 + 1e-10)), candidates)
      if(cert$gap <= tol*cert$objective && cert$kkt_max <= tol*lambda){
        if(length(viol) == 0){
          bool_done <- TRUE
        } else {
          candidates <- sort(union(candidates, viol))
        }
      } else {
        # the certificate is over all groups, so a wrongly screened group keeps
        # the gap open while the screened subproblem is already stationary; the
        # re-admission must not wait for it (estimator plan R7b). A candidate
        # set is stationary when its KKT residual is within tol (or it is
        # empty); the periodic fallback covers a subproblem that crawls.
        if(length(viol) > 0 &&
           (length(candidates) == 0 || max(cert$kkt[candidates]) <= tol*lambda || iter %% 25 == 0)){
          if(verbose > 1){
            print(paste0("    sweep ", iter, ": re-admitting ", length(viol), " screened-out group",
                         if(length(viol) != 1) "s", " violating KKT"))
          }
          candidates <- sort(union(candidates, viol))
        }
        active <- which(group_norm(b) > 0)
        stable_count <- if(identical(active, active_prev)) stable_count + 1L else 0L
        active_prev <- active
        # polish once the active set has settled, and periodically regardless
        # (a flickering near-zero group would otherwise postpone it forever)
        if((stable_count >= 3 || iter %% 25 == 0) && length(active) > 0){
          polished <- .newton_polish(Q, g_cols, w_vec, lambda, active, b, y)
          if(verbose > 1 && .verbose_tick(iter)){
            print(paste0("    sweep ", iter, ": Newton polish on ", length(active), " active groups ",
                         if(is.null(polished)) "declined (singular Hessian or no decrease)" else
                           paste0("applied, active groups now ", sum(group_norm(polished$b) > 0))))
          }
          if(!is.null(polished)){
            b <- polished$b
            r <- polished$r
            stable_count <- 0L
          }
        }
      }
      if(bool_done) break
      if(iter >= max_iter){
        break
      }
    }
    cert <- .gap_certificate(Q, g_cols, w_vec, lambda, b, r, y)
    if(verbose > 0){
      print(paste0("  lambda ", l, " of ", L, " = ", signif(lambda, 4), ": ", iter, " sweep", if(iter != 1) "s",
                   ", gap/objective = ", signif(.ratio(cert$gap, cert$objective), 3),
                   ", kkt/lambda = ", signif(cert$kkt_max/lambda, 3),
                   ", active groups = ", sum(group_norm(b) > 0), " of ", G, ", ",
                   if(bool_done) "converged" else paste0("not converged (max_iter = ", max_iter, " sweeps reached)")))
    }
    b_mat[, l] <- b
    objective[l] <- cert$objective
    gap[l] <- cert$gap
    kkt_max[l] <- cert$kkt_max
    iterations[l] <- iter
    converged[l] <- bool_done
    score_prev <- cert$score_partial
  }
  list(b = b_mat, objective = objective, gap = gap, kkt_max = kkt_max, iterations = iterations,
       converged = converged)
}

#' Which sweeps the level-3 trace prints: every sweep up to 10, every tenth
#' to 100, then every hundredth
#' @noRd
.verbose_tick <- function(iter){
  iter <= 10 || (iter <= 100 && iter %% 10 == 0) || iter %% 100 == 0
}

#' A ratio that is 0 when its denominator is 0 (a zero response column)
#' @noRd
.ratio <- function(numerator,
                   denominator){
  if(denominator > 0) numerator/denominator else 0
}

#' Safeguarded extrapolation along a displacement
#'
#' Doubles the step along `direction` from `b` while the objective decreases
#' and returns the best point found (the input when no step helps).
#' @noRd
.accelerate <- function(Q,
                        g_cols,
                        w_vec,
                        lambda,
                        b,
                        r,
                        y,
                        direction,
                        max_doublings = 20){
  if(all(direction == 0)) return(list(b = b, r = r))
  objective_of <- function(bb, rr){
    pen <- sum(vapply(seq_along(g_cols), function(g) w_vec[g]*sqrt(sum(bb[g_cols[[g]]]^2)), numeric(1)))
    0.5*sum(rr^2) + lambda*pen
  }
  Q_dir <- c(Q %*% direction)
  f_best <- objective_of(b, r)
  b_best <- b
  r_best <- r
  alpha <- 1
  for(k in seq_len(max_doublings)){
    b_try <- b + alpha*direction
    r_try <- r - alpha*Q_dir
    f_try <- objective_of(b_try, r_try)
    if(!(f_try < f_best - 1e-15*abs(f_best))) break
    f_best <- f_try
    b_best <- b_try
    r_best <- r_try
    alpha <- 2*alpha
  }
  list(b = b_best, r = r_best)
}

#' One sweep of exact block updates over the given groups
#' @noRd
.cd_sweep <- function(Q,
                      g_cols,
                      w_vec,
                      lambda,
                      which_groups,
                      b,
                      r){
  max_change <- 0
  for(g in which_groups){
    cols <- g_cols[[g]]
    Q_g <- Q[, cols, drop = FALSE]
    b_old <- b[cols]
    u <- c(crossprod(Q_g, r)) + b_old
    nu <- sqrt(sum(u^2))
    b_new <- if(nu > lambda*w_vec[g]) (1 - lambda*w_vec[g]/nu)*u else rep(0, length(cols))
    delta <- b_new - b_old
    if(any(delta != 0)){
      r <- r - c(Q_g %*% delta)
      b[cols] <- b_new
      max_change <- max(max_change, sqrt(sum(delta^2)))
    }
  }
  list(b = b, r = r, max_change = max_change)
}

#' Objective, duality gap and KKT violation at an iterate
#' @noRd
.gap_certificate <- function(Q,
                             g_cols,
                             w_vec,
                             lambda,
                             b,
                             r,
                             y){
  z <- c(crossprod(Q, r))
  G <- length(g_cols)
  score <- numeric(G)
  score_partial <- numeric(G)
  kkt <- numeric(G)
  penalty <- 0
  for(g in seq_len(G)){
    cols <- g_cols[[g]]
    z_g <- z[cols]
    b_g <- b[cols]
    nb <- sqrt(sum(b_g^2))
    score[g] <- sqrt(sum(z_g^2))
    score_partial[g] <- sqrt(sum((z_g + b_g)^2))
    penalty <- penalty + w_vec[g]*nb
    kkt[g] <- if(nb > 0) sqrt(sum((z_g - lambda*w_vec[g]*b_g/nb)^2)) else max(0, score[g] - lambda*w_vec[g])
  }
  objective <- 0.5*sum(r^2) + lambda*penalty
  ratio_vec <- ifelse(score > 0, lambda*w_vec/score, Inf)
  scale_val <- min(1, min(ratio_vec))
  nu <- scale_val*r
  dual <- 0.5*sum(y^2) - 0.5*sum((y - nu)^2)
  list(objective = objective,
       gap = max(objective - dual, 0),
       kkt = kkt,
       kkt_max = max(kkt),
       score = score,
       score_partial = score_partial)
}

#' Newton polish of the stationarity conditions on the active set
#'
#' Solves `Q_A'(y - Q_A b_A) = lambda W_A b_g / ||b_g||` by Newton's method
#' from the current iterate. A group whose Newton step flips its direction
#' (for single-column groups, its sign) or shrinks it below `1e-10` of the
#' iterate is dropped from the active set (set exactly to zero) and the
#' remaining groups are re-solved, which is the active-set method of the lasso
#' when `q = 1`; a dropped group that violates its KKT condition is picked up
#' again by the next coordinate-descent sweep. Returns `NULL` when the Hessian
#' is singular or the objective does not decrease, in which case the
#' coordinate descent continues unchanged.
#' @noRd
.newton_polish <- function(Q,
                           g_cols,
                           w_vec,
                           lambda,
                           active,
                           b,
                           y,
                           max_steps = 50){
  objective_of <- function(bb){
    pen <- sum(vapply(seq_along(g_cols), function(g) w_vec[g]*sqrt(sum(bb[g_cols[[g]]]^2)), numeric(1)))
    0.5*sum((y - Q %*% bb)^2) + lambda*pen
  }
  f_old <- objective_of(b)
  b_cur <- b
  for(step in seq_len(max_steps)){
    if(length(active) == 0) break
    a_cols <- unlist(g_cols[active])
    Q_A <- Q[, a_cols, drop = FALSE]
    b_A <- b_cur[a_cols]
    local_cols <- lapply(active, function(g) match(g_cols[[g]], a_cols))
    r_A <- y - c(Q_A %*% b_A)
    grad <- -c(crossprod(Q_A, r_A))
    H <- crossprod(Q_A)
    for(i in seq_along(active)){
      cols <- local_cols[[i]]
      b_g <- b_A[cols]
      nb <- sqrt(sum(b_g^2))
      if(nb == 0) return(NULL)
      grad[cols] <- grad[cols] + lambda*w_vec[active[i]]*b_g/nb
      H[cols, cols] <- H[cols, cols] + lambda*w_vec[active[i]]/nb*(diag(length(cols)) - tcrossprod(b_g)/nb^2)
    }
    delta <- tryCatch(solve(H, -grad), error = function(e) NULL)
    if(is.null(delta) || any(!is.finite(delta))) return(NULL)
    b_try <- b_A + delta
    # groups whose direction flips or collapses leave the active set
    drop <- vapply(seq_along(active), function(i){
      cols <- local_cols[[i]]
      sum(b_try[cols]*b_A[cols]) <= 0 || sqrt(sum(b_try[cols]^2)) < 1e-10*(1 + sqrt(sum(b_A^2)))
    }, logical(1))
    b_new <- numeric(length(b))
    if(any(drop)){
      keep_cols <- unlist(local_cols[!drop])
      b_try[unlist(local_cols[drop])] <- 0
      b_new[a_cols] <- b_try
      active <- active[!drop]
      b_cur <- b_new
      next
    }
    b_new[a_cols] <- b_try
    b_cur <- b_new
    if(sqrt(sum(delta^2)) <= 1e-15*(1 + sqrt(sum(b_try^2)))) break
  }
  f_new <- objective_of(b_cur)
  if(!is.finite(f_new) || f_new > f_old*(1 + 1e-12) + 1e-300) return(NULL)
  list(b = b_cur, r = y - c(Q %*% b_cur))
}

#' Nuisance coefficients of one protein and cell from the profile
#' @noRd
.nuisance_coef <- function(prof,
                           j,
                           theta_vec){
  if(is.na(j)) return(c(NA_real_, NA_real_))
  c(qr.coef(prof$qr_Z, prof$y[, j] - c(prof$Psi %*% theta_vec)))
}

#' Map the columns of `design$y` to the columns of a cell's profile
#'
#' Position when the counts agree; otherwise each response column is matched
#' to the profile's weighted profiled columns by value (NA when no unique match).
#' @noRd
.profile_columns <- function(design,
                             cell,
                             rows){
  prof <- design$profile[[cell]]
  y <- as.matrix(design$y)
  if(is.null(prof) || is.null(prof$y_tilde)) return(rep(NA_integer_, ncol(y)))
  if(ncol(prof$y_tilde) == ncol(y)) return(seq_len(ncol(y)))
  omega_c <- if(is.null(design$omega)) 1 else unname(design$omega[[if(cell == "M") 1 else 2]])
  candidates <- sqrt(omega_c)*as.matrix(prof$y_tilde)
  vapply(seq_len(ncol(y)), function(j){
    hits <- which(vapply(seq_len(ncol(candidates)), function(k){
      isTRUE(all.equal(unname(candidates[, k]), unname(y[rows, j]), tolerance = 1e-10))
    }, logical(1)))
    if(length(hits) == 1) hits else NA_integer_
  }, integer(1))
}

#' The fitted-contribution edge rule
#' @noRd
.edge_present <- function(x_block,
                          theta_vec,
                          y_vec,
                          edge_tol){
  contribution <- c(x_block %*% theta_vec)
  sqrt(sum(contribution^2)) > edge_tol*sqrt(sum(y_vec^2))
}

#' The discounts, per regulator and cell, clamped at 1
#' @noRd
.discount <- function(design,
                      mu){
  # p regulators from the groups (not from `y`, which may hold a subset of
  # the proteins), q from the group size
  p <- length(unique(design$index))/3
  q <- ncol(design$x)/(3*p)
  pq <- p*q
  norm_of <- function(cols) sqrt(sum(design$x[, cols, drop = FALSE]^2))
  eta_M <- eta_D <- numeric(p)
  for(k in seq_len(p)){
    k_cols <- (k - 1)*q + seq_len(q)
    W_S <- mu[["S"]]*norm_of(k_cols)
    W_M <- mu[["M"]]*norm_of(pq + k_cols)
    W_D <- mu[["D"]]*norm_of(2*pq + k_cols)
    eta_D[k] <- if(W_D > 0) min(1, 1 - min(W_D, W_S - W_M)/W_D) else NA_real_
    eta_M[k] <- if(W_M > 0) min(1, 1 - min(W_M, W_S - W_D)/W_M) else NA_real_
  }
  list(M = eta_M, D = eta_D)
}

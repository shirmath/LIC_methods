#' Rewire a regulatory network at the level of its edge units
#'
#' A daughter editor for simulating pairs at a chosen level of inheritance. The off-diagonal positions of `ode_pars` are grouped into **units**:
#' an antisymmetric pair `(j, k), (k, j)` when both orientations carry a
#' coefficient (the damped rotations `sim_md_network()` draws), and a single
#' edge otherwise. Each unit is kept with probability `keep_prob` and dropped
#' otherwise; a dropped unit is, with probability `replace_prob`, replaced by
#' a fresh unit of the same kind at an unordered position that none of the
#' kept units occupies, so that at `replace_prob = 1` **the density is
#' preserved** and a dropped position may be redrawn by chance (at
#' `keep_prob = 0` the result's topology is an independent draw at the input's
#' density, not a disjoint one), while at `replace_prob = 0` the result is a
#' thinned copy of the input (the daughter loses edges and gains none). A kept
#' unit is, with probability
#' `flip_prob`, **sign-flipped** (both orientations multiplied by minus a
#' log-uniform draw in `unequal_factor`); otherwise, with probability
#' `unequal_prob`, **reweighted** by a positive log-uniform draw; otherwise
#' copied exactly. Both orientations of a pair get one factor, so a pair stays
#' an exact rotation (a flipped pair rotates the other way) and a system of
#' antisymmetric pairs stays stable at `q = 1`. The diagonal is never changed.
#'
#' A replacement's coefficients are drawn as `sim_md_network()` draws them:
#' one magnitude uniform in `edge_magnitude` shared by the two orientations of
#' a pair (opposite signs, random orientation), the higher powers uniform in
#' `edge_magnitude` per position, each `m`-th power divided by `factorial(m)`.
#' Donnat and Holmes (2018) call this the rewiring model; the weighted Jaccard
#' distance is a one-to-one function of the rewiring rate whatever the
#' density.
#'
#' Only reproducibility is promised for the random stream: the same
#' `seed_number` gives the same output, and `seed_number = NULL` consumes the
#' caller's stream deterministically. The number of draws is an implementation
#' detail.
#'
#' @param ode_pars `p x pq` coefficient matrix, row = equation, column
#' `(k - 1) q + m` = the coefficient of regulator `k`'s `m`-th power (the
#' `m_ode_pars` of a `sim_md_network()` truth).
#' @param keep_prob probability that a unit is kept, in `[0, 1]`.
#' @param edge_magnitude range of the uniform draw of a replacement's
#' coefficient magnitudes.
#' @param flip_prob probability that a kept unit is sign-flipped, in `[0, 1]`.
#' @param q degree of the polynomial basis; `NULL` derives it from the
#' dimensions.
#' @param replace_prob probability that a dropped unit is replaced by a fresh
#' one, in `[0, 1]`; `1` preserves the density.
#' @param unequal_factor range of the log-uniform factor of a reweighted or
#' flipped unit; both ends positive. The default halves or doubles at most,
#' so a rotation's frequency stays within what 48 frames resolve.
#' @param unequal_prob probability that a kept, unflipped unit is reweighted,
#' in `[0, 1]`.
#' @param seed_number seed set before the first draw; `NULL` leaves the stream
#' to the caller.
#' @returns a list with `pars` (the rewired `p x pq` matrix) and `units_df`
#' (a data frame with one row per unit of the input and one per replacement:
#' `j`, `k` (the unordered position, `j < k`), `kind` (`"pair"` or
#' `"single"`), `orientation` (`1` for `(j, k)`, `-1` for `(k, j)`, `0` for a
#' pair), `status` (`"kept_equal"`, `"kept_unequal"`, `"kept_flipped"`,
#' `"dropped"` or `"added"`) and `factor` (the multiplier of a kept unit, `NA`
#' otherwise)); a dropped unit that was not replaced has no `"added"` row.
#' @examples
#' truth <- sim_md_network(p = 6, num_shared = 8, num_mother_only = 0,
#'                         num_daughter_only = 0, frac_unequal = 0,
#'                         seed_number = 1)
#' rewired <- rewire_network(truth$m_ode_pars, keep_prob = 0.5,
#'                           unequal_prob = 0.5, flip_prob = 0.2,
#'                           seed_number = 2)
#' rewired$units_df
#' # the density is preserved
#' c(sum(truth$m_ode_pars != 0), sum(rewired$pars != 0))
#' @export
rewire_network <- function(ode_pars,
                           keep_prob,
                           edge_magnitude = c(4*pi, 10*pi),
                           flip_prob = 0,
                           q = NULL,
                           replace_prob = 1,
                           unequal_factor = c(1/2, 2),
                           unequal_prob = 0,
                           seed_number = 10){
  stopifnot(is.matrix(ode_pars), is.numeric(ode_pars), all(is.finite(ode_pars)),
            length(edge_magnitude) == 2, length(unequal_factor) == 2, all(unequal_factor > 0),
            length(replace_prob) == 1, replace_prob >= 0, replace_prob <= 1)
  p <- nrow(ode_pars)
  if(is.null(q)) q <- ncol(ode_pars)/p
  if(p < 2 || q < 1 || q != round(q) || ncol(ode_pars) != p*q){
    stop("`ode_pars` must be p x pq with p >= 2 and an integer q, not ", p, " x ", ncol(ode_pars))
  }
  for(arg in c("keep_prob", "flip_prob", "unequal_prob")){
    val <- get(arg)
    if(length(val) != 1 || !is.finite(val) || val < 0 || val > 1){
      stop("`", arg, "` must be one number in [0, 1], not ", paste(val, collapse = ", "))
    }
  }
  if(!is.null(seed_number)) set.seed(seed_number)

  block_cols <- function(k) (k - 1)*q + seq_len(q)
  has_edge <- function(pars, j, k) any(pars[j, block_cols(k)] != 0)
  # the units of the input: an unordered position {j < k} is a pair when both
  # orientations are nonzero and a single when exactly one is
  pos_mat <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  units_list <- list()
  for(r in seq_len(nrow(pos_mat))){
    j <- pos_mat[r, 1]
    k <- pos_mat[r, 2]
    jk <- has_edge(ode_pars, j, k)
    kj <- has_edge(ode_pars, k, j)
    if(jk && kj){
      units_list[[length(units_list) + 1]] <- data.frame(j = j, k = k, kind = "pair", orientation = 0,
                                                         stringsAsFactors = FALSE)
    } else if(jk || kj){
      units_list[[length(units_list) + 1]] <- data.frame(j = j, k = k, kind = "single",
                                                         orientation = if(jk) 1 else -1,
                                                         stringsAsFactors = FALSE)
    }
  }
  units_df <- if(length(units_list) > 0) do.call(rbind, units_list) else
    data.frame(j = integer(0), k = integer(0), kind = character(0), orientation = numeric(0),
               stringsAsFactors = FALSE)
  units_df$status <- rep("kept_equal", nrow(units_df))
  units_df$factor <- rep(1, nrow(units_df))

  # the positions (row, column) a unit occupies
  unit_positions <- function(unit){
    if(unit$kind == "pair") return(rbind(c(unit$j, unit$k), c(unit$k, unit$j)))
    if(unit$orientation == 1) rbind(c(unit$j, unit$k)) else rbind(c(unit$k, unit$j))
  }
  log_factor_draw <- function() exp(stats::runif(1, log(unequal_factor[1]), log(unequal_factor[2])))

  # keep, reweight or drop each unit of the input
  new_pars <- ode_pars
  for(u in seq_len(nrow(units_df))){
    unit <- units_df[u, ]
    positions <- unit_positions(unit)
    if(stats::runif(1) < keep_prob){
      if(stats::runif(1) < flip_prob){
        factor_val <- -log_factor_draw()
        units_df$status[u] <- "kept_flipped"
      } else if(stats::runif(1) < unequal_prob){
        factor_val <- log_factor_draw()
        units_df$status[u] <- "kept_unequal"
      } else {
        factor_val <- 1
      }
      units_df$factor[u] <- factor_val
      for(r in seq_len(nrow(positions))){
        cols <- block_cols(positions[r, 2])
        new_pars[positions[r, 1], cols] <- factor_val*new_pars[positions[r, 1], cols]
      }
    } else {
      units_df$status[u] <- "dropped"
      units_df$factor[u] <- NA_real_
      for(r in seq_len(nrow(positions))) new_pars[positions[r, 1], block_cols(positions[r, 2])] <- 0
    }
  }

  # replace a dropped unit (with probability replace_prob) at a position none
  # of the kept units occupies; a dropped position is available again, which
  # is what makes keep_prob = 0 an independent draw rather than a disjoint one.
  # The replacement coin is tossed for every dropped unit before any
  # replacement is placed, so the number of draws does not depend on the
  # positions chosen.
  dropped_idx <- which(units_df$status == "dropped")
  bool_replace <- stats::runif(length(dropped_idx)) < replace_prob
  for(u in dropped_idx[bool_replace]){
    kind <- units_df$kind[u]
    available <- which(vapply(seq_len(nrow(pos_mat)), function(r){
      !has_edge(new_pars, pos_mat[r, 1], pos_mat[r, 2]) && !has_edge(new_pars, pos_mat[r, 2], pos_mat[r, 1])
    }, logical(1)))
    if(length(available) == 0){
      stop("no empty position is left to place a replacement unit; `ode_pars` is too dense to rewire")
    }
    # `sample(x)` on a length-1 vector permutes seq_len(x), so draw by index
    chosen <- available[sample.int(length(available), 1)]
    j <- pos_mat[chosen, 1]
    k <- pos_mat[chosen, 2]
    magnitude <- stats::runif(1, edge_magnitude[1], edge_magnitude[2])
    orient <- sample(c(-1, 1), 1)
    draw_block <- function(sign_val){
      higher <- if(q > 1) stats::runif(q - 1, edge_magnitude[1], edge_magnitude[2]) else numeric(0)
      sign_val*c(magnitude, higher)/factorial(seq_len(q))
    }
    if(kind == "pair"){
      new_pars[j, block_cols(k)] <- draw_block(orient)
      new_pars[k, block_cols(j)] <- draw_block(-orient)
      orientation <- 0
    } else {
      # a single keeps one orientation, chosen at random, with a random sign
      orientation <- sample(c(-1, 1), 1)
      if(orientation == 1) new_pars[j, block_cols(k)] <- draw_block(orient) else
        new_pars[k, block_cols(j)] <- draw_block(orient)
    }
    units_df <- rbind(units_df, data.frame(j = j, k = k, kind = kind, orientation = orientation,
                                           status = "added", factor = NA_real_, stringsAsFactors = FALSE))
  }
  rownames(units_df) <- NULL

  # the density is preserved exactly when every dropped unit is replaced, never
  # exceeded, and the diagonal is untouched
  off_in <- row(matrix(0, p, p)) != col(matrix(0, p, p))
  density_in <- sum(.network_of(ode_pars, q)[off_in])
  density_out <- sum(.network_of(new_pars, q)[off_in])
  stopifnot(density_out <= density_in,
            replace_prob < 1 || density_out == density_in,
            all(diag(.network_of(new_pars, q)) == diag(.network_of(ode_pars, q))))
  list(pars = new_pars,
       units_df = units_df)
}

#' Assemble a mother-daughter truth from two coefficient matrices
#'
#' The inverse of what `sim_md_network()` does on its way out: given the two
#' cells' systems, it recovers the scaffold-plus-deviation decomposition and
#' the edge types the package's scorer, plots and dataset builder read
#' (`build_md_dataset(truth = )`, `score_md_networks()`). Position by position
#' (the `q` columns of regulator `k` in equation `j`): a position nonzero in
#' one cell only is `mother_only` or `daughter_only` (that cell's deviation
#' holds it); a position with identical coefficients is `shared_equal` (the
#' scaffold holds it); a position present in both with the **same sign** of
#' the first-power coefficient is `shared_unequal`, the scaffold holds the
#' cell with the smaller first-power magnitude and the other cell's deviation
#' the difference (the generator's convention); a position present in both
#' with **opposite signs** is `shared_unequal` with a zero scaffold and each
#' cell's deviation equal to its full coefficients (nothing is common to a
#' flipped edge). The diagonal is `self`: the scaffold holds the mother's
#' self-decay and the daughter's deviation any difference. Because a flipped
#' edge is `shared_unequal` like any reweighted one, the returned
#' `flip_network` marks it (the scorer's confusion rows do not see it).
#'
#' @param m_ode_pars the mother's `p x pq` coefficient matrix.
#' @param d_ode_pars the daughter's, same dimensions.
#' @param d_intercept the daughter's intercepts (length `p`); `NULL` copies
#' the mother's.
#' @param m_intercept the mother's intercepts (length `p`); `NULL` is zeros.
#' @param q degree of the polynomial basis; `NULL` derives it from the
#' dimensions.
#' @returns a list in the shape of `sim_md_network()`'s: `d_delta`,
#' `d_intercept`, `d_network`, `d_ode_pars`, `edge_type`, `m_delta`,
#' `m_intercept`, `m_network`, `m_ode_pars`, `p`, `q`, `s_network`, `s_pars`,
#' with `m_ode_pars == s_pars + m_delta` and `d_ode_pars == s_pars + d_delta`
#' to floating-point rounding, plus `flip_network` (`p x p` logical, `TRUE`
#' at a position present in both cells with opposite first-power signs).
#' `s_network` is `m_network & d_network`, the package's convention (the
#' edges both cells have), so a flipped edge is in `s_network` although its
#' `s_pars` block is zero: the scaffold MCC of `score_md_networks()` is not
#' a meaningful score on a truth with flips.
#' @examples
#' truth <- sim_md_network(p = 6, num_shared = 8, num_mother_only = 0,
#'                         num_daughter_only = 0, frac_unequal = 0,
#'                         seed_number = 1)
#' rewired <- rewire_network(truth$m_ode_pars, keep_prob = 0.5,
#'                           unequal_prob = 0.5, flip_prob = 0.2,
#'                           seed_number = 2)
#' pair <- md_truth_from_pars(truth$m_ode_pars, rewired$pars)
#' table(pair$edge_type)
#' ds <- build_md_dataset(seed = 3, num_d = 24, noise_level_ratio = 0.5,
#'                        truth = pair, steps_per_frame = 100)
#' @export
md_truth_from_pars <- function(m_ode_pars,
                               d_ode_pars,
                               d_intercept = NULL,
                               m_intercept = NULL,
                               q = NULL){
  stopifnot(is.matrix(m_ode_pars), is.matrix(d_ode_pars), is.numeric(m_ode_pars), is.numeric(d_ode_pars),
            all(dim(m_ode_pars) == dim(d_ode_pars)), all(is.finite(m_ode_pars)), all(is.finite(d_ode_pars)))
  p <- nrow(m_ode_pars)
  if(is.null(q)) q <- ncol(m_ode_pars)/p
  if(p < 1 || q < 1 || q != round(q) || ncol(m_ode_pars) != p*q){
    stop("`m_ode_pars` must be p x pq with an integer q, not ", p, " x ", ncol(m_ode_pars))
  }
  if(is.null(m_intercept)) m_intercept <- rep(0, p)
  if(is.null(d_intercept)) d_intercept <- m_intercept
  stopifnot(length(m_intercept) == p, length(d_intercept) == p)

  block_cols <- function(k) (k - 1)*q + seq_len(q)
  s_pars <- matrix(0, p, p*q)
  m_delta <- matrix(0, p, p*q)
  d_delta <- matrix(0, p, p*q)
  edge_type <- matrix("absent", p, p)
  flip_network <- matrix(FALSE, p, p)
  for(j in seq_len(p)) for(k in seq_len(p)){
    cols <- block_cols(k)
    m_blk <- m_ode_pars[j, cols]
    d_blk <- d_ode_pars[j, cols]
    m_on <- any(m_blk != 0)
    d_on <- any(d_blk != 0)
    if(j == k){
      edge_type[j, k] <- "self"
      s_pars[j, cols] <- m_blk
      d_delta[j, cols] <- d_blk - m_blk
    } else if(m_on && !d_on){
      edge_type[j, k] <- "mother_only"
      m_delta[j, cols] <- m_blk
    } else if(d_on && !m_on){
      edge_type[j, k] <- "daughter_only"
      d_delta[j, cols] <- d_blk
    } else if(m_on && d_on){
      if(all(m_blk == d_blk)){
        edge_type[j, k] <- "shared_equal"
        s_pars[j, cols] <- m_blk
      } else {
        edge_type[j, k] <- "shared_unequal"
        # the first power carries the sign of the edge; equal signs share the
        # smaller cell, opposite signs share nothing
        bool_same_sign <- m_blk[1] != 0 && d_blk[1] != 0 && sign(m_blk[1]) == sign(d_blk[1])
        if(bool_same_sign && abs(m_blk[1]) <= abs(d_blk[1])){
          s_pars[j, cols] <- m_blk
          d_delta[j, cols] <- d_blk - m_blk
        } else if(bool_same_sign){
          s_pars[j, cols] <- d_blk
          m_delta[j, cols] <- m_blk - d_blk
        } else {
          m_delta[j, cols] <- m_blk
          d_delta[j, cols] <- d_blk
          flip_network[j, k] <- m_blk[1] != 0 && d_blk[1] != 0
        }
      }
    }
  }
  m_network <- .network_of(m_ode_pars, q)
  d_network <- .network_of(d_ode_pars, q)
  # s + (d - s) can differ from d in the last bit, so the identities hold to
  # rounding, not bit for bit
  stopifnot(isTRUE(all.equal(m_ode_pars, s_pars + m_delta, check.attributes = FALSE)),
            isTRUE(all.equal(d_ode_pars, s_pars + d_delta, check.attributes = FALSE)))
  list(d_delta = d_delta,
       d_intercept = as.numeric(d_intercept),
       d_network = d_network,
       d_ode_pars = d_ode_pars,
       edge_type = edge_type,
       flip_network = flip_network,
       m_delta = m_delta,
       m_intercept = as.numeric(m_intercept),
       m_network = m_network,
       m_ode_pars = m_ode_pars,
       p = p,
       q = q,
       s_network = m_network & d_network,
       s_pars = s_pars)
}

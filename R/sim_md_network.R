#' Draw a mother-daughter pair of regulatory networks
#'
#' The truth generator of the PRNODE simulations. The two cells share a scaffold and each has a
#' deviation, so that every edge type the estimator must tell apart is
#' present: shared with equal coefficients, shared with unequal coefficients,
#' mother-only, daughter-only and absent. Every diagonal entry is a negative
#' self-degradation rate, equal in both cells, so that no linear combination
#' of the state is conserved and the profiled design has full column rank
#' (a generator that conserved a triplet's sum would violate this).
#'
#' @details
#' **Edges come in antisymmetric pairs.** Off-diagonal edges are drawn as
#' unordered pairs of positions `(j, k)` and `(k, j)` with coefficients of
#' opposite sign, `+m` and `-m`, so that each pair is a damped rotation of
#' the two proteins at angular frequency `m`; a count that is
#' odd leaves one single position. With the default `edge_magnitude = c(4 pi,
#' 10 pi)` a pair completes two to five cycles per unit of time (the movie; a
#' yeast movie of nine hours holds about five cell cycles), which keeps every
#' observation window, the short daughter's included, on a curved trajectory. This is what makes the integrated designs
#' well-conditioned: with edges of random sign the dynamics are damped
#' exponentials, whose integrals are nearly collinear (condition numbers of
#' 1e7 to 1e13 on the 48-frame designs, and no network recovery at any noise
#' level), while sums of damped sinusoids at a few frequencies are not (Phase
#' 4 study, 2026-10-21). With antisymmetric coupling `A + A'` is twice the
#' negative diagonal, so the mother is always stable; the daughter's unequal
#' edges break the antisymmetry slightly.
#'
#' Pairs are sampled without replacement, first for the shared edges, then
#' the mother-only, then the daughter-only; the two positions of a pair share
#' one magnitude (the first-power coefficient), so a pair is an exact rotation
#' and each cell's off-diagonal coupling is antisymmetric. Of the shared pairs (and the odd
#' single, if any), `round(frac_unequal * number of units)` are unequal
#' (chosen uniformly); both positions of an unequal pair are scaled in the
#' same cell by one factor drawn in `unequal_factor`, so that each cell's
#' coupling stays antisymmetric and stable. Each position's magnitude is drawn
#' uniformly in `edge_magnitude` and, for each of its `q` coefficients,
#' divided by `factorial(m)` for the `m`-th power. On an unequal shared edge
#' the scaffold holds the smaller-magnitude cell's coefficients, so that the
#' larger cell's deviation has the scaffold's sign. The diagonal is `-runif(self_decay)` in
#' both cells. Intercepts are `rnorm(0, intercept_sd)` per equation, shared by
#' both cells unless `bool_shared_intercept = FALSE`.
#'
#' Only reproducibility is promised for the random stream: the same
#' `seed_number` gives identical truths, and `seed_number = NULL` consumes the
#' caller's stream deterministically and leaves it advanced. The
#' number of draws is an implementation detail.
#'
#' @param p number of proteins (at least 2 when any off-diagonal edge is
#' requested).
#' @param q degree of the polynomial basis: each edge carries `q` coefficients,
#' one per power `x, ..., x^q`.
#' @param num_shared number of off-diagonal edges present in both cells.
#' @param num_mother_only number of off-diagonal edges present only in the
#' mother.
#' @param num_daughter_only number of off-diagonal edges present only in the
#' daughter. The three counts must sum to at most `p (p - 1)`.
#' @param bool_shared_intercept if `TRUE`, the daughter has the mother's
#' intercepts; if `FALSE`, her own draw.
#' @param edge_magnitude range of the uniform draw of each coefficient's
#' magnitude; see Details for the choice of the default.
#' @param frac_unequal fraction of the shared edges whose coefficients differ
#' between the cells, in `[0, 1]`.
#' @param intercept_sd standard deviation of the intercept draws; `0` gives
#' zero intercepts.
#' @param self_decay range of the uniform draw of the self-degradation rates;
#' the diagonal of both cells is minus that draw.
#' @param unequal_factor range of the factor between the two cells'
#' coefficients on an unequal shared edge.
#' @param seed_number seed set before the first draw; `NULL` leaves the
#' stream to the caller.
#' @returns a list with `d_delta`, `d_intercept` (length `p`), `d_network`
#' (`p x p` logical, row = equation, column = regulator), `d_ode_pars`
#' (`p x pq`, column `(k - 1) q + m` is the coefficient of regulator `k`'s
#' `m`-th power), `edge_type` (`p x p` character in `absent`, `shared_equal`,
#' `shared_unequal`, `mother_only`, `daughter_only`, `self`), `m_delta`,
#' `m_intercept`, `m_network`, `m_ode_pars`, `p`, `q`, `s_network` (`m_network
#' & d_network`) and `s_pars`, with `m_ode_pars == s_pars + m_delta` and
#' `d_ode_pars == s_pars + d_delta` exactly.
#' @examples
#' truth <- sim_md_network(p = 4, num_shared = 4, num_mother_only = 1,
#'                         num_daughter_only = 1, seed_number = 7)
#' truth$edge_type
#' round(truth$m_ode_pars, 2)
#' round(truth$d_ode_pars, 2)
#' # the daughter's deviation from the scaffold
#' round(truth$d_delta, 2)
#' @export
sim_md_network <- function(p,
                           q = 1,
                           num_shared,
                           num_mother_only,
                           num_daughter_only,
                           bool_shared_intercept = TRUE,
                           edge_magnitude = c(4*pi, 10*pi),
                           frac_unequal = 0.5,
                           intercept_sd = 0,
                           self_decay = c(0.5, 2),
                           unequal_factor = c(0.5, 2),
                           seed_number = NULL){
  stopifnot(length(p) == 1, p >= 1, p == round(p),
            length(q) == 1, length(edge_magnitude) == 2, length(unequal_factor) == 2,
            length(self_decay) == 2, intercept_sd >= 0)
  if(q < 1 || q != round(q)) stop("`q` must be a positive integer, not ", q)
  for(arg in c("num_shared", "num_mother_only", "num_daughter_only")){
    val <- get(arg)
    if(length(val) != 1 || is.na(val) || val < 0 || val != round(val)){
      stop("`", arg, "` must be a single non-negative integer, not ", val)
    }
  }
  if(frac_unequal < 0 || frac_unequal > 1){
    stop("`frac_unequal` must lie in [0, 1], not ", frac_unequal)
  }
  num_total <- num_shared + num_mother_only + num_daughter_only
  if(num_total > p*(p - 1)){
    stop("`num_shared` + `num_mother_only` + `num_daughter_only` = ", num_total,
         " exceeds the number of off-diagonal positions p(p - 1) = ", p*(p - 1))
  }
  if(!is.null(seed_number)) set.seed(seed_number)

  # positions as unordered pairs {j < k} in a random order; a full pair takes
  # both orientations, a single one orientation (two singles of different
  # types may share a pair). `sample(x)` on a length-1 vector would permute
  # seq_len(x), so pairs are drawn by index.
  pair_mat <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  pair_order <- if(nrow(pair_mat) > 0) sample(nrow(pair_mat), nrow(pair_mat)) else integer(0)
  next_pair <- 1L
  single_slot <- NULL
  counts <- c(shared = num_shared, mother_only = num_mother_only, daughter_only = num_daughter_only)
  edge_type <- matrix("absent", p, p)
  diag(edge_type) <- "self"
  # `pair_id` records which unordered pair a position belongs to and
  # `orientation` its sign (+1 for (j, k), -1 for (k, j)) so that paired
  # positions get opposite signs
  sign_mat <- matrix(0, p, p)
  pair_magnitude <- matrix(NA_real_, p, p)
  position_list <- list(shared = integer(0), mother_only = integer(0), daughter_only = integer(0))
  shared_unit_list <- list()
  for(type in names(counts)){
    num_pairs <- counts[[type]] %/% 2
    for(i in seq_len(num_pairs)){
      jk <- pair_mat[pair_order[next_pair], ]
      next_pair <- next_pair + 1L
      pos <- c(jk[1] + (jk[2] - 1)*p, jk[2] + (jk[1] - 1)*p)
      orient <- sample(c(-1, 1), 1)
      sign_mat[pos] <- c(orient, -orient)
      # one magnitude for both orientations: the pair is then an exact rotation
      # and each cell's off-diagonal coupling stays antisymmetric (stable)
      pair_magnitude[pos] <- stats::runif(1, edge_magnitude[1], edge_magnitude[2])
      position_list[[type]] <- c(position_list[[type]], pos)
      if(type == "shared") shared_unit_list[[length(shared_unit_list) + 1]] <- pos
    }
    if(counts[[type]] %% 2 == 1){
      if(is.null(single_slot)){
        jk <- pair_mat[pair_order[next_pair], ]
        next_pair <- next_pair + 1L
        pos <- jk[1] + (jk[2] - 1)*p
        single_slot <- jk[2] + (jk[1] - 1)*p
      } else {
        pos <- single_slot
        single_slot <- NULL
      }
      sign_mat[pos] <- sample(c(-1, 1), 1)
      position_list[[type]] <- c(position_list[[type]], pos)
      if(type == "shared") shared_unit_list[[length(shared_unit_list) + 1]] <- pos
    }
  }
  shared_idx <- position_list$shared
  # unequal shared edges are chosen per pair (both orientations scaled by one
  # factor in one cell), so that each cell's coupling stays antisymmetric and
  # both cells stay stable; an odd single is chosen on its own
  shared_units <- shared_unit_list
  num_unequal <- round(frac_unequal*length(shared_units))
  unequal_units <- shared_units[sample(length(shared_units), num_unequal)]
  unequal_idx <- unlist(unequal_units)
  edge_type[shared_idx] <- "shared_equal"
  edge_type[unequal_idx] <- "shared_unequal"
  unit_factor <- stats::setNames(rep(NA_real_, p*p), seq_len(p*p))
  unit_larger <- stats::setNames(rep(NA, p*p), seq_len(p*p))
  for(unit in unequal_units){
    unit_factor[unit] <- stats::runif(1, unequal_factor[1], unequal_factor[2])
    unit_larger[unit] <- sample(c(TRUE, FALSE), 1)
  }
  edge_type[position_list$mother_only] <- "mother_only"
  edge_type[position_list$daughter_only] <- "daughter_only"
  chosen <- c(position_list$shared, position_list$mother_only, position_list$daughter_only)

  s_pars <- matrix(0, p, p*q)
  m_delta <- matrix(0, p, p*q)
  d_delta <- matrix(0, p, p*q)
  block_cols <- function(k) (k - 1)*q + seq_len(q)
  row_of <- row(matrix(0, p, p))
  col_of <- col(matrix(0, p, p))
  for(idx in chosen){
    j <- row_of[idx]
    k <- col_of[idx]
    cols <- block_cols(k)
    magnitude <- if(is.na(pair_magnitude[idx])) stats::runif(1, edge_magnitude[1], edge_magnitude[2]) else pair_magnitude[idx]
    higher <- if(q > 1) stats::runif(q - 1, edge_magnitude[1], edge_magnitude[2]) else numeric(0)
    coef_vec <- sign_mat[idx]*c(magnitude, higher)/factorial(seq_len(q))
    type <- edge_type[idx]
    if(type == "shared_equal"){
      s_pars[j, cols] <- coef_vec
    } else if(type == "shared_unequal"){
      factor_val <- unit_factor[[idx]]
      bool_mother_larger <- unit_larger[[idx]]
      small_vec <- coef_vec*min(1, factor_val)
      large_vec <- coef_vec*max(1, factor_val)
      s_pars[j, cols] <- small_vec
      if(bool_mother_larger){
        m_delta[j, cols] <- large_vec - small_vec
      } else {
        d_delta[j, cols] <- large_vec - small_vec
      }
    } else if(type == "mother_only"){
      m_delta[j, cols] <- coef_vec
    } else {
      d_delta[j, cols] <- coef_vec
    }
  }
  decay_vec <- -stats::runif(p, self_decay[1], self_decay[2])
  for(k in seq_len(p)){
    s_pars[k, block_cols(k)] <- decay_vec[k]/factorial(seq_len(q))
  }
  m_intercept <- stats::rnorm(p, 0, intercept_sd)
  d_intercept <- if(bool_shared_intercept) m_intercept else stats::rnorm(p, 0, intercept_sd)

  m_ode_pars <- s_pars + m_delta
  d_ode_pars <- s_pars + d_delta
  m_network <- .network_of(m_ode_pars, q)
  d_network <- .network_of(d_ode_pars, q)

  stopifnot(all(diag(m_network)), all(diag(d_network)),
            sum(edge_type != "absent" & edge_type != "self") == num_total)
  list(d_delta = d_delta,
       d_intercept = d_intercept,
       d_network = d_network,
       d_ode_pars = d_ode_pars,
       edge_type = edge_type,
       m_delta = m_delta,
       m_intercept = m_intercept,
       m_network = m_network,
       m_ode_pars = m_ode_pars,
       p = p,
       q = q,
       s_network = m_network & d_network,
       s_pars = s_pars)
}

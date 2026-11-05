# The daughter editor of the inheritance bank: rewire_network() and
# md_truth_from_pars() (R/rewire_network.R, package 0.6.0; design in
# the inheritance study's design notes, section 3.2).

.mother_pars <- function(p = 6, q = 1, num_units = 4, seed = 1){
  sim_md_network(p = p, q = q, num_shared = 2*num_units, num_mother_only = 0, num_daughter_only = 0,
                 frac_unequal = 0, seed_number = seed)$m_ode_pars
}

.off_adjacency <- function(pars, q){
  p <- nrow(pars)
  adj <- matrix(FALSE, p, p)
  for(k in seq_len(p)) adj[, k] <- apply(pars[, (k - 1)*q + seq_len(q), drop = FALSE] != 0, 1, any)
  diag(adj) <- FALSE
  adj
}

## RW1: the clone corner. With every unit kept, none reweighted and none
## flipped, the editor must return its input bit for bit, so that a pair built
## from it has a zero daughter deviation (the anchor every definition must put
## at its maximum).
test_that("RW1: keep_prob = 1 with no reweighting returns the input unchanged", {
  grid <- expand.grid(q = c(1, 3), num_units = c(3, 5), seed = 1:3)
  for(r in seq_len(nrow(grid))){
    g <- grid[r, ]
    label <- paste(names(g), unlist(g), collapse = " ")
    pars <- .mother_pars(q = g$q, num_units = g$num_units, seed = g$seed)
    out <- rewire_network(pars, keep_prob = 1, unequal_prob = 0, flip_prob = 0, seed_number = g$seed)
    expect_identical(out$pars, pars, info = label)
    expect_equal(nrow(out$units_df), g$num_units, info = label)
    expect_true(all(out$units_df$status == "kept_equal"), info = label)
    expect_true(all(out$units_df$factor == 1), info = label)
    expect_true(all(out$units_df$kind == "pair"), info = label)
  }
})

## RW2: density preservation and the pair structure. Every dropped unit is
## replaced by one of the same kind, so the number of off-diagonal edges is
## unchanged, every pair stays a pair (both orientations or neither), and the
## diagonal is never touched. Swept over keep_prob including the two corners.
test_that("RW2: the off-diagonal density, the pair structure and the diagonal are preserved", {
  grid <- expand.grid(q = c(1, 3), keep_prob = c(0, 0.3, 0.7), num_units = c(3, 5), seed = 1:4)
  for(r in seq_len(nrow(grid))){
    g <- grid[r, ]
    label <- paste(names(g), unlist(g), collapse = " ")
    pars <- .mother_pars(q = g$q, num_units = g$num_units, seed = g$seed)
    out <- rewire_network(pars, keep_prob = g$keep_prob, unequal_prob = 0.5, flip_prob = 0.2,
                          seed_number = 100 + g$seed)
    adj_in <- .off_adjacency(pars, g$q)
    adj_out <- .off_adjacency(out$pars, g$q)
    expect_equal(sum(adj_out), sum(adj_in), info = label)
    expect_equal(sum(adj_out), 2*g$num_units, info = label)
    expect_true(all(adj_out == t(adj_out)), info = label)
    diag_cols <- (seq_len(nrow(pars)) - 1)*g$q + 1
    for(m in seq_len(g$q)){
      expect_identical(out$pars[cbind(seq_len(nrow(pars)), diag_cols + m - 1)],
                       pars[cbind(seq_len(nrow(pars)), diag_cols + m - 1)], info = label)
    }
    num_dropped <- sum(out$units_df$status == "dropped")
    expect_equal(sum(out$units_df$status == "added"), num_dropped, info = label)
    if(g$keep_prob == 0) expect_equal(num_dropped, g$num_units, info = label)
  }
})

## RW3: a kept unit is the input times its factor on both orientations and
## every power, a flipped unit has a negative factor in the stated range, and
## a kept or added pair stays an exact rotation at the first power (opposite
## coefficients), which is what keeps the daughter stable at q = 1.
test_that("RW3: kept units are scaled by one factor per unit and pairs stay antisymmetric", {
  grid <- expand.grid(q = c(1, 3), seed = 1:6)
  for(r in seq_len(nrow(grid))){
    g <- grid[r, ]
    label <- paste(names(g), unlist(g), collapse = " ")
    pars <- .mother_pars(q = g$q, num_units = 5, seed = g$seed)
    out <- rewire_network(pars, keep_prob = 0.6, unequal_prob = 0.5, flip_prob = 0.5,
                          unequal_factor = c(1/4, 4), seed_number = 200 + g$seed)
    units_df <- out$units_df
    for(u in seq_len(nrow(units_df))){
      unit <- units_df[u, ]
      j <- unit$j
      k <- unit$k
      cols_k <- (k - 1)*g$q + seq_len(g$q)
      cols_j <- (j - 1)*g$q + seq_len(g$q)
      if(unit$status %in% c("kept_equal", "kept_unequal", "kept_flipped")){
        expect_equal(out$pars[j, cols_k], unit$factor*pars[j, cols_k], info = paste(label, u))
        expect_equal(out$pars[k, cols_j], unit$factor*pars[k, cols_j], info = paste(label, u))
        if(unit$status == "kept_flipped") expect_true(unit$factor < 0 && abs(unit$factor) >= 1/4 && abs(unit$factor) <= 4,
                                                      info = paste(label, u))
        if(unit$status == "kept_unequal") expect_true(unit$factor > 0 && unit$factor != 1, info = paste(label, u))
        if(unit$status == "kept_equal") expect_equal(unit$factor, 1, info = paste(label, u))
      }
      if(unit$status != "dropped"){
        expect_equal(out$pars[j, cols_k[1]], -out$pars[k, cols_j[1]], info = paste(label, u))
      }
    }
  }
})

## RW4: keep_prob = 0 is an independent draw, not a disjoint one. A dropped
## position is available again for a replacement, so over many seeds the
## rewired network overlaps the input at some positions and not at others.
test_that("RW4: at keep_prob = 0 the replacements overlap the input by chance", {
  pars <- .mother_pars(q = 1, num_units = 4, seed = 2)
  adj_in <- .off_adjacency(pars, 1)
  overlap_vec <- vapply(1:60, function(s){
    out <- rewire_network(pars, keep_prob = 0, seed_number = s)
    sum(adj_in & .off_adjacency(out$pars, 1))
  }, numeric(1))
  expect_true(any(overlap_vec > 0))
  expect_true(any(overlap_vec < sum(adj_in)))
  # the expected overlap of two independent draws of 4 unordered pairs among 15
  # is 4 * 4 / 15 pairs; two edges each
  expect_true(mean(overlap_vec) > 0.5 && mean(overlap_vec) < 5)
})

## RW7: replace_prob < 1 thins the daughter. At replace_prob = 0 every dropped
## unit vanishes, so the output is a subset of the input (mother_retained and
## Jaccard then differ, which the bank's density knob relies on); at 1 the
## density is preserved (RW2); in between the number of "added" rows is at
## most the number of "dropped" rows.
test_that("RW7: replace_prob thins the daughter instead of rewiring her", {
  grid <- expand.grid(q = c(1, 3), seed = 1:5)
  for(r in seq_len(nrow(grid))){
    g <- grid[r, ]
    label <- paste(names(g), unlist(g), collapse = " ")
    pars <- .mother_pars(q = g$q, num_units = 5, seed = g$seed)
    adj_in <- .off_adjacency(pars, g$q)
    out0 <- rewire_network(pars, keep_prob = 0.4, replace_prob = 0, seed_number = 400 + g$seed)
    adj0 <- .off_adjacency(out0$pars, g$q)
    expect_true(all(adj0[adj0] == adj_in[adj0]), info = label)
    expect_true(sum(adj0) < sum(adj_in), info = label)
    expect_false(any(out0$units_df$status == "added"), info = label)
    out_half <- rewire_network(pars, keep_prob = 0.2, replace_prob = 0.5, seed_number = 400 + g$seed)
    expect_true(sum(out_half$units_df$status == "added") <= sum(out_half$units_df$status == "dropped"), info = label)
    expect_true(sum(.off_adjacency(out_half$pars, g$q)) <= sum(adj_in), info = label)
  }
  pars <- .mother_pars(q = 1, num_units = 3, seed = 1)
  expect_error(rewire_network(pars, keep_prob = 0.5, replace_prob = 1.5), "replace_prob")
})

## RW5: reproducibility by seed_number, and seed_number = NULL leaves the
## stream to the caller (pitfalls: a NULL seed must consume the stream, so two
## consecutive calls differ, and the same ambient seed reproduces them).
test_that("RW5: the same seed reproduces the draw and a NULL seed uses the caller's stream", {
  pars <- .mother_pars(q = 3, num_units = 4, seed = 5)
  a <- rewire_network(pars, keep_prob = 0.5, unequal_prob = 0.5, flip_prob = 0.3, seed_number = 7)
  b <- rewire_network(pars, keep_prob = 0.5, unequal_prob = 0.5, flip_prob = 0.3, seed_number = 7)
  c <- rewire_network(pars, keep_prob = 0.5, unequal_prob = 0.5, flip_prob = 0.3, seed_number = 8)
  expect_identical(a, b)
  expect_false(identical(a$pars, c$pars))
  set.seed(3)
  d1 <- rewire_network(pars, keep_prob = 0.5, unequal_prob = 0.5, flip_prob = 0.3, seed_number = NULL)
  d2 <- rewire_network(pars, keep_prob = 0.5, unequal_prob = 0.5, flip_prob = 0.3, seed_number = NULL)
  set.seed(3)
  e1 <- rewire_network(pars, keep_prob = 0.5, unequal_prob = 0.5, flip_prob = 0.3, seed_number = NULL)
  expect_identical(d1, e1)
  expect_false(identical(d1$pars, d2$pars))
})

## RW6: the preamble names the offending argument.
test_that("RW6: bad inputs are refused by name", {
  pars <- .mother_pars(q = 1, num_units = 3, seed = 1)
  expect_error(rewire_network(pars, keep_prob = 1.5), "`keep_prob`")
  expect_error(rewire_network(pars, keep_prob = 0.5, flip_prob = -0.1), "`flip_prob`")
  expect_error(rewire_network(pars, keep_prob = 0.5, unequal_prob = 2), "`unequal_prob`")
  expect_error(rewire_network(pars[, 1:5], keep_prob = 0.5), "p x pq")
  expect_error(rewire_network(pars, keep_prob = 0.5, unequal_factor = c(-1, 2)))
})

## MT1: md_truth_from_pars() recovers sim_md_network()'s own decomposition.
## The generator's scaffold on a shared-unequal edge is the smaller-magnitude
## cell, which is the convention the assembler uses, so the round trip must be
## exact for every element (a drift here would make the two truths score
## differently under score_md_networks()).
test_that("MT1: a sim_md_network() truth round-trips through md_truth_from_pars()", {
  grid <- expand.grid(q = c(1, 3), seed = 1:5)
  for(r in seq_len(nrow(grid))){
    g <- grid[r, ]
    label <- paste(names(g), unlist(g), collapse = " ")
    truth <- sim_md_network(p = 6, q = g$q, num_shared = 6, num_mother_only = 2, num_daughter_only = 2,
                            frac_unequal = 0.5, seed_number = g$seed)
    back <- md_truth_from_pars(truth$m_ode_pars, truth$d_ode_pars, m_intercept = truth$m_intercept,
                               d_intercept = truth$d_intercept)
    for(nm in c("s_pars", "m_delta", "d_delta", "edge_type", "m_network", "d_network", "s_network", "p", "q",
                "m_intercept", "d_intercept")){
      expect_equal(back[[nm]], truth[[nm]], info = paste(label, nm))
    }
  }
})

## MT2: the assembler on rewired pairs. The identities hold exactly, the edge
## type of every position matches the two adjacencies, a sign-flipped shared
## edge has a zero scaffold with both deltas full, a same-sign unequal edge
## has the smaller cell in the scaffold, and a clone has no daughter
## deviation at all.
test_that("MT2: rewired pairs decompose consistently, flips have a zero scaffold, clones a zero deviation", {
  grid <- expand.grid(q = c(1, 3), keep_prob = c(0.2, 0.7, 1), seed = 1:4)
  for(r in seq_len(nrow(grid))){
    g <- grid[r, ]
    label <- paste(names(g), unlist(g), collapse = " ")
    pars <- .mother_pars(q = g$q, num_units = 4, seed = g$seed)
    out <- rewire_network(pars, keep_prob = g$keep_prob, unequal_prob = 0.5, flip_prob = 0.4,
                          seed_number = 300 + g$seed)
    pair <- md_truth_from_pars(pars, out$pars)
    expect_equal(pair$m_ode_pars, pair$s_pars + pair$m_delta, info = label)
    expect_equal(pair$d_ode_pars, pair$s_pars + pair$d_delta, info = label)
    adj_m <- .off_adjacency(pars, g$q)
    adj_d <- .off_adjacency(out$pars, g$q)
    off <- row(adj_m) != col(adj_m)
    expect_true(all(diag(pair$edge_type) == "self"), info = label)
    expect_identical(pair$edge_type[off] %in% c("shared_equal", "shared_unequal"), (adj_m & adj_d)[off], info = label)
    expect_identical(pair$edge_type[off] == "mother_only", (adj_m & !adj_d)[off], info = label)
    expect_identical(pair$edge_type[off] == "daughter_only", (!adj_m & adj_d)[off], info = label)
    expect_identical(pair$edge_type[off] == "absent", (!adj_m & !adj_d)[off], info = label)
    for(u in seq_len(nrow(out$units_df))){
      unit <- out$units_df[u, ]
      if(!grepl("^kept", unit$status)) next
      cols_k <- (unit$k - 1)*g$q + seq_len(g$q)
      s_blk <- pair$s_pars[unit$j, cols_k]
      if(unit$status == "kept_flipped"){
        expect_true(all(s_blk == 0), info = paste(label, u))
        expect_true(pair$flip_network[unit$j, unit$k] || pair$flip_network[unit$k, unit$j], info = paste(label, u))
        expect_equal(pair$m_delta[unit$j, cols_k], pars[unit$j, cols_k], info = paste(label, u))
        expect_equal(pair$d_delta[unit$j, cols_k], out$pars[unit$j, cols_k], info = paste(label, u))
      } else if(unit$status == "kept_unequal"){
        expect_false(pair$flip_network[unit$j, unit$k] || pair$flip_network[unit$k, unit$j], info = paste(label, u))
        smaller <- if(unit$factor < 1) out$pars[unit$j, cols_k] else pars[unit$j, cols_k]
        expect_equal(s_blk, smaller, info = paste(label, u))
        expect_equal(pair$edge_type[unit$j, unit$k], "shared_unequal", info = paste(label, u))
      } else {
        expect_equal(pair$edge_type[unit$j, unit$k], "shared_equal", info = paste(label, u))
      }
    }
  }
  pars <- .mother_pars(q = 1, num_units = 4, seed = 9)
  clone <- md_truth_from_pars(pars, rewire_network(pars, keep_prob = 1, seed_number = 1)$pars)
  expect_true(all(clone$d_delta == 0))
  expect_true(all(clone$edge_type %in% c("absent", "shared_equal", "self")))
  expect_false(any(clone$flip_network))
})

## MT3: the assembled truth is accepted by build_md_dataset(truth = ) and by
## score_md_networks(), which is the whole reason it has the generator's shape:
## the bank is built through build_md_dataset(), the truth scores itself at
## MCC 1 and a clone's true inheritance (the flow-map similarity) is exactly 1.
test_that("MT3: an assembled truth builds a dataset and scores", {
  pars <- .mother_pars(q = 1, num_units = 4, seed = 11)
  out <- rewire_network(pars, keep_prob = 0.5, unequal_prob = 0.5, flip_prob = 0.3, seed_number = 12)
  pair <- md_truth_from_pars(pars, out$pars)
  ds <- build_md_dataset(seed = 13, num_d = 24, noise_level_ratio = 0.5, truth = pair,
                         obs_time = seq(0, 1, length.out = 48), steps_per_frame = 100)
  expect_equal(ds$truth$edge_type, pair$edge_type)
  expect_true(is.finite(ds$truth$inheritance))
  expect_lt(ds$truth$inheritance, 1)
  self_fit <- list(m_network_est = pair$m_network, d_network_est = pair$d_network, s_network_est = pair$s_network)
  sc <- score_md_networks(self_fit, ds$truth)
  expect_equal(unname(sc$mcc[["D"]]), 1)
  expect_equal(unname(sc$mcc[["M"]]), 1)
  clone <- md_truth_from_pars(pars, pars)
  ds_clone <- build_md_dataset(seed = 13, num_d = 24, noise_level_ratio = 0.5, truth = clone,
                               obs_time = seq(0, 1, length.out = 48), steps_per_frame = 100)
  expect_identical(ds_clone$truth$inheritance, 1)
})

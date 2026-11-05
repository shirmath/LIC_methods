# Plan section 15: map_to() (U1). U2 (the retired names stay out of the
# namespace) and P1 (the fork's S3 methods) were removed on 2026-10-23: the fork now lives in tests/testthat/helper-grplasso_fork.R
# as the G1/F6(a) cross-check and is not tested itself.

## U1: map_to() keeps its nearest-neighbour contract as a snapshot written as
## code; its one remaining caller is make_daughter_grid() (B4).
test_that("U1: map_to() maps grid points to their positions, midpoints to the first neighbour, and clamps the ends", {
  grid <- seq(0, 1, by = 0.25)
  expect_identical(map_to(c(0, 0.25, 0.5, 1), grid), c(1L, 2L, 3L, 5L))
  expect_identical(map_to(c(0.125, 0.375), grid), c(1L, 2L))
  expect_identical(map_to(c(-1, 5), grid), c(1L, 5L))
  expect_identical(map_to(0.3, grid), 2L)
})

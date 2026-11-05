#' Nearest grid position of each value
#'
#' For each element of `vec_from`, the index of the closest value in `vec_to`
#' (ties go to the first). A nearest-neighbour lookup that the estimator no
#' longer uses for frames (those are read by exact index); its remaining caller
#' is `make_daughter_grid()`.
#'
#' @param vec_from values to look up.
#' @param vec_to the grid.
#' @returns an integer vector of positions in `vec_to`.
#' @export
map_to <- function(vec_from, vec_to){
  n <- length(vec_from)
  res <- sapply(seq_len(n), function(i){
    which.min(abs(vec_to - vec_from[i]))
  })
  res
}

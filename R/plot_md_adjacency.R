#' Draw a network as an adjacency image
#'
#' A `p x p` grid with equation `j` on the rows (equation 1 at the top) and
#' regulator `k` on the columns, a filled cell where the edge `j <- k` is
#' present. Works on a true network (`truth$m_network`, a coefficient matrix
#' such as `truth$m_ode_pars` at `q = 1`) and on a fitted one
#' (`fit$d_network_est` at one `lambda`): any nonzero entry counts as an edge.
#' The function draws one panel and sets no `par(mfrow)`, so the caller
#' arranges several panels. The drawn matrix is returned so that tests check
#' numbers, not pixels.
#'
#' @section The same edges as `plot_md_network()`:
#' A filled cell in row `j`, column `k` of this image is the arrow from node
#' `k` to node `j` in [plot_md_network()]: the image shows presence only, the
#' graph adds each edge's sign and size. The one difference by default is the
#' diagonal. A true network has every self edge `j <- j` (self-degradation),
#' the image draws them, and the graph leaves them out, as
#' [score_md_networks()] leaves them out of every score. With `bool_diagonal
#' = FALSE` the diagonal is blanked and the two plots show exactly the same
#' edge set.
#'
#' @param network a `p x p` logical or numeric matrix, rows = equation `j`,
#' columns = regulator `k`; nonzero (or `TRUE`) is an edge.
#' @param bool_diagonal draw the self edges `j <- j`; `FALSE` blanks the
#' diagonal (and sets it to `FALSE` in the returned matrix).
#' @param col two colours, absent then present.
#' @param main the panel title.
#' @param protein_names axis labels, length `p`; defaults to the row names of
#' `network` or to `1:p`.
#' @param ... further graphical parameters for `graphics::image()`.
#' @returns invisibly, the `p x p` logical matrix drawn (`network != 0`, with
#' the diagonal `FALSE` when `bool_diagonal = FALSE`), with `protein_names`
#' as row and column names.
#' @examples
#' truth <- sim_md_network(p = 4, num_shared = 4, num_mother_only = 1,
#'                         num_daughter_only = 1, seed_number = 3)
#' old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
#' plot_md_adjacency(truth$m_network, main = "mother")
#' drawn <- plot_md_adjacency(truth$d_network, main = "daughter")
#' par(old_par)
#' drawn
#' @export
plot_md_adjacency <- function(network,
                              bool_diagonal = TRUE,
                              col = c("white", "grey30"),
                              main = "",
                              protein_names = NULL,
                              ...){
  if(!is.matrix(network) || nrow(network) != ncol(network)){
    stop("`network` must be a square matrix (p x p); got ",
         paste(dim(as.matrix(network)), collapse = " x "))
  }
  if(!is.logical(network) && !is.numeric(network)){
    stop("`network` must be a logical or numeric matrix")
  }
  p <- nrow(network)
  if(is.null(protein_names)) protein_names <- rownames(network)
  if(is.null(protein_names)) protein_names <- as.character(seq_len(p))
  if(length(protein_names) != p){
    stop("`protein_names` has length ", length(protein_names), " but the network has p = ", p)
  }
  stopifnot(length(col) == 2, length(bool_diagonal) == 1)
  if(is.na(bool_diagonal) || !is.logical(bool_diagonal)) stop("`bool_diagonal` must be TRUE or FALSE")

  drawn <- network != 0
  drawn[is.na(drawn)] <- FALSE
  if(!bool_diagonal) diag(drawn) <- FALSE
  dimnames(drawn) <- list(protein_names, protein_names)
  # image() puts row 1 of z at the bottom, so the rows are reversed to put
  # equation 1 at the top; zlim is fixed so an empty network is all `col[1]`
  z_mat <- t(drawn[p:1, , drop = FALSE] + 0)
  graphics::image(seq_len(p), seq_len(p), z_mat, col = col, zlim = c(0, 1), axes = FALSE,
                  xlab = "regulator k", ylab = "equation j", main = main, ...)
  graphics::axis(1, seq_len(p), protein_names)
  graphics::axis(2, seq_len(p), rev(protein_names), las = 1)
  graphics::box()
  invisible(drawn)
}

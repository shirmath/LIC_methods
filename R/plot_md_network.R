#' Draw a regulatory network as a graph, with the direction and size of each edge
#'
#' Each protein is a node and each edge `j <- k` an arrow from regulator `k`
#' to equation `j`, drawn with `igraph`. The edge set is `network` when
#' given (for example `fit$d_network_est` at one `lambda`), otherwise every
#' regulator block of `coef` with a nonzero coefficient. The diagonal
#' (self-degradation) is left out by default, as in `score_md_networks()`.
#' The arrow `k -> j` here is the filled cell in row `j`, column `k` of
#' [plot_md_adjacency()] (which draws the diagonal unless `bool_diagonal =
#' FALSE`); the two functions draw the same edge set, this one with each
#' edge's sign and size added.
#'
#' @section How to read the plot:
#' Every drawn edge carries one signed **weight** `w_jk`, the slope of
#' regulator `k`'s contribution to equation `j`, `g_jk(x) = sum_m theta_jkm
#' x^m`, at the cell's **mean state**: `w_jk = g_jk'(xbar_k)`, the `(j, k)`
#' entry of the Jacobian of the fitted field at `xbar`, the column means of
#' `X` (the observations of that cell, `raw$Y_M` or `raw$Y_D`) unless `x_ref`
#' is given. This is the linearization the flow-map inheritance reads. For
#' `q = 1` the weight is the coefficient `theta_jk` itself and neither `X` nor
#' `x_ref` is needed; for `q > 1` `X` is required. The sign of the weight is
#' the **direction** at the mean state: positive is activating (dark green),
#' negative inhibiting (purple); a zero weight is grey. The **line width** is
#' linear in `|w_jk|`, from 0.5 at zero to 4 at `width_max`, which is the
#' largest `|w_jk|` in this plot unless supplied; pass the same `width_max` to
#' two calls (mother and daughter, truth and estimate) to put them on one
#' scale. When `q > 1` and the derivative of `g_jk` changes sign on the
#' observed range of `x_k` in `X` (so the regulator activates at some levels
#' and inhibits at others) the edge is **dashed**: its colour is its direction
#' at the mean state, but not at every level the cell visits. Nodes sit on a
#' circle in protein order by default, so the same protein is at the same
#' place in every plot of the same `p`.
#'
#' @param coef a `p x pq` coefficient matrix, rows = equation `j`, columns =
#' the `q` powers of regulator 1, then of regulator 2, and so on: a fit's
#' `m_ode_est` / `d_ode_est` / `s_ode_est` at one `lambda`, or a truth's
#' `m_ode_pars` / `d_ode_pars` / `s_pars`.
#' @param bool_diagonal draw the self edges `j <- j` as loops.
#' @param layout `"circle"`, or a `p x 2` matrix of node coordinates.
#' @param main the title.
#' @param network `NULL` (every nonzero block is an edge) or a `p x p` logical
#' matrix of the edges to draw.
#' @param protein_names node labels, length `p`; defaults to the row names of
#' `coef` or to `1:p`.
#' @param q the basis degree; `NULL` reads it as `ncol(coef) / nrow(coef)`.
#' @param width_max the `|w_jk|` drawn at the largest width; `NULL` uses the
#' largest weight in this plot.
#' @param X `NULL`, or an `n x p` matrix of this cell's observations (rows =
#' frames, columns = proteins), required when `q > 1`: its column means are
#' the default mean state and its column ranges decide the dashing.
#' @param x_ref `NULL` (the column means of `X`), or a length-`p` state at
#' which the weights are read when `q > 1`; ignored at `q = 1`.
#' @param ... further arguments for `igraph::plot.igraph()`; one that this
#' function also sets (`vertex.size`, `edge.arrow.size`, `edge.curved`,
#' `vertex.color`, ...) overrides the default here.
#' @returns invisibly, a list with `edges` (data frame, one row per drawn edge
#' in equation-major order, columns `from` and `to` (protein names, regulator
#' then equation), `weight`, `sign`, `monotone`, `lwd`, `colour`) and `graph`
#' (the `igraph` object drawn).
#' @examples
#' if(requireNamespace("igraph", quietly = TRUE)){
#'   truth <- sim_md_network(p = 4, num_shared = 4, num_mother_only = 1,
#'                           num_daughter_only = 1, seed_number = 3)
#'   width_max <- max(abs(c(truth$m_ode_pars, truth$d_ode_pars)))
#'   old_par <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
#'   plot_md_network(truth$m_ode_pars, main = "mother", width_max = width_max)
#'   out <- plot_md_network(truth$d_ode_pars, main = "daughter",
#'                          width_max = width_max)
#'   par(old_par)
#'   out$edges
#' }
#' @export
plot_md_network <- function(coef,
                            bool_diagonal = FALSE,
                            layout = "circle",
                            main = "",
                            network = NULL,
                            protein_names = NULL,
                            q = NULL,
                            width_max = NULL,
                            X = NULL,
                            x_ref = NULL,
                            ...){
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("plot_md_network() needs the `igraph` package; install it with install.packages(\"igraph\")")
  }
  coef <- as.matrix(coef)
  p <- nrow(coef)
  if(is.null(q)){
    q <- ncol(coef)/p
    if(q != round(q) || q < 1) stop("`coef` has ", ncol(coef), " columns, not a multiple of its p = ", p, " rows")
  } else if(length(q) != 1 || q < 1 || q != round(q) || ncol(coef) != p*q){
    stop("`q` = ", paste(q, collapse = ","), " does not match `coef` (", p, " x ", ncol(coef), ")")
  }
  if(anyNA(coef)) stop("`coef` has NA entries; the plot needs every coefficient")
  if(is.null(protein_names)) protein_names <- rownames(coef)
  if(is.null(protein_names)) protein_names <- as.character(seq_len(p))
  if(length(protein_names) != p){
    stop("`protein_names` has length ", length(protein_names), " but `coef` has p = ", p, " rows")
  }
  if(anyDuplicated(protein_names)){
    stop("`protein_names` must be distinct (igraph names its vertices by them); duplicated: ",
         paste(unique(protein_names[duplicated(protein_names)]), collapse = ", "))
  }
  if(is.null(network)){
    network <- .network_of(coef, q)
  } else {
    network <- as.matrix(network) != 0
    if(nrow(network) != p || ncol(network) != p) stop("`network` must be a p x p matrix with p = ", p)
    if(anyNA(network)) stop("`network` has NA entries")
  }
  if(!bool_diagonal) diag(network) <- FALSE
  if(q > 1){
    if(is.null(X)) stop("`X` (this cell's n x p observations) is required when q > 1 to set each regulator's range")
    X <- as.matrix(X)
    if(ncol(X) != p) stop("`X` must have p = ", p, " columns (one per protein); it has ", ncol(X))
    if(is.null(x_ref)){
      x_ref <- colMeans(X, na.rm = TRUE)
    } else if(length(x_ref) != p || anyNA(x_ref)){
      stop("`x_ref` must be a length-p state with no NA, p = ", p)
    }
  }
  if(!is.null(width_max) && (length(width_max) != 1 || width_max <= 0)) stop("`width_max` must be one positive number")

  # the edges in equation-major order (row j, then column k within it)
  edge_idx <- which(t(network), arr.ind = TRUE)[, c(2, 1), drop = FALSE]
  num_edges <- nrow(edge_idx)
  weight_vec <- numeric(num_edges)
  monotone_vec <- rep(TRUE, num_edges)
  # the weight of edge j <- k: the slope of g_jk at the mean state x_ref (the
  # Jacobian entry), the coefficient itself at q = 1; monotone when the
  # derivative keeps one sign on the regulator's observed range
  for(e in seq_len(num_edges)){
    j <- edge_idx[e, 1]
    k <- edge_idx[e, 2]
    theta_vec <- coef[j, (k - 1)*q + seq_len(q)]
    if(q == 1){
      weight_vec[e] <- unname(theta_vec)
      next
    }
    x_range <- suppressWarnings(range(X[, k], na.rm = TRUE))
    if(!all(is.finite(x_range))){
      stop("column ", k, " of `X` (", protein_names[k], ") has no finite value, so regulator ", k, "'s range is undefined")
    }
    dg_fun <- function(x) sum(seq_len(q)*theta_vec*x^(seq_len(q) - 1))
    weight_vec[e] <- dg_fun(x_ref[k])
    dg_vec <- vapply(seq(x_range[1], x_range[2], length.out = 200), dg_fun, numeric(1))
    monotone_vec[e] <- all(dg_vec >= 0) || all(dg_vec <= 0)
  }
  edges <- data.frame(from = protein_names[edge_idx[, 2]], to = protein_names[edge_idx[, 1]],
                      weight = weight_vec, sign = sign(weight_vec), monotone = monotone_vec,
                      stringsAsFactors = FALSE)
  if(is.null(width_max)) width_max <- if(num_edges > 0 && any(weight_vec != 0)) max(abs(weight_vec)) else 1
  edges$lwd <- 0.5 + 3.5*abs(edges$weight)/width_max
  colours <- .md_edge_colours()
  edges$colour <- unname(colours[c("neutral", "activating", "inhibiting")[match(edges$sign, c(0, 1, -1))]])

  graph <- igraph::graph_from_data_frame(edges[, c("from", "to"), drop = FALSE], directed = TRUE,
                                         vertices = data.frame(name = protein_names, stringsAsFactors = FALSE))
  if(is.character(layout)){
    if(layout != "circle") stop("`layout` must be \"circle\" or a p x 2 matrix of coordinates")
    layout_mat <- igraph::layout_in_circle(graph)
  } else {
    layout_mat <- as.matrix(layout)
    if(nrow(layout_mat) != p || ncol(layout_mat) != 2) stop("a matrix `layout` must be p x 2")
  }
  legend_text <- if(q == 1) c("activating (w > 0)", "inhibiting (w < 0)") else
    c("activating at the mean state", "inhibiting at the mean state")
  legend_col <- unname(colours[c("activating", "inhibiting")])
  legend_lty <- c(1, 1)
  if(any(!edges$monotone)){
    legend_text <- c(legend_text, "dashed: direction changes over the range")
    legend_col <- c(legend_col, "grey30")
    legend_lty <- c(legend_lty, 2)
  }
  # the layout is rescaled to [-1, 1]; the lower y-limit leaves a band below
  # the circle for the legend, which no node reaches (one row per entry)
  y_low <- -1.15 - 0.17*length(legend_text)
  plot_args <- list(layout = layout_mat, main = main, ylim = c(y_low, 1.05),
                    edge.color = edges$colour, edge.width = edges$lwd,
                    edge.lty = ifelse(edges$monotone, 1, 2), edge.arrow.size = 0.5, edge.curved = 0.2,
                    vertex.color = "grey92", vertex.frame.color = "grey40", vertex.label.color = "black",
                    vertex.size = 28)
  plot_args <- utils::modifyList(plot_args, list(...))
  do.call(igraph::plot.igraph, c(list(graph), plot_args))
  graphics::legend("bottom", legend_text, col = legend_col, lty = legend_lty, lwd = 2, bty = "n", cex = 0.8)
  invisible(list(edges = edges,
                 graph = graph))
}

#' The mother-daughter colour convention
#'
#' `plot_md_pair()` and `plot_md_pca3d()` draw from these ramps, so the
#' convention (mother blue, daughter red) is defined once: `hcl.colors(n,
#' "Blues 3", rev = TRUE)` and `hcl.colors(n, "Reds 3", rev = TRUE)`, light at
#' the first entry and dark at the last. `plot_md_pair()` uses the darkest
#' entry of each three-colour ramp.
#'
#' @param n_M,n_D ramp lengths.
#' @returns a list with `M` and `D`, each a vector of colours.
#' @noRd
.md_palette <- function(n_M,
                        n_D){
  list(M = grDevices::hcl.colors(n_M, "Blues 3", rev = TRUE),
       D = grDevices::hcl.colors(n_D, "Reds 3", rev = TRUE))
}

#' The edge colours of `plot_md_network()`
#'
#' Activating (positive weight) is dark green and inhibiting (negative weight)
#' purple, the two ends of the PRGn diverging palette, so that neither clashes
#' with the mother-blue / daughter-red convention of `.md_palette()`. A zero
#' weight is grey.
#'
#' @returns a named character vector `activating`, `inhibiting`, `neutral`.
#' @noRd
.md_edge_colours <- function(){
  c(activating = "#1B7837", inhibiting = "#762A83", neutral = "grey50")
}

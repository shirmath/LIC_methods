# LICmethods

PRNODE (protein regulatory networks via ordinary differential equations): a
fused group-lasso ODE method that jointly estimates the protein regulatory
networks of a mother yeast cell and one of her daughters from live-cell
time-lapse microscopy, and quantifies how much of the network the daughter
inherits.

## Install

```r
# install.packages("remotes")
remotes::install_github("shirmath/LIC_methods", build_vignettes = TRUE)
```

Optional packages: `glmnet`, `igraph` (network plots), `deSolve` (ODE
simulation), `rmarkdown` and `knitr` (the vignette).

## Start here

```r
library(LICmethods)
vignette("prnode_overview")
```

The vignette shows how to simulate mother-daughter pairs (`build_md_dataset()`),
fit PRNODE to a pair (`fit_prnode()`) and estimate and test the pair's
inheritance (`inheritance_flow()`, `inheritance_test()`).

## Development

Run the tests with `NOT_CRAN=true`:

```r
devtools::test()
```

To enable the large-file guard, run `git config core.hooksPath .githooks` once
after cloning.

## Authors

Shirley Mathur (maintainer) and Kevin Lin. MIT license.

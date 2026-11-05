# LICmethods 0.13.0

First release of this repository. It holds the package, its tests and the
overview vignette.

The overview vignette is a user guide to three tasks: simulating
mother-daughter pairs, fitting PRNODE, and the inheritance permutation test.

State of the method at this version:

- **Estimator.** `make_prnode_design()` builds the two-block stacked design of
  a mother and one daughter (scaffold, mother-deviation and daughter-deviation
  blocks, each cell's intercept and slope profiled out, one common weight);
  `fit_prnode_path()` is the reference group-lasso solver and `fit_prnode()` the
  user-facing wrapper.
- **Tuning.** `tune_prnode()` picks the penalty by the known-variance BIC and
  `tune_prnode_spar()` chooses the smoothing level jointly with it.
- **Inheritance.** `flow_map_similarity()`, `inheritance_flow()` and the
  wrong-mother permutation test `inheritance_test()`. At `q > 1` each cell's
  fitted field is linearized at that cell's own mean state (`mean_state` in the
  design); nothing changes at `q = 1`.
- **Simulation and scoring.** `sim_md_network()`, `sim_ode_trajectory()`,
  `sim_md_trajectory()`, `build_md_dataset()`, `rewire_network()`,
  `md_truth_from_pars()` and `score_md_networks()`.
- **Plots.** `plot_md_pair()`, `plot_md_pca3d()`, `plot_md_adjacency()`,
  `plot_md_network()` and `plot_md_inheritance()`.

# Dynamo API quick reference

Dynamo uses short module aliases. `import dynamo as dyn`, then:

| Alias | Full module | Role |
|-------|-------------|------|
| `dyn.pp` | `dynamo.preprocessing` | QC, normalization, gene selection, PCA |
| `dyn.tl` | `dynamo.tools` | dynamics/velocity, dimension reduction, clustering |
| `dyn.vf` | `dynamo.vectorfield` | vector-field reconstruction, differential geometry, topography |
| `dyn.pd` | `dynamo.prediction` | fate, least-action paths, perturbation |
| `dyn.pl` | `dynamo.plot` | plotting for every stage |
| `dyn.sample_data` | — | bundled datasets (download on first use) |

## Data I/O
- `dyn.read_h5ad(path)`, `dyn.read_loom(path)`, `dyn.read(path)` — load AnnData.
- `adata.write_h5ad(path)` — save.
- `dyn.sample_data.zebrafish()`, `.hematopoiesis()`, `.scNT_seq_neuron_labeling()`,
  `.scEU_seq_rpe1()`, `.scEU_seq_organoid()`, `.DentateGyrus()`, `.pancreatic_endocrinogenesis()`, ...

## Preprocessing (`dyn.pp`)
- `Preprocessor(cell_cycle_score_enable=...)` — pipeline object.
  - `.config_monocle_recipe(adata, n_top_genes=2000)` — configure the monocle recipe.
  - `.preprocess_adata(adata, recipe="monocle"|"seurat"|"sctransform"|"pearson_residuals", tkey=None, experiment_type=None)`.
- Individual steps: `calc_sz_factor`, `normalize`, `select_genes_monocle`,
  `filter_cells_by_outliers`, `filter_genes_by_outliers`, `pca`.
- `cell_cycle_scores(adata)` — cell-cycle phase.

## Tools (`dyn.tl`)
- `dynamics(adata, model="stochastic"|"deterministic"|"auto", est_method=..., tkey=None, group=None, cores=1)`
  — kinetics + RNA velocity. **The core velocity step.**
- `reduceDimension(adata, reduction_method="umap"|"tsne"|"psl", n_pca_components=30, n_neighbors=30)`.
- `cell_velocities(adata, basis="umap", ekey=None, vkey=None)` — project velocity onto embedding.
- `neighbors`, `leiden`, `louvain`, `hdbscan`, `scc` — graph + clustering.
- Labeling recipes (all-in-one): `recipe_kin_data`, `recipe_deg_data`,
  `recipe_one_shot_data`, `recipe_mix_kin_deg_data`.

## Vector field (`dyn.vf`)
- `VectorField(adata, basis="umap", map_topography=False, pot_curl_div=False, cores=1)`
  — reconstruct continuous field → `adata.uns['VecFld_<basis>']`.
- Differential geometry: `jacobian(adata, regulators=, effectors=, basis="pca")`,
  `acceleration`, `curvature`, `divergence`, `curl`, `speed`, `torsion`, `hessian`, `sensitivity`, `laplacian`.
- Gene ranking: `rank_velocity_genes`, `rank_jacobian_genes(return_df=True)`,
  `rank_acceleration_genes`, `rank_curvature_genes`, `rank_divergence_genes`, `rank_sensitivity_genes`.
- Topography: `topography(adata, basis="umap")`, `assign_fixedpoints`, `FixedPoints`.
- Potential landscape: `Potential`, `Ao_pot_map`.

## Prediction (`dyn.pd`)
- `fate(adata, init_cells, basis=None, direction="forward"|"backward"|"both", interpolation_num=250)`.
- `least_action(adata, init_cells, target_cells, basis="pca", n_points=25)` — LAP.
- `perturbation(adata, genes, expression=10, basis="pca", emb_basis="umap")` — in-silico perturbation.
- `KO(adata, KO_genes, basis="pca", emb_basis="umap")` — knockout.
- `rank_perturbation_genes`, `rank_perturbation_cells`, `rank_perturbation_cell_clusters`.
- `state_graph`, `Trajectory`, `GeneTrajectory`, `GeneLeastActionPath`.

## Plotting (`dyn.pl`)
- Embedding/velocity: `streamline_plot`, `cell_wise_vectors`, `grid_vectors`, `umap`, `pca`, `tsne`.
- Phase/kinetics: `phase_portraits`, `dynamics`, `kinetic_curves`, `kinetic_heatmap`, `jacobian_kinetics`.
- Field: `topography`, `plot_fixed_points`, `plot_flow_field`, `plot_nullclines`, `show_landscape`.
- Differential geometry scalars: `acceleration`, `curvature`, `divergence`, `speed`, `jacobian`.
- Prediction: `fate`, `fate_bias`, `least_action`, `lap_min_time`, `state_graph`.

**Save convention used by the toolkit scripts:** call a `dyn.pl.*` function with
`save_show_or_return="return"`, then `plt.savefig(...)` (helper `save_fig` in `_common.py`).

# Differential geometry — extracting regulatory function from the vector field

Once a continuous vector field `f(x)` is reconstructed, dynamo differentiates it
to read out regulatory/dynamical structure. **Do this in the `pca` basis** (the
Jacobian is defined in the high-dim state space), then rank genes.

Script: `differential_geometry.py` (defaults to `--basis pca`).

## Prerequisite
Build the field in PCA first:
```bash
python vector_field.py red.h5ad -o vf_pca.h5ad --basis pca
```
`vector_field.py` will project velocity onto `pca` automatically if only `umap`
was projected earlier.

## The quantities (`--quantities`)
- **jacobian** `∂fᵢ/∂xⱼ` — how a change in regulator *j* changes the velocity of
  effector *i*. Sign = activation (+) / repression (−). The core of regulatory
  inference. Stored in `uns['jacobian_pca']`.
- **acceleration** `df/dt` — where cells speed up / commit. `obsm['acceleration_*']`.
- **curvature** — how sharply trajectories bend (fate-decision regions).
- **divergence** — net expansion (sources) vs contraction (sinks); per-cell scalar
  in `obs['divergence_pca']`.
- **curl** / **speed** — rotation / magnitude of the flow (2D bases).

```bash
python differential_geometry.py vf_pca.h5ad -o dg.h5ad --basis pca \
    --quantities jacobian acceleration curvature divergence
```

## Ranking genes → CSVs in `--out-dir` (default `results/diffgeo`)
- **acceleration / curvature** produce per-gene rankings (`rank_acceleration.csv`,
  `rank_curvature.csv`) — top driver genes by that quantity.
- **jacobian** ranking requires explicit regulators/effectors — pass `--genes`.
  The gene-wise Jacobian only exists when `jacobian()` is given genes:
  ```bash
  python differential_geometry.py vf_pca.h5ad -o dg.h5ad --basis pca \
      --quantities jacobian --genes GATA1 KLF1 SPI1 CEBPA
  ```
  Output: `rank_jacobian_all.csv` (per-regulator top interacting genes).
- **divergence / curl / speed** are per-cell scalars, not gene-ranked — computed
  for plotting/coloring only (`plot.py --kind divergence`).

Skip ranking with `--no-rank`. Limit rows with `--top-n` (default 25).

## Reading the outputs
- Jacobian entry `J[effector, regulator] > 0` → activation; `< 0` → repression.
  `dyn.vf.rank_jacobian_genes` sorts, for each regulator, the effectors it most
  activates/represses — a data-driven GRN edge list.
- High **acceleration** genes drive commitment; high **curvature** genes mark
  where fate decisions bend trajectories.
- **Divergence**: positive = source (progenitor-like), negative = sink (terminal).

## Underlying calls
```python
dyn.vf.jacobian(adata, regulators=genes, effectors=genes, basis="pca")
dyn.vf.acceleration(adata, basis="pca")
dyn.vf.curvature(adata, basis="pca")
dyn.vf.divergence(adata, basis="pca")
dyn.vf.rank_jacobian_genes(adata, return_df=True)   # -> {group: DataFrame}
dyn.vf.rank_acceleration_genes(adata)               # -> uns['rank_acceleration']
dyn.vf.rank_curvature_genes(adata)                  # -> uns['rank_curvature']
```

## Graph potential / vector-field pseudotime — `--quantities potential`

Unlike the others, `potential` is not a derivative of the field in a basis. It
comes from a Hodge decomposition of a **cell-cell transition matrix in `obsp`**,
and yields dynamo's vector-field pseudotime.

```bash
python differential_geometry.py vf.h5ad -o dg.h5ad --quantities potential \
    --transition-key cosine_transition_matrix
python differential_geometry.py dg.h5ad -o dg.h5ad --quantities potential \
    --transition-key fp_transition_rate
```

Writes `obs['potential_<key>']` and `obs['pseudotime_<key>']` (= −potential).
Underlying calls:

```python
from dynamo.tools.graph_operators import build_graph, div, potential
g = build_graph(adata.obsp[key]); pot = potential(g, -div(g))
```

**The sign of `div` depends on the matrix.** dynamo's own tutorials use
`potential(g, -div(g))` for a cosine transition matrix but `potential(g, +div(g))`
for `fp_transition_rate`. Getting it backwards silently inverts pseudotime — an
error that looks like a biological result rather than a bug. `--potential-sign
auto` (default) picks `pos` for `fp_transition*` keys and `neg` otherwise, and
logs the choice; override with `--potential-sign neg|pos`.

Sanity check: the progenitor population should sit at **low** potential / early
pseudotime. If it doesn't, suspect the sign before believing the biology.

## Visualize
```bash
python plot.py dg.h5ad --kind acceleration --basis umap
python plot.py dg.h5ad --kind jacobian --genes GATA1 KLF1          # one regulator->effector
python plot.py dg.h5ad --kind jacobian_heatmap --genes GATA1 KLF1 SPI1  # pairwise matrix
python plot.py dg.h5ad --kind umap --color potential_cosine_transition_matrix
# compare a scalar ACROSS groups — an embedding plot cannot answer "which is highest?"
python plot.py dg.h5ad --kind scalar_by_group --scalar speed_pca --group cell_type
```

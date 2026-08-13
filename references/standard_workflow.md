# Standard workflow — conventional splicing-based velocity

End-to-end for droplet scRNA-seq with **spliced/unspliced** layers (10x + velocyto
/ scVelo loom, etc.). The toolkit scripts wrap each step; this shows the exact
`dyn.*` calls so you can customize beyond the flags.

## 0. Input requirement
`adata.layers["spliced"]` and `adata.layers["unspliced"]` must exist. Confirm with
`python inspect_data.py data.h5ad` (reports `Detected modality: splicing`).

## 1. Preprocess — `preprocess.py`
```python
import dynamo as dyn
adata = dyn.read_h5ad("raw.h5ad")

preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
preprocessor.config_monocle_recipe(adata, n_top_genes=2000)
preprocessor.preprocess_adata(adata, recipe="monocle")
```
Recipes: `monocle` (default), `seurat`, `sctransform`, `pearson_residuals`. This
does QC, size-factor normalization, feature-gene selection, and PCA
(`adata.obsm["X_pca"]`).

## 2. Dynamics / velocity — `dynamics.py`
```python
dyn.tl.dynamics(adata, model="stochastic")   # or "deterministic"
```
Fits per-gene transcription/splicing/degradation rates and writes velocity
(`adata.layers["velocity_S"]`, smoothed expression `M_s`). `model="stochastic"`
uses second moments and is usually more robust; `"deterministic"` is faster.

## 3. Dimension reduction + projection — `reduce_dimensions.py`
```python
dyn.tl.reduceDimension(adata, reduction_method="umap")   # adata.obsm["X_umap"]
dyn.tl.cell_velocities(adata, basis="umap")              # project velocity onto umap
```
`cell_velocities` builds the low-dim velocity used to fit the field. For splicing
data `ekey`/`vkey` auto-resolve to `M_s`/`velocity_S`.

## 4. Vector field — `vector_field.py`
```python
dyn.vf.VectorField(adata, basis="umap", map_topography=True)
```
Reconstructs the continuous field → `adata.uns["VecFld_umap"]`. `map_topography`
adds fixed points for topography plots. **For differential geometry / perturbation,
also build the field in PCA:** `dyn.vf.VectorField(adata, basis="pca")`.

## 5. Visualize
```python
dyn.pl.streamline_plot(adata, basis="umap", color="cell_type", save_show_or_return="return")
dyn.pl.topography(adata, basis="umap", color="cell_type", save_show_or_return="return")
```
Or `python plot.py vf.h5ad --kind streamline --color cell_type`.

## 6. Downstream
See `differential_geometry.md` (regulatory genes) and `prediction.md` (fate, LAP,
perturbation).

## One command
```bash
python run_pipeline.py raw.h5ad -o vf.h5ad --model stochastic --basis umap --diff-geometry
```
Runs steps 1–4 (+ optional differential geometry). Equivalent to the chain above.

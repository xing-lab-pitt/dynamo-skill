# Plotting guide

All figures come from `plot.py` (one `--kind` per call, saved to `--figdir` as
`<kind>.png`) or directly from `dyn.pl.*`. The toolkit convention is to call a
`dyn.pl` function with `save_show_or_return="return"` then save via matplotlib
(helper `save_fig` in `_common.py`), so figures land in `./figures/` headlessly.

## Kinds (`plot.py --kind`)

| Kind | Needs | Notes |
|------|-------|-------|
| `streamline` | `cell_velocities` (any basis) | velocity streamlines on the embedding; `--color obs/var` |
| `umap` | an embedding | plain embedding colored by `--color` |
| `phase` | dynamics + `--genes` | per-gene phase portraits (spliced/unspliced or new/total) |
| `acceleration` | `differential_geometry.py` | scalar field on the embedding |
| `curvature` | `differential_geometry.py` | scalar field |
| `divergence` | `differential_geometry.py` | sources (+) / sinks (−) |
| `jacobian` | Jacobian computed | interaction field |
| `speed` | `dyn.tl.speed(adata, basis=)` first | otherwise raises "speed_… not in .obs" |
| `topography` | `vector_field.py --map-topography` | streamlines + fixed points; without fixed points it raises "'x' values must be equally spaced" |
| `kinetic_heatmap` | dynamics + `--genes` | expression kinetics ordered along the field |

```bash
python plot.py vf.h5ad  --kind streamline   --color cell_type
python plot.py dg.h5ad  --kind acceleration --basis umap
python plot.py dyn.h5ad --kind phase        --genes GATA1 KLF1
python plot.py vf.h5ad  --kind topography   --color cell_type   # build field with --map-topography
```

## Common direct calls (beyond plot.py)
```python
dyn.pl.streamline_plot(adata, basis="umap", color="cell_type", save_show_or_return="return")
dyn.pl.cell_wise_vectors(adata, basis="umap", color="cell_type", save_show_or_return="return")
dyn.pl.grid_vectors(adata, basis="umap", save_show_or_return="return")
dyn.pl.topography(adata, basis="umap", color="cell_type", save_show_or_return="return")
dyn.pl.phase_portraits(adata, genes=["GATA1"], save_show_or_return="return")
dyn.pl.kinetic_heatmap(adata, genes=[...], save_show_or_return="return")
```

## Tips
- **Headless**: `_common.configure_dynamo` sets the matplotlib Agg backend, so
  everything works without a display.
- **Fixed points / topography**: only available if the field was reconstructed
  with `vector_field.py --map-topography` (or `dyn.vf.VectorField(..., map_topography=True)`).
- **Basis**: pass the basis the field/velocity was computed in (`umap` for most
  visuals; scalar differential-geometry fields can be shown on `umap` after being
  computed in `pca`).
- **Color** accepts obs columns (cell type, cluster) or var names (gene expression).

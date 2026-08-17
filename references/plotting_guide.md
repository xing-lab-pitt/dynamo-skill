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
| `jacobian` | Jacobian computed | interaction field; pass `--genes` for a named regulator→effector pair, omit for the overall field |
| `jacobian_heatmap` | Jacobian + `--genes` | pairwise regulator × effector matrix |
| `speed` | `dyn.tl.speed(adata, basis=)` first | otherwise raises "speed_… not in .obs" |
| `topography` | `vector_field.py --map-topography` | streamlines + fixed points; without fixed points it raises "'x' values must be equally spaced". Curate them with `topography.py --list/--keep` |
| `scalar_by_group` | any per-cell scalar in `obs` + `--scalar --group` | violin+box of the scalar per group, sorted by median; `--palette` takes JSON |
| `lap_paths` | `lap_matrix.py --save-paths` + `--paths --pairs` | least-action paths overlaid on the field, colored by action |
| `kinetic_heatmap` | dynamics + `--genes` | expression kinetics ordered along the field |

```bash
python plot.py vf.h5ad  --kind streamline   --color cell_type
python plot.py dg.h5ad  --kind acceleration --basis umap
python plot.py dyn.h5ad --kind phase        --genes GATA1 KLF1
python plot.py vf.h5ad  --kind topography   --color cell_type   # build field with --map-topography
python plot.py dg.h5ad  --kind scalar_by_group --scalar speed_pca --group cell_type
python plot.py vf.h5ad  --kind lap_paths --paths results/lap/laps.pkl \
    --pairs "HSC->Meg" "HSC->Ery" --out-name lap_development
```

**Embedding vs. group plots.** Painting a scalar on the embedding shows *where*
it is high; it cannot support a claim like "Meg is the fastest cell type". That
comparison needs `scalar_by_group`, which prints the per-group medians alongside
the figure so the ordering is stated, not eyeballed.

**`--out-name`.** One `--kind` can produce several figures (three groupings of
`lap_paths`, several scalars). Without `--out-name` they all save to
`<kind>.png` and overwrite each other.

**`--compute-basis` — the usual differential-geometry gotcha.** Scalars are
computed in `pca` but you want to see them on `umap`. `dyn.pl.speed`,
`divergence`, `acceleration` and `curvature` use one `basis` for *both* the obs
key and the embedding, so `--kind speed --basis umap` looks for `speed_umap` and
raises `speed_umap is not existed in .obs`. Say where it was computed:

```bash
python differential_geometry.py vf.h5ad -o dg.h5ad --basis pca --quantities speed
python plot.py dg.h5ad --kind speed --basis umap --compute-basis pca   # correct
```

When the two differ, the scalar is painted on the requested embedding via
`dyn.pl.scatters`. `--kind jacobian` passes it as dynamo's own `j_basis`, which
already separates the two.

**`--kind jacobian_heatmap` draws one panel per cell.** Left unbounded that is
thousands of panels. The default here is the population-averaged matrix
(`average=True` over all cells); pass `--cell-idx 0 1 2` for specific cells, and
add `--average` to average just those.

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

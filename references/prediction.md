# Prediction — fate, least-action paths, in-silico perturbation

Three predictive analyses built on the reconstructed field (`dyn.pd`). All need a
vector field first (`vector_field.py`).

---

## 1. Cell fate — `fate.py`
Integrate cells forward/backward along the field to trace where they go.

```bash
python fate.py vf.h5ad -o fate.h5ad --init-group Cell_type:HSC --basis umap --color Cell_type
```
- `--init-group col:value` or `--init-cells NAME ...` — starting cells.
- `--direction forward|backward|both`; `--interpolation-num` (path resolution).
- Runs in the **umap** basis for visualization. Produces `fate_streamline.png`
  and a topography plot when the field has fixed points
  (`vector_field.py --map-topography`).
- Underlying: `dyn.pd.fate(adata, init_cells=..., basis="umap", direction=...)`.

---

## 2. Least-action path (LAP) — `least_action.py`
The most probable transition/reprogramming route between two states — dynamo's
optimal-control answer to "what path (and which TFs) converts state A → B".

```bash
python least_action.py vf_pca.h5ad -o lap.h5ad \
    --init-group Cell_type:HSC --target-group Cell_type:Meg --basis pca
```
- Run in the **pca** basis.
- **LAP solves one optimization per init cell**, so it is expensive on whole
  populations. `--n-cells` (default **1**) caps each endpoint — the canonical
  "representative cell → representative cell" LAP. Raise it (e.g. `--n-cells 5`)
  or pass explicit `--init-cells`/`--target-cells`.
- **`--endpoint` picks which cell represents a group**: `centroid` (default,
  closest to the group's centroid in the field's basis) or `first` (first in row
  order, the pre-`--endpoint` behaviour). `first` makes the answer depend on row
  ordering, which is why it is no longer the default. `lap_matrix.py` uses the
  same `_common.representative_cell` helper, so a single LAP and the action
  matrix refer to the same cells.
- Output: `least_action_path.png`; the path + action stored in `uns['LAP_pca']`.
- Genes along the path: `dyn.pd.GeneLeastActionPath` + `dyn.pl.kinetic_heatmap`.

---

## 2b. Action matrix across states — `lap_matrix.py`
One LAP is a route; the **matrix of actions over every ordered pair** is the
comparison that makes LAP interpretable — it ranks how hard each conversion is.
In Qiu et al. (Cell 2022) this is the Fig. 7 result: developmental transitions
cost less action than de-differentiation or transdifferentiation.

```bash
python lap_matrix.py vf_pca.h5ad -o lap.h5ad --basis pca \
    --group-col cell_type --source HSC --targets Meg Ery Bas Mon Neu \
    --out-dir results/lap
```
- Outputs `lap_action_matrix.csv` / `.png` (rows = source, cols = target) and,
  when `--source` is given, `lap_forward_vs_reverse.csv` +
  `lap_action_barplot.png` with a `ratio_reverse_over_forward` column. A ratio
  near 1 means that state is **not** acting as an attractor — worth checking
  against `divergence` before trusting it.
- **Endpoints are the cell closest to each group's centroid** (`--endpoint
  centroid`, default), not the first N cells, so the matrix does not depend on
  row order. `--endpoint first` restores `least_action.py`'s behaviour.
- **Cost grows quadratically**: N groups = N·(N−1) optimizations, each a few
  minutes on a few thousand cells. 6 groups (30 pairs) takes hours — start with
  3–4 groups, or a `--source` plus a couple of `--targets`.
- A failing pair is recorded and skipped rather than aborting the matrix; check
  the log for `FAILED` lines before reading the result as complete.
- Low-n groups are flagged (`--min-cells`, default 50) — a terminal state with a
  handful of cells gives a weakly constrained path.

### Keeping the paths, not just the actions — `--save-paths`

By default the matrix keeps only the final action per pair and **discards the
path geometry**, even though solving for it is the expensive part. `--save-paths`
persists each pair's coordinates and per-point action so they can be re-plotted
without re-solving:

```bash
python lap_matrix.py vf_umap.h5ad -o lapm.h5ad --basis umap \
    --group-col cell_type --groups HSC Meg Ery Bas Mon Neu \
    --save-paths results/lap/laps.pkl

# then group them however you like — no re-optimization
python plot.py vf_umap.h5ad --kind lap_paths --paths results/lap/laps.pkl \
    --pairs "HSC->Meg" "HSC->Ery" "HSC->Bas" --out-name lap_development
python plot.py vf_umap.h5ad --kind lap_paths --paths results/lap/laps.pkl \
    --pairs-file transdiff_pairs.txt --out-name lap_transdifferentiation
```

Keys are `"<source>-><target>"`. Use the **umap** basis when the point is the
figure (paths are drawn in the embedding's first two dimensions); use **pca**
when the point is the action values. Paths are colored by action via
`dynamo.plot.utils.map2color`, so a lighter path is a costlier transition.

---

## 3. In-silico perturbation / knockout — `perturbation.py`
Predict how activating, suppressing, or knocking out gene(s) reshapes the flow —
the effect propagates through the Jacobian.

```bash
# continuous perturbation (activate with +, suppress with -)
python perturbation.py dg.h5ad -o pert.h5ad --genes GATA1 --expression 100
python perturbation.py dg.h5ad -o pert.h5ad --genes SPI1  --expression -100
# knockout
python perturbation.py dg.h5ad -o ko.h5ad   --genes GATA1 --ko
```
- Run with the field in **pca** (`--basis pca`); effect projected to **umap**
  (`--emb-basis umap`) for visualization.
- **Continuous perturbation** (`dyn.pd.perturbation`): computes the delta-velocity
  field, ranks the most-responding genes and cells → `perturbation_rank_genes.csv`
  / `perturbation_rank_cells.csv` in `--out-dir`, and plots the perturbed
  streamlines on `umap_perturbation`.
- **Knockout** (`--ko`, `dyn.pd.KO`): re-simulates the field with the gene zeroed,
  writing `X_umap_KO` / `velocity_umap_KO`, and plots streamlines on `umap_KO`.
  KO does not produce a delta-response matrix, so gene/cell ranking is skipped
  (use continuous perturbation for ranked responders).

### Reading perturbation output
- The perturbed streamline plot shows how fate flow redirects (e.g. suppressing a
  lineage TF diverts cells away from that fate).
- `perturbation_rank_genes.csv`: genes whose velocity changes most — downstream
  targets of the perturbed gene.

---

## Basis cheat-sheet
| Analysis | Field basis | Why |
|----------|-------------|-----|
| fate | umap | trajectories drawn on the embedding |
| LAP | pca | action defined in high-dim state space |
| LAP action matrix | pca | same — actions are only comparable within one basis |
| perturbation / KO | pca (+ umap embedding) | Jacobian propagation in PCA, view in umap |

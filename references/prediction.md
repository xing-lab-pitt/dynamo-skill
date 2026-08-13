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
  populations. `--n-cells` (default **1**) caps each endpoint to the first N
  selected cells — the canonical "representative cell → representative cell" LAP.
  Raise it (e.g. `--n-cells 5`) or pass explicit `--init-cells`/`--target-cells`.
- Output: `least_action_path.png`; the path + action stored in `uns['LAP_pca']`.
- Genes along the path: `dyn.pd.GeneLeastActionPath` + `dyn.pl.kinetic_heatmap`.

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
| perturbation / KO | pca (+ umap embedding) | Jacobian propagation in PCA, view in umap |

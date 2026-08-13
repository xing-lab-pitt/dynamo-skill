# Metabolic labeling workflow (scNT-seq, scEU-seq, scSLAM-seq)

Metabolic-labeling scRNA-seq measures **newly transcribed (labeled) vs total**
RNA at known labeling times, giving a more direct velocity than splicing. Dynamo
fits an explicit kinetic model (transcription/degradation, and splicing when
present). Distinguishing feature: **layers `new` and `total`** plus an obs column
holding the **labeling time**.

## 0. Input requirement
- `adata.layers["new"]` and `adata.layers["total"]` (labeled and total counts).
- An obs column with the labeling duration per cell, e.g. `adata.obs["time"]`.
- `inspect_data.py` reports `Detected modality: labeling`.

If your object has the four splice+label layers (`ul, sl, uu, su`), combine them
into `new`/`total` first (see the scEU-seq tutorial), e.g.:
```python
adata.layers["new"]   = adata.layers["ul"] + adata.layers["sl"]
adata.layers["total"] = adata.layers["uu"] + adata.layers["ul"] + adata.layers["su"] + adata.layers["sl"]
```

## Experiment types
- `kin`  — **kinetics** (pulse): labeling time increases; measures synthesis.
- `deg`  — **degradation** (chase): measures decay after labeling.
- `one-shot` — a single labeling time point.
- `mix_std_stm` — mixed steady-state + stimulation designs.

## 1. Preprocess — `preprocess.py --modality labeling`
```python
preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
preprocessor.config_monocle_recipe(adata, n_top_genes=1000)
preprocessor.preprocess_adata(adata, recipe="monocle", tkey="time", experiment_type="kin")
```
CLI: `python preprocess.py labeled.h5ad -o pp.h5ad --modality labeling --tkey time --experiment-type kin`

## 2. Dynamics — `dynamics.py --tkey time`
```python
dyn.tl.dynamics(adata, model="deterministic", tkey="time", est_method="twostep", cores=8)
```
Passing `tkey` switches dynamo to the **kinetic** model (not steady state).
`est_method="twostep"` is the usual choice for scEU/scNT data. Velocity is written
to `velocity_T` with smoothed expression `M_t`.

## 3. Reduce + project — `reduce_dimensions.py`
```python
dyn.tl.reduceDimension(adata, reduction_method="umap")
dyn.tl.cell_velocities(adata, basis="umap", ekey="M_t", vkey="velocity_T")
```
For labeling data pass `--ekey M_t --vkey velocity_T` (or let auto-detect handle it).

## 4. Vector field + downstream
Identical to the splicing workflow: `vector_field.py`, then
`differential_geometry.py` / `fate.py` / `least_action.py` / `perturbation.py`.

## All-in-one recipes (alternative)
Dynamo also ships bundled recipes that do preprocessing + moments + reduction +
velocity in one call — handy for quick starts:
```python
dyn.tl.recipe_kin_data(adata, tkey="time", n_top_genes=1000)      # kinetics
dyn.tl.recipe_deg_data(adata, tkey="time")                        # degradation
dyn.tl.recipe_one_shot_data(adata, tkey="time")                   # one-shot
dyn.tl.recipe_mix_kin_deg_data(adata, tkey="time")                # mixed
```
The step scripts use the modular `Preprocessor` + `dynamics` path instead, so you
can inspect and iterate between stages.

## One command
```bash
python run_pipeline.py labeled.h5ad -o vf.h5ad --modality labeling --tkey time --experiment-type kin
```

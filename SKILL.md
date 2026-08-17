---
name: dynamo
description: RNA velocity and continuous vector-field analysis of single-cell RNA-seq with dynamo. Use for velocity estimation (splicing-based and metabolic-labeling scRNA-seq such as scNT-seq/scEU-seq), continuous vector-field reconstruction, differential geometry (Jacobian/acceleration/curvature/divergence to find regulatory genes), cell-fate prediction, least-action reprogramming paths, and in-silico perturbation. For standard QC/clustering/marker analysis use scanpy; for deep-learning models use scvi-tools.
license: BSD-3-Clause
metadata:
  version: "1.0"
  dynamo-version: "1.5.3"
---

> **THIS MACHINE (xing lab cluster) — local activation note (not from upstream).**
> Activate the shared venv before running any script here:
> ```
> source /net/capricorn/home/xing/lul176/mskcc/blood_combined/.venv/bin/activate
> ```
> Uses the already-present **dynamo 1.5.3** (imports as `import dynamo as dyn`;
> anndata/scanpy also present — no extra install).
> The `scripts/` here import a local `_common`, so run them from this skill's
> `scripts/` dir (e.g. `python run_pipeline.py raw.h5ad -o vf.h5ad`).
> Upstream of this skill: get an `.h5ad` with the right layers first — **spliced/
> unspliced** for splicing velocity, or **new/total** for labeling velocity (see
> the **scanpy** / **scrnaseq-fetch-process** skills for acquisition). This skill
> is the **vector-field analysis** stage.

# Dynamo: Vector-Field scRNA-seq Analysis

## Overview

Dynamo maps the **transcriptomic vector field** of single cells: it estimates
RNA velocity, reconstructs a *continuous* vector field over cell state space,
and then uses differential geometry and dynamical-systems tools on that field to
predict cell fate, optimal reprogramming paths, and the effect of genetic
perturbations. Current version here: **dynamo 1.5.3**.

Dynamo goes **beyond** velocity (which is scvelo's scope): the reconstructed
field is what enables Jacobian-based regulatory inference, least-action paths,
and in-silico perturbation. It supports both conventional splicing-based data
and **metabolic-labeling** data (scNT-seq, scEU-seq, scSLAM-seq, ...).

Everything is `anndata.AnnData`-centric (like scanpy):
`.X` expression · `.layers[]` (spliced/unspliced or new/total, plus velocity) ·
`.obs`/`.var` metadata · `.obsm[]` embeddings · `.uns['VecFld_*']` vector fields.

## When to Use This Skill

- Estimating RNA velocity from splicing (spliced/unspliced) **or** metabolic
  labeling (new/total, with a labeling time).
- Reconstructing a continuous vector field and drawing streamline/topography plots.
- Finding regulatory / driver genes via **differential geometry** (Jacobian,
  acceleration, curvature, divergence).
- Predicting **cell fate**, **least-action reprogramming paths**, or **in-silico
  perturbation / knockout** effects.

For QC, clustering, and marker genes, use **scanpy** first; for probabilistic
integration models, use **scvi-tools**.

## Script Toolkit (prefer these over writing code from scratch)

This skill bundles ready-to-run CLI scripts in `scripts/` for every stage. **Run
these instead of hand-writing dynamo code** — they handle loading by extension,
modality detection, figure setup, sensible defaults, and progress logging. Each
reads and writes `.h5ad`, so they chain together, and each has its own `--help`.
Drop down to writing dynamo code only for steps a script doesn't cover.

All scripts import the shared `scripts/_common.py` (loading, saving, figure
config, modality detection) — keep it alongside the others and run from the
`scripts/` dir. Figures default to `./figures/`.

| Script | Purpose | Typical call |
|--------|---------|--------------|
| `run_pipeline.py` | **Full core workflow in one command**: load → preprocess → dynamics → reduceDimension → cell_velocities → VectorField (+ optional diff-geometry) | `python run_pipeline.py raw.h5ad -o vf.h5ad` |
| `inspect_data.py` | Summarize a dataset: shape, **detected modality**, layers, and which dynamo steps already ran | `python inspect_data.py data.h5ad` |
| `preprocess.py` | QC + gene selection + normalize + PCA via `dyn.pp.Preprocessor` (recipe-based; handles labeling with `--tkey`) | `python preprocess.py raw.h5ad -o pp.h5ad` |
| `dynamics.py` | Kinetics + RNA velocity (`dyn.tl.dynamics`); `--tkey` enables the labeling kinetic model | `python dynamics.py pp.h5ad -o dyn.h5ad --model stochastic` |
| `reduce_dimensions.py` | Embedding (`reduceDimension`) + project velocity onto it (`cell_velocities`) | `python reduce_dimensions.py dyn.h5ad -o red.h5ad` |
| `vector_field.py` | Reconstruct the continuous field (`dyn.vf.VectorField`); `--map-topography` for fixed points | `python vector_field.py red.h5ad -o vf.h5ad --basis umap` |
| `differential_geometry.py` | Jacobian/acceleration/curvature/divergence + **ranked regulatory-gene CSVs** | `python differential_geometry.py vf.h5ad -o dg.h5ad --basis pca` |
| `fate.py` | Cell-fate prediction (`dyn.pd.fate`) from init cells + streamline/topography plots | `python fate.py vf.h5ad -o fate.h5ad --init-group cell_type:HSC` |
| `least_action.py` | Least-action path between two states (`dyn.pd.least_action`) | `python least_action.py vf.h5ad -o lap.h5ad --init-group cell_type:HSC --target-group cell_type:Meg` |
| `lap_matrix.py` | **Action matrix over ALL ordered pairs** of states + forward-vs-reverse table: which conversions are hard | `python lap_matrix.py vf_pca.h5ad -o lap.h5ad --group-col cell_type --source HSC --targets Meg Ery` |
| `perturbation.py` | In-silico perturbation / KO (`dyn.pd.perturbation`/`KO`) + ranked response | `python perturbation.py dg.h5ad -o pert.h5ad --genes GATA1 --expression 100` |
| `plot.py` | Any standard figure: streamline/topography/phase/umap/kinetic_heatmap/scalar-field | `python plot.py vf.h5ad --kind streamline --color cell_type` |

### One-shot end-to-end run

```bash
# Splicing data: raw counts (spliced/unspliced) -> reconstructed vector field + figures
python run_pipeline.py raw.h5ad -o vf.h5ad --model stochastic --basis umap
# Also compute differential geometry:
python run_pipeline.py raw.h5ad -o vf.h5ad --diff-geometry
# Metabolic labeling (scNT/scEU): pass the labeling time + experiment type
python run_pipeline.py labeled.h5ad -o vf.h5ad --modality labeling --tkey time --experiment-type kin
# Reproducible parameters via JSON (keys mirror flag names with underscores):
python run_pipeline.py raw.h5ad -o vf.h5ad --config ../assets/pipeline_config.json
```

### Step-by-step chain (inspect/iterate between stages)

```bash
python inspect_data.py         raw.h5ad
python preprocess.py           raw.h5ad -o pp.h5ad  --recipe monocle --n-top-genes 2000
python dynamics.py             pp.h5ad  -o dyn.h5ad --model stochastic
python reduce_dimensions.py    dyn.h5ad -o red.h5ad --reduction umap --basis umap
python vector_field.py         red.h5ad -o vf.h5ad  --basis umap --map-topography
# downstream (any subset):
python differential_geometry.py vf.h5ad -o dg.h5ad  --basis pca --genes GATA1 KLF1
python fate.py                  vf.h5ad -o fate.h5ad --init-group cell_type:HSC
python least_action.py          vf.h5ad -o lap.h5ad  --init-group cell_type:HSC --target-group cell_type:Meg
python lap_matrix.py            vf.h5ad -o lapm.h5ad --basis pca --group-col cell_type --source HSC
python perturbation.py          dg.h5ad -o pert.h5ad --genes GATA1 --ko
```

## Modality guide (splicing vs labeling)

Scripts auto-detect modality from layers (`--modality auto`): **spliced/unspliced**
→ `splicing`; **new/total** → `labeling`. Override with `--modality`.

- **Splicing**: no `--tkey`. Velocity uses spliced/unspliced (`M_s`/`velocity_S`).
- **Labeling** (scNT-seq, scEU-seq): set `--tkey <obs time column>` and
  `--experiment-type kin|deg|one-shot|mix_std_stm`. The kinetic model is used and
  velocity is `velocity_T` (`M_t`). See `references/labeling_workflow.md`.

The **basis** matters downstream: build the field in **umap** for streamline/fate
visualization, and in **pca** for Jacobian-based differential geometry and
perturbation (run `vector_field.py --basis pca` for those). See
`references/differential_geometry.md`.

## Key Parameters to Adjust

- **Preprocess**: `--recipe` (monocle/seurat/sctransform/pearson_residuals),
  `--n-top-genes` (2000–3000).
- **Dynamics**: `--model` (deterministic/stochastic), `--est-method` (`twostep`
  for labeling), `--tkey` (labeling time), `--group` (per-cell-type kinetics).
- **Reduce**: `--reduction` (umap/tsne), `--basis`.
- **Vector field**: `--basis`, `--map-topography`, `--pot-curl-div`.
- **Differential geometry**: `--quantities`, `--genes` (Jacobian regulators/effectors).

## Common Pitfalls and Best Practices

1. **Right layers first**: splicing needs spliced/unspliced; labeling needs
   new/total. `inspect_data.py` reports the detected modality — check it.
2. **Labeling needs `--tkey`**: without a labeling time, dynamo can't fit the
   kinetic model and falls back to steady state.
3. **PCA field for differential geometry / perturbation**: the Jacobian and
   `perturbation.py` expect the field in the `pca` basis. Build it accordingly.
4. **`cell_velocities` before `VectorField`**: the field is fit from the
   *projected* velocity, so dimension reduction + projection must run first.
5. **Fixed points / topography plots** need `vector_field.py --map-topography`.
6. **Save intermediates**: each stage writes `.h5ad`; long runs can fail partway.
7. **Reproducibility**: dynamo seeds default to 19491001; record recipe, model,
   est_method, and basis for each run.

## Bundled Resources

### scripts/ (CLI toolkit)
Composable `.h5ad`-in/`.h5ad`-out scripts plus `run_pipeline.py`. See the
**Script Toolkit** table above. Each has `--help`. `_common.py` is the shared
helper (not a CLI) — keep it alongside the others. **Default to these scripts
before writing dynamo code from scratch.**

### references/
- `standard_workflow.md` — conventional splicing end-to-end, exact `dyn.*` calls per stage.
- `labeling_workflow.md` — metabolic labeling (scNT/scEU): layer setup, recipes, `tkey`, `est_method`.
- `differential_geometry.md` — Jacobian/acceleration/curvature/divergence + ranking regulatory genes.
- `prediction.md` — fate, least-action paths, in-silico perturbation: inputs and how to read outputs.
- `api_reference.md` — quick `dyn.pp/tl/vf/pd/pl` function map.
- `plotting_guide.md` — streamline / topography / phase-portrait / kinetic-heatmap recipes.

### assets/
- `pipeline_config.json` — parameter template for `run_pipeline.py --config`.
- `analysis_template.py` — copy-and-edit full end-to-end script.

## Additional Resources

- Dynamo docs: https://dynamo-release.readthedocs.io/
- Paper: Qiu et al., *Cell* (2022) — "Mapping transcriptomic vector fields of single cells".
- Tutorials (in the repo): `docs/tutorials/` — conventional, labeling, differential geometry.

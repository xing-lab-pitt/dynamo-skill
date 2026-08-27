---
name: dynamo
description: RNA velocity and continuous vector-field analysis of single-cell RNA-seq with dynamo. Use for velocity estimation (splicing-based and metabolic-labeling scRNA-seq such as scNT-seq/scEU-seq), continuous vector-field reconstruction, differential geometry (Jacobian/acceleration/curvature/divergence to find regulatory genes), cell-fate prediction, least-action reprogramming paths, and in-silico perturbation. For standard QC/clustering/marker analysis use scanpy; for deep-learning models use scvi-tools.
license: BSD-3-Clause
metadata:
  version: "1.0"
  dynamo-version: "1.5.3"
---

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

## ⚠ Before you run: confirm the dataset's parameters with the user

**These scripts' defaults are sensible, not universal.** Velocity and
preprocessing choices are dataset-specific, and a wrong one does not error — it
produces a plausible-looking vector field pointing the **wrong way**, which then
silently corrupts every downstream figure, LAP and Jacobian.

Before running `run_pipeline.py` (or `preprocess.py` → `dynamics.py`) on
unfamiliar data, **ask the user, or check the paper / tutorial / methods section,
for**:

| Question | Why it matters |
|---|---|
| Which **genes**? An author-curated list, or select HVGs? | The paper's list is often forced (`force_gene_list`); HVG selection gives a different field |
| Splicing or labeling **velocity**, if the object has both? | Objects with spliced/unspliced *and* new/total support either; the choice can invert the field |
| `--tkey` and `--experiment-type` (`one-shot` / `kin` / `deg` / `mix_std_stm`)? | Wrong experiment type fits the wrong kinetic model |
| `--model`, `--est-method` (and `one_shot_method`, which has no flag)? | Papers pin these for reproducibility; defaults "give similar results", not identical |
| `--group` — and is it the **cell type** or the **collection timepoint**? | Kinetics are often estimated per group rather than pooled, but the two choices answer different questions and you only get one column. On an already-smoothed input this also needs `--re-smooth` — both below |
| Reuse the published embedding, or recompute? | A recomputed UMAP will not match the paper's figures |
| If `var_names` are Ensembl IDs, which **annotation release** names them? | dynamo names IDs offline from **Ensembl 77 (2014)**, so a modern curated list silently loses every gene renamed since — pin `--ensembl-release` |

**Gene-symbol vintage.** When `var_names` are Ensembl IDs, `Preprocessor` converts
them to symbols before anything can match a gene list by name. dynamo 1.5.3 does
that against Ensembl release 77 (its docstring says 109, the code says 77), giving
SEPT3, FYB, ATP5D where current HGNC says SEPTIN3, FYB1, ATP5F1D. Matching a
2022-era list against those names dropped 323 of 1956 genes on the hematopoiesis
data — no error, just a quietly smaller gene set. `preprocess.py --ensembl-release`
makes the vintage explicit; the database is downloaded and indexed once (minutes,
~1.6 GB) into `~/.cache/pyensembl`. Match the release to the list you compare
against, and record which one you used.

**`--group`: cell type or collection timepoint?** `dyn.tl.dynamics` takes a single
column, and the two natural choices answer different questions:

| `--group` is… | What you get | Reach for it when |
|---|---|---|
| a **cell type** column | Kinetic rates fitted separately per cell type — "do progenitors and neutrophils turn RNA over at different speeds?" | Rates are the result you are after, and all cells were collected together |
| a **collection timepoint** column | Each timepoint's cells smoothed only against their own timepoint, so collection times cannot bleed into one another | Cells were harvested at several times and you do not want that structure averaged away before the analysis starts |

You cannot have both from one column. Building a joint `timepoint_celltype` column
works in principle but multiplies the group count, and dynamo warns below 50 cells
per group because such groups can return all-NaN velocities, asking you to
coarse-grain instead — so the joint column can push groups under that floor. `dynamics.py` prints the
per-group counts before it runs, which is the number to look at before committing
to a grouping.

**`--group` on an already-smoothed object needs `--re-smooth`.** Smoothing replaces
each cell's value with a kNN average, and neighbours are picked by expression
similarity alone. `--group` instead builds one neighbour graph per group, which is
how cells from different collection timepoints are kept out of each other's
averages. But `dyn.tl.dynamics` rebuilds those `M_*` layers only when fewer than
two exist or `re_smooth=True`, so on an input that already carries them — a
published object (`sample_data.hematopoiesis` ships five), anything a
`recipe_*_data` produced, or a second `dynamics` run — `--group` is accepted and
then ignored, leaving a warning in the log and pooled averages under every later
number. `dynamics.py` and `run_pipeline.py` now refuse that combination and name
the fix, and print the per-group cell counts, flagging groups under 50 cells where
dynamo can return NaN velocities:

```bash
python dynamics.py hsc.h5ad -o dyn.h5ad --group time --re-smooth
```

A fresh `preprocess.py` → `dynamics.py` chain is unaffected: preprocessing writes
no `M_*` layers, so `dynamics` builds them respecting `--group` on the first run.

**Worked example — dynamo's own hematopoiesis data.** Tutorial 301 pins
`force_gene_list=adata.uns["genes_to_use"]`, sets
`adata.uns["pp"]["has_splicing"] = False` to force the labeling-only model, and
uses `group="time", one_shot_method="sci_fate", model="deterministic"`. Running
this dataset on defaults instead — HVG selection, splicing velocity — yields a
field flowing from mature monocytes *back into* progenitors.

Some of these (`force_gene_list`, `has_splicing`, `one_shot_method`) are
deliberately **not** CLI flags, as they are per-dataset rather than general. Set
them around the scripts when a source requires them:

```python
import dynamo as dyn
adata = dyn.read_h5ad("raw.h5ad")
pp = dyn.pp.Preprocessor(force_gene_list=adata.uns["genes_to_use"])
pp.config_monocle_recipe(adata, n_top_genes=len(adata.uns["genes_to_use"]))
pp.preprocess_adata_monocle(adata, tkey="time", experiment_type="one-shot")
adata.uns["pp"]["has_splicing"] = False          # force the labeling model
dyn.tl.dynamics(adata, group="time", one_shot_method="sci_fate", model="deterministic")
adata.write("pp.h5ad")                            # then resume with the scripts
```

**Reproducing a published figure?** Prefer the processed object
(`sample_data.py --dataset hematopoiesis`) over re-deriving it. It ships the
paper's embedding, velocity, vector field and Jacobians, so the field cannot
disagree with the publication. Re-derive only when the derivation itself is the
question — and then verify the flow direction against known biology before
trusting anything downstream.

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
| `sample_data.py` | Fetch a **bundled published dataset** (`dyn.sample_data`) — most tutorials start from a processed object rather than re-deriving one | `python sample_data.py --list` · `python sample_data.py --dataset hematopoiesis -o hsc.h5ad` |
| `run_pipeline.py` | **Full core workflow in one command**: load → preprocess → dynamics → reduceDimension → cell_velocities → VectorField (+ optional diff-geometry) | `python run_pipeline.py raw.h5ad -o vf.h5ad` |
| `inspect_data.py` | Summarize a dataset: shape, **detected modality**, layers, and which dynamo steps already ran | `python inspect_data.py data.h5ad` |
| `preprocess.py` | QC + gene selection + normalize + PCA via `dyn.pp.Preprocessor` (recipe-based; handles labeling with `--tkey`) | `python preprocess.py raw.h5ad -o pp.h5ad` |
| `dynamics.py` | Kinetics + RNA velocity (`dyn.tl.dynamics`); `--tkey` enables the labeling kinetic model; `--group` estimates per group and needs `--re-smooth` on an already-smoothed input | `python dynamics.py pp.h5ad -o dyn.h5ad --model stochastic` |
| `reduce_dimensions.py` | Embedding (`reduceDimension`) + project velocity onto it (`cell_velocities`) | `python reduce_dimensions.py dyn.h5ad -o red.h5ad` |
| `vector_field.py` | Reconstruct the continuous field (`dyn.vf.VectorField`); `--map-topography` for fixed points | `python vector_field.py red.h5ad -o vf.h5ad --basis umap` |
| `topography.py` | **Inspect and curate fixed points**: `--list` what was detected, `--keep` the real ones. Raw topography over-detects, so every published figure is curated | `python topography.py vf.h5ad --basis umap --n 750 --list` |
| `differential_geometry.py` | Jacobian/acceleration/curvature/divergence/**graph potential** + **ranked regulatory-gene CSVs** | `python differential_geometry.py vf.h5ad -o dg.h5ad --basis pca` |
| `fate.py` | Cell-fate prediction (`dyn.pd.fate`) from init cells + streamline/topography plots | `python fate.py vf.h5ad -o fate.h5ad --init-group cell_type:HSC` |
| `least_action.py` | Least-action path between two states (`dyn.pd.least_action`) | `python least_action.py vf.h5ad -o lap.h5ad --init-group cell_type:HSC --target-group cell_type:Meg` |
| `lap_matrix.py` | **Action matrix over ALL ordered pairs** of states + forward-vs-reverse table: which conversions are hard. `--save-paths` also keeps the path geometry for plotting | `python lap_matrix.py vf_pca.h5ad -o lap.h5ad --group-col cell_type --source HSC --targets Meg Ery` |
| `perturbation.py` | In-silico perturbation / KO (`dyn.pd.perturbation`/`KO`) + ranked response | `python perturbation.py dg.h5ad -o pert.h5ad --genes GATA1 --expression 100` |
| `plot.py` | Any standard figure: streamline/topography/phase/umap/kinetic_heatmap/scalar-field, plus `scalar_by_group`, `jacobian_heatmap`, `lap_paths` | `python plot.py vf.h5ad --kind streamline --color cell_type` |

**`plot.py --kind` cheat sheet for the less obvious ones:**

```bash
# compare a per-cell scalar ACROSS groups (painting it on the embedding cannot show this)
python plot.py dg.h5ad --kind scalar_by_group --scalar speed_pca --group cell_type
# a named regulator -> effector Jacobian, and the full pairwise matrix
python plot.py dg.h5ad --kind jacobian         --genes SPI1 GATA1
python plot.py dg.h5ad --kind jacobian_heatmap --genes SPI1 GATA1 KLF1 CEBPA
# overlay least-action paths saved by lap_matrix.py --save-paths
python plot.py vf.h5ad --kind lap_paths --paths results/lap/laps.pkl \
    --pairs "HSC->Meg" "HSC->Ery" --out-name lap_development
```

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
  for labeling), `--tkey` (labeling time), `--group` (per-group kinetics: cell
  type or collection timepoint), `--re-smooth` (rebuild the `M_*` layers; required
  alongside `--group` on an already-smoothed input).
- **Reduce**: `--reduction` (umap/tsne), `--basis`.
- **Vector field**: `--basis`, `--map-topography`, `--pot-curl-div`.
- **Differential geometry**: `--quantities`, `--genes` (Jacobian regulators/effectors).

## Common Pitfalls and Best Practices

0. **Confirm dataset-specific parameters before running** — see the warning
   section above. A wrong gene list or modality choice inverts the field silently.
   **Always sanity-check the streamline direction against known biology before
   reading anything downstream**; progenitors should flow toward mature states.
1. **Right layers first**: splicing needs spliced/unspliced; labeling needs
   new/total. `inspect_data.py` reports the detected modality — check it. An
   object carrying **both** supports either velocity, and the two can disagree;
   the source, not the default, decides which is right.
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

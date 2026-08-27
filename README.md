# dynamo-skill

A [Claude Code](https://claude.com/claude-code) **skill** — a toolkit of
composable command-line scripts for RNA velocity and continuous **vector-field**
analysis of single-cell RNA-seq with [dynamo](https://github.com/aristoteleo/dynamo-release)
(v1.5.3). It works equally as a plain script toolkit you run yourself.

Dynamo goes beyond velocity estimation: it reconstructs a *continuous* vector
field over cell-state space and uses differential geometry and dynamical-systems
tools on that field to find regulatory genes, predict cell fate, compute
least-action reprogramming paths, and simulate in-silico perturbations. This
skill covers both **splicing** (spliced/unspliced) and **metabolic-labeling**
(new/total; scNT-seq, scEU-seq) data.

## Pipeline

```
load → preprocess → dynamics (velocity) → reduceDimension → cell_velocities
     → VectorField → { differential geometry | fate | least-action path | perturbation } → plots
```

Every stage is a script that reads and writes `.h5ad`, so they chain together;
`run_pipeline.py` runs the core chain in one command.

## Requirements

- Python with [`dynamo-release`](https://pypi.org/project/dynamo-release/) 1.5.3
  (`pip install dynamo-release`), plus `anndata`.
- Input: an `.h5ad` with the right layers — **spliced/unspliced** for splicing
  velocity, or **new/total** (plus a labeling-time column) for labeling velocity.

## Quickstart

```bash
cd scripts

# One-shot: raw counts -> reconstructed vector field (+ figures)
python run_pipeline.py raw.h5ad -o vf.h5ad --model stochastic --basis umap

# Metabolic labeling (scNT/scEU): pass the labeling time + experiment type
python run_pipeline.py labeled.h5ad -o vf.h5ad --modality labeling --tkey time --experiment-type kin
```

Or run stage by stage, inspecting between steps:

```bash
cd scripts
python inspect_data.py          raw.h5ad
python preprocess.py            raw.h5ad -o pp.h5ad  --recipe monocle --n-top-genes 2000
python dynamics.py              pp.h5ad  -o dyn.h5ad --model stochastic
python reduce_dimensions.py     dyn.h5ad -o red.h5ad --reduction umap --basis umap
python vector_field.py          red.h5ad -o vf.h5ad  --basis umap --map-topography
# downstream (any subset):
python differential_geometry.py vf.h5ad  -o dg.h5ad  --basis pca --genes GATA1 KLF1
python fate.py                  vf.h5ad  -o fate.h5ad --init-group cell_type:HSC
python least_action.py          vf.h5ad  -o lap.h5ad  --init-group cell_type:HSC --target-group cell_type:Meg
python perturbation.py          dg.h5ad  -o pert.h5ad --genes GATA1 --ko
```

Every script has `--help`. Scripts import the shared `scripts/_common.py`, so run
them from the `scripts/` directory; figures default to `./figures/`.

## A note on parameters

Velocity and preprocessing choices are dataset-specific, and a wrong one usually
does not raise — it yields a plausible field pointing the wrong way, which then
corrupts every figure, least-action path and Jacobian downstream. `SKILL.md`
carries a checklist of what to confirm before a run. One case is enforced in code:
estimating kinetics per group (`--group`, e.g. per cell type or collection
timepoint) also needs `--re-smooth` when the input was already smoothed, because
dynamo would otherwise accept `--group` and quietly keep the pooled averages. The
scripts refuse that combination rather than proceeding.

## Layout

```
SKILL.md            Full skill documentation (read this for details)
scripts/            CLI toolkit (.h5ad in -> .h5ad out); _common.py is the shared helper
references/         standard & labeling workflows, differential geometry, prediction, API map, plotting
assets/             pipeline_config.json (run_pipeline --config template) + analysis_template.py
```

## Documentation

- **[`SKILL.md`](SKILL.md)** — full reference: script table, modality guide, key
  parameters, and common pitfalls.
- `references/` — per-topic guides with the exact `dyn.*` calls.
- Dynamo docs: https://dynamo-release.readthedocs.io/
- Paper: Qiu et al., *Cell* (2022), "Mapping transcriptomic vector fields of single cells".

## License

BSD 3-Clause — see [`LICENSE`](LICENSE).

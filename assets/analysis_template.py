#!/usr/bin/env python3
"""
Copy-and-edit dynamo analysis template — the full vector-field workflow as plain
`dyn.*` calls, for when you need to go beyond what the toolkit scripts expose.

Run inside the venv:
    source /net/capricorn/home/xing/lul176/mskcc/blood_combined/.venv/bin/activate
    python analysis_template.py

Splicing (spliced/unspliced) is shown; the labeling switches are marked LABELING.
"""

import dynamo as dyn

# Headless figures.
dyn.configuration.set_figure_params("dynamo", background="white", dpi=100)

INPUT = "raw.h5ad"          # must have spliced/unspliced (splicing) or new/total (labeling)
OUTPUT = "vf.h5ad"
FIGDIR = "figures"

# LABELING: set these for scNT/scEU data (else leave TKEY=None).
TKEY = None                 # e.g. "time"
EXPERIMENT_TYPE = None      # "kin" | "deg" | "one-shot" | "mix_std_stm"

# ---------------------------------------------------------------- load
adata = dyn.read_h5ad(INPUT)

# ---------------------------------------------------------------- preprocess
pp = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
pp.config_monocle_recipe(adata, n_top_genes=2000)
pp.preprocess_adata(adata, recipe="monocle", tkey=TKEY, experiment_type=EXPERIMENT_TYPE)

# ---------------------------------------------------------------- dynamics / velocity
# Splicing: model="stochastic". LABELING: pass tkey + est_method="twostep".
if TKEY:
    dyn.tl.dynamics(adata, model="deterministic", est_method="twostep", tkey=TKEY, cores=4)
    EKEY, VKEY = "M_t", "velocity_T"
else:
    dyn.tl.dynamics(adata, model="stochastic", cores=4)
    EKEY, VKEY = None, None  # auto -> M_s / velocity_S

# ---------------------------------------------------------------- embed + project
dyn.tl.reduceDimension(adata, reduction_method="umap")
dyn.tl.cell_velocities(adata, basis="umap", ekey=EKEY, vkey=VKEY)

# ---------------------------------------------------------------- vector field (umap for viz)
dyn.vf.VectorField(adata, basis="umap", map_topography=True)
dyn.pl.streamline_plot(adata, basis="umap", color="Cell_type", save_show_or_return="return")

# ---------------------------------------------------------------- differential geometry (pca)
dyn.tl.cell_velocities(adata, basis="pca")          # project velocity onto pca
dyn.vf.VectorField(adata, basis="pca")              # field in state space
GENES = ["GATA1", "KLF1", "SPI1"]                   # regulators/effectors of interest
dyn.vf.jacobian(adata, regulators=GENES, effectors=GENES, basis="pca")
dyn.vf.acceleration(adata, basis="pca")
dyn.vf.curvature(adata, basis="pca")
jac_ranked = dyn.vf.rank_jacobian_genes(adata, return_df=True)   # {group: DataFrame}

# ---------------------------------------------------------------- prediction
# Fate (umap):
hsc = list(adata.obs_names[adata.obs["Cell_type"] == "HSC"])[:1]
dyn.pd.fate(adata, init_cells=hsc, basis="umap", direction="forward")

# Least-action path (pca) — one representative cell per endpoint:
meg = list(adata.obs_names[adata.obs["Cell_type"] == "Meg"])[:1]
dyn.pd.least_action(adata, init_cells=hsc, target_cells=meg, basis="pca", n_points=25)

# In-silico perturbation / KO (pca -> umap):
dyn.pd.perturbation(adata, genes="GATA1", expression=100, basis="pca", emb_basis="umap")
dyn.pl.streamline_plot(adata, basis="umap_perturbation", save_show_or_return="return")
# dyn.pd.KO(adata, KO_genes="GATA1", basis="pca", emb_basis="umap")  # -> umap_KO

# ---------------------------------------------------------------- save
adata.write_h5ad(OUTPUT)
print("done ->", OUTPUT)

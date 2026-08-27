#!/usr/bin/env python3
"""
End-to-end dynamo vector-field pipeline in one command.

Runs the core workflow from raw counts to a reconstructed vector field:
    load -> preprocess (Preprocessor) -> dynamics (velocity)
    -> reduceDimension -> cell_velocities -> VectorField
    -> (optional) differential geometry.

Handles both modalities: pass --modality (or let it auto-detect from layers) and,
for labeling data, --tkey / --experiment-type. Tune knobs via flags or a JSON
--config (keys mirror the flag names with underscores). Use the individual step
scripts when you need to inspect or iterate on a single stage.

Examples:
    python run_pipeline.py raw.h5ad -o vf.h5ad
    python run_pipeline.py raw.h5ad -o vf.h5ad --model stochastic --basis umap
    python run_pipeline.py labeled.h5ad -o vf.h5ad --modality labeling --tkey time --experiment-type kin
    python run_pipeline.py raw.h5ad -o vf.h5ad --config params.json --diff-geometry
"""

import argparse
import json

from _common import (check_group_smoothing, configure_dynamo, info, load_adata,
                     resolve_modality, save_adata)


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Input raw-counts .h5ad")
    p.add_argument("-o", "--output", default="vf.h5ad", help="Output .h5ad")
    p.add_argument("--figdir", default="figures", help="Figure directory")
    p.add_argument("--config", default=None, help="JSON overriding any option below")
    # modality / preprocessing
    p.add_argument("--modality", default="auto", choices=["auto", "splicing", "labeling"])
    p.add_argument("--recipe", default="monocle",
                   choices=["monocle", "seurat", "sctransform", "pearson_residuals",
                            "monocle_pearson_residuals"])
    p.add_argument("--n-top-genes", type=int, default=2000)
    p.add_argument("--tkey", default=None)
    p.add_argument("--experiment-type", default=None,
                   choices=["conventional", "kin", "deg", "one-shot", "mix_std_stm"])
    p.add_argument("--no-cell-cycle", action="store_true")
    # dynamics
    p.add_argument("--model", default="auto", choices=["auto", "deterministic", "stochastic"])
    p.add_argument("--est-method", default="auto",
                   choices=["auto", "ols", "rlm", "ransac", "gmm", "negbin", "twostep", "direct"])
    p.add_argument("--group", default=None,
                   help="obs column for per-group kinetics: a cell type, or the "
                        "collection timepoint to keep timepoints from mixing")
    p.add_argument("--re-smooth", action="store_true",
                   help="Rebuild the smoothed M_* layers; needed with --group on an "
                        "input that was already smoothed")
    p.add_argument("--cores", type=int, default=1)
    # reduction / projection / field
    p.add_argument("--reduction", default="umap", choices=["umap", "tsne", "psl", "sude"])
    p.add_argument("--basis", default="umap", help="Basis for embedding, projection and field")
    p.add_argument("--map-topography", action="store_true")
    # optional downstream
    p.add_argument("--diff-geometry", action="store_true",
                   help="Also compute Jacobian/acceleration/curvature/divergence in --basis")
    return p


def apply_config(args):
    if args.config:
        with open(args.config) as fh:
            cfg = json.load(fh)
        for key, value in cfg.items():
            setattr(args, key.replace("-", "_"), value)
    return args


def main():
    args = apply_config(build_parser().parse_args())
    dyn = configure_dynamo(figdir=args.figdir)

    info("[1/6] Loading data")
    adata = load_adata(args.input)
    modality = resolve_modality(adata, args.modality)
    info(f"      {adata.n_obs} cells x {adata.n_vars} genes  (modality={modality})")

    info(f"[2/6] Preprocess (recipe={args.recipe})")
    pp = dyn.pp.Preprocessor(cell_cycle_score_enable=not args.no_cell_cycle)
    if args.recipe == "monocle":
        pp.config_monocle_recipe(adata, n_top_genes=args.n_top_genes)
    else:
        pp.select_genes_kwargs["n_top_genes"] = args.n_top_genes
    pp.preprocess_adata(adata, recipe=args.recipe, tkey=args.tkey,
                        experiment_type=args.experiment_type)

    # Labeling kinetics: the 'auto' estimator can hit numerical bounds; twostep is
    # dynamo's robust default for scEU/scNT data.
    est_method = args.est_method
    if args.tkey and est_method == "auto":
        est_method = "twostep"
        info("      labeling kinetics: est_method auto -> twostep")
    check_group_smoothing(adata, args.group, args.re_smooth)
    info(f"[3/6] Dynamics (model={args.model}, est_method={est_method})")
    dyn.tl.dynamics(adata, model=args.model, est_method=est_method,
                    tkey=args.tkey, group=args.group,
                    re_smooth=args.re_smooth, cores=args.cores)

    info(f"[4/6] Dimension reduction ({args.reduction}) + cell velocities (basis={args.basis})")
    dyn.tl.reduceDimension(adata, reduction_method=args.reduction, cores=args.cores)
    dyn.tl.cell_velocities(adata, basis=args.basis)

    info(f"[5/6] Vector field (basis={args.basis})")
    dyn.vf.VectorField(adata, basis=args.basis, map_topography=args.map_topography,
                       cores=args.cores)

    if args.diff_geometry:
        info("[6/6] Differential geometry (jacobian, acceleration, curvature, divergence)")
        dyn.vf.jacobian(adata, basis=args.basis)
        for q in ("acceleration", "curvature", "divergence"):
            getattr(dyn.vf, q)(adata, basis=args.basis)
    else:
        info("[6/6] Skipping differential geometry (pass --diff-geometry to enable)")

    save_adata(adata, args.output)
    info("Pipeline complete. Next: fate.py / least_action.py / perturbation.py / plot.py")


if __name__ == "__main__":
    main()

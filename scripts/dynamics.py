#!/usr/bin/env python3
"""
Estimate RNA kinetics and velocity with dyn.tl.dynamics.

This is the core velocity step. It fits transcription / splicing / degradation
(and, for labeling data, RNA labeling) rates per gene and writes velocity into
the layers. Run after preprocess.py.

  * splicing data -> steady-state / stochastic model on spliced/unspliced.
  * labeling data -> pass --tkey (labeling time) so the kinetic model is used;
                     --est-method twostep is the usual choice for scEU/scNT data.

Examples:
    python dynamics.py pp.h5ad -o dyn.h5ad
    python dynamics.py pp.h5ad -o dyn.h5ad --model stochastic
    python dynamics.py pp.h5ad -o dyn.h5ad --tkey time --est-method twostep --cores 8
"""

import argparse

from _common import (add_io_args, check_group_smoothing, configure_dynamo, info,
                     load_adata, save_adata)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="dyn.h5ad")
    p.add_argument("--model", default="auto",
                   choices=["auto", "deterministic", "stochastic"],
                   help="Kinetics model (default: auto)")
    p.add_argument("--est-method", default="auto",
                   choices=["auto", "ols", "rlm", "ransac", "gmm", "negbin", "twostep", "direct"],
                   help="Estimation method (default: auto). twostep suits labeling data.")
    p.add_argument("--tkey", default=None,
                   help="obs column with labeling time (enables the kinetic model)")
    p.add_argument("--group", default=None,
                   help="obs column to estimate kinetics per group. A cell-type column "
                        "gives cell-type-specific rates; a collection-timepoint column "
                        "instead keeps timepoints out of each other's kNN averages. "
                        "Only one column, so pick the question you are asking.")
    p.add_argument("--re-smooth", action="store_true",
                   help="Rebuild the smoothed M_* layers. Required with --group when "
                        "the input was already smoothed, or dynamo ignores --group.")
    p.add_argument("--cores", type=int, default=1, help="Parallel cores (default: 1)")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)
    check_group_smoothing(adata, args.group, args.re_smooth)

    # For labeling kinetics (--tkey) the 'auto' estimator can hit numerical bounds;
    # twostep is dynamo's robust default for scEU/scNT data.
    est_method = args.est_method
    if args.tkey and est_method == "auto":
        est_method = "twostep"
        info("[dynamics] labeling kinetics: est_method auto -> twostep")

    info(f"[dynamics] model={args.model}  est_method={est_method}  "
         f"tkey={args.tkey}  group={args.group}  re_smooth={args.re_smooth}")
    dyn.tl.dynamics(
        adata,
        model=args.model,
        est_method=est_method,
        tkey=args.tkey,
        group=args.group,
        re_smooth=args.re_smooth,
        cores=args.cores,
    )

    save_adata(adata, args.output)
    info("Next: reduce_dimensions.py (embedding + project velocity).")


if __name__ == "__main__":
    main()

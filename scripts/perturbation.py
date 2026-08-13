#!/usr/bin/env python3
"""
In-silico genetic perturbation: predict how up/down-regulating (or knocking out)
genes reshapes the vector field and cell-fate flow (dyn.pd.perturbation / KO),
then rank the most-affected genes and cells.

Perturbation propagates through the Jacobian, so compute it first in the PCA
basis:  vector_field.py --basis pca  ->  differential_geometry.py --quantities jacobian --basis pca.

  --expression > 0 activates, < 0 suppresses a gene. --ko forces a knockout.

Examples:
    python perturbation.py dg.h5ad -o pert.h5ad --genes GATA1 --expression 100
    python perturbation.py dg.h5ad -o pert.h5ad --genes SPI1 --ko
"""

import argparse
import os

from _common import add_io_args, configure_dynamo, info, load_adata, save_adata, save_fig


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="pert.h5ad")
    p.add_argument("--genes", nargs="+", required=True, help="Gene(s) to perturb")
    p.add_argument("--expression", type=float, default=10,
                   help="Target expression (>0 activate, <0 suppress; default: 10)")
    p.add_argument("--ko", action="store_true", help="Knockout instead of continuous perturbation")
    p.add_argument("--basis", default="pca", help="Field basis for propagation (default: pca)")
    p.add_argument("--emb-basis", default="umap", help="Embedding to project effect onto (default: umap)")
    p.add_argument("--top-n", type=int, default=25, help="Top genes/cells per rank table")
    p.add_argument("--out-dir", default="results/perturbation", help="CSV output directory")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)

    if args.ko:
        # KO re-simulates the field with the gene(s) zeroed, writing X/velocity_<emb>_KO.
        # It does not produce the delta-X perturbation matrix that rank_perturbation_* consume.
        info(f"[perturbation] KO {args.genes}")
        dyn.pd.KO(adata, KO_genes=args.genes, basis=args.basis, emb_basis=args.emb_basis)
        plot_basis = f"{args.emb_basis}_KO"
    else:
        info(f"[perturbation] {args.genes} -> expression={args.expression}")
        dyn.pd.perturbation(adata, genes=args.genes, expression=args.expression,
                            basis=args.basis, emb_basis=args.emb_basis)
        plot_basis = f"{args.emb_basis}_perturbation"

        # Rank the perturbation response (only meaningful for continuous perturbation).
        os.makedirs(args.out_dir, exist_ok=True)
        rank_tables = {
            "genes": ("rank_perturbation_genes", "rank_j_delta_x_perturbation"),
            "cells": ("rank_perturbation_cells", "rank_j_delta_x_perturbation_cells"),
        }
        for kind, (fn_name, uns_key) in rank_tables.items():
            fn = getattr(dyn.pd, fn_name, None)
            if fn is None:
                continue
            try:
                fn(adata)  # writes adata.uns[uns_key]
                df = adata.uns.get(uns_key)
                if df is not None and hasattr(df, "to_csv"):
                    out = os.path.join(args.out_dir, f"perturbation_rank_{kind}.csv")
                    df.head(args.top_n).to_csv(out)
                    info(f"      wrote {out}")
            except Exception as e:
                info(f"      ranking {kind} skipped ({e})")

    # Streamline of the perturbed / KO'd field on the embedding.
    try:
        dyn.pl.streamline_plot(adata, basis=plot_basis, save_show_or_return="return")
        save_fig(f"{args.figdir}/perturbation_streamline.png")
    except Exception as e:
        info(f"perturbation streamline skipped ({e})")

    save_adata(adata, args.output)


if __name__ == "__main__":
    main()

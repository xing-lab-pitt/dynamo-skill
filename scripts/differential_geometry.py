#!/usr/bin/env python3
"""
Differential geometry of the reconstructed vector field: compute Jacobian,
acceleration, curvature, divergence, curl, and speed, then rank genes by each
to surface regulatory / driver genes. This is dynamo's signature analysis for
extracting regulatory function from the vector field.

Requires a vector field (run vector_field.py first). The Jacobian and most
ranking need the field in the PCA basis, so build the field with
`vector_field.py --basis pca` when you plan to run this.

Ranked-gene tables are written as CSVs to --out-dir.

Examples:
    python differential_geometry.py vf.h5ad -o dg.h5ad --basis pca
    python differential_geometry.py vf.h5ad -o dg.h5ad \\
        --quantities jacobian acceleration curvature --genes GATA1 KLF1 SPI1
"""

import argparse
import os

from _common import add_io_args, configure_dynamo, die, info, load_adata, save_adata


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="dg.h5ad")
    p.add_argument("--basis", default="pca",
                   help="Vector-field basis to differentiate (default: pca)")
    p.add_argument("--quantities", nargs="+",
                   default=["jacobian", "acceleration", "curvature", "divergence"],
                   choices=["jacobian", "acceleration", "curvature", "divergence",
                            "curl", "speed"],
                   help="Geometric quantities to compute")
    p.add_argument("--genes", nargs="+", default=None,
                   help="Regulators/effectors for the Jacobian (default: all velocity genes)")
    p.add_argument("--no-rank", action="store_true",
                   help="Skip ranking genes after computing quantities")
    p.add_argument("--top-n", type=int, default=25, help="Top genes per rank table")
    p.add_argument("--out-dir", default="results/diffgeo", help="CSV output directory")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)

    # Compute each requested quantity.
    for q in args.quantities:
        info(f"[diffgeo] computing {q} (basis={args.basis})")
        fn = getattr(dyn.vf, q)
        if q == "jacobian":
            fn(adata, regulators=args.genes, effectors=args.genes, basis=args.basis)
        else:
            fn(adata, basis=args.basis)

    if not args.no_rank:
        os.makedirs(args.out_dir, exist_ok=True)
        # Per-gene rankable scalars: rank_*_genes stores a DataFrame in uns['rank_<q>'].
        # (divergence/curl/speed are per-cell scalars — computed for plotting/coloring,
        # not gene-ranked — so they are intentionally absent here.)
        scalar_rankers = {
            "acceleration": "rank_acceleration_genes",
            "curvature": "rank_curvature_genes",
        }
        for q in args.quantities:
            if q not in scalar_rankers:
                continue
            rank_fn = getattr(dyn.vf, scalar_rankers[q], None)
            if rank_fn is None:
                continue
            try:
                info(f"[diffgeo] ranking genes by {q}")
                rank_fn(adata)  # writes adata.uns['rank_<q>']
                df = adata.uns.get(f"rank_{q}")
                if df is not None and hasattr(df, "to_csv"):
                    out = os.path.join(args.out_dir, f"rank_{q}.csv")
                    df.head(args.top_n).to_csv(out, index=False)
                    info(f"      wrote {out}")
            except Exception as e:  # ranking is best-effort; don't abort the run
                info(f"      ranking {q} skipped ({e})")

        # Jacobian ranking needs the gene-wise Jacobian, which only exists when
        # jacobian() was given explicit regulators/effectors (--genes).
        if "jacobian" in args.quantities:
            if not args.genes:
                info("      ranking jacobian skipped (needs --genes: regulators/effectors)")
            else:
                try:
                    info("[diffgeo] ranking genes by jacobian")
                    # return_df=True -> dict {group: DataFrame}
                    tables = dyn.vf.rank_jacobian_genes(adata, return_df=True)
                    for group, df in (tables or {}).items():
                        if hasattr(df, "to_csv"):
                            out = os.path.join(args.out_dir, f"rank_jacobian_{group}.csv")
                            df.head(args.top_n).to_csv(out, index=False)
                            info(f"      wrote {out}")
                except Exception as e:
                    info(f"      ranking jacobian skipped ({e})")

    save_adata(adata, args.output)
    info("Tip: plot.py --kind jacobian / acceleration to visualize these on the embedding.")


if __name__ == "__main__":
    main()

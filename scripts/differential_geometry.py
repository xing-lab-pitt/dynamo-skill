#!/usr/bin/env python3
"""
Differential geometry of the reconstructed vector field: compute Jacobian,
acceleration, curvature, divergence, curl, speed, and the graph (ddhodge)
potential, then rank genes by each to surface regulatory / driver genes. This
is dynamo's signature analysis for extracting regulatory function from the
vector field.

Requires a vector field (run vector_field.py first). The Jacobian and most
ranking need the field in the PCA basis, so build the field with
`vector_field.py --basis pca` when you plan to run this.

`potential` is different from the others: it is a graph quantity, derived by
Hodge decomposition of a cell-cell transition matrix in obsp (--transition-key)
rather than by differentiating the field in a basis. It yields the vector-field
pseudotime — obs['potential_<key>'] and obs['pseudotime_<key>'].

Ranked-gene tables are written as CSVs to --out-dir.

Examples:
    python differential_geometry.py vf.h5ad -o dg.h5ad --basis pca
    python differential_geometry.py vf.h5ad -o dg.h5ad \\
        --quantities jacobian acceleration curvature --genes GATA1 KLF1 SPI1
    python differential_geometry.py vf.h5ad -o dg.h5ad \\
        --quantities potential --transition-key fp_transition_rate
"""

import argparse
import os

from _common import add_io_args, configure_dynamo, die, info, load_adata, save_adata


def compute_potential(adata, key, sign="auto"):
    """Hodge-decomposition potential / pseudotime from an obsp transition matrix.

    Writes obs['potential_<key>'] and obs['pseudotime_<key>'] (= -potential).

    The sign of the divergence fed to potential() depends on the matrix's
    orientation, and dynamo's own tutorials differ: a cosine transition matrix
    is used as potential(g, -div(g)), while the Fokker-Planck rate matrix
    fp_transition_rate is used as potential(g, +div(g)). "auto" picks by key and
    logs the choice; --potential-sign overrides.

    The two conventions orient the potential OPPOSITELY, so pseudotime cannot be
    a blanket -potential:
      * neg (cosine-like): potential is LOW at the source, matching dynamo's own
        obs['*_ddhodge_potential'] -> pseudotime = +potential
      * pos (fp-like):     potential is HIGH at the source -> pseudotime = -potential
    Either way pseudotime increases away from the progenitor. Getting this
    backwards inverts the trajectory silently, which reads as a biological
    result rather than a bug — so the orientation is asserted here, not left to
    the caller.
    """
    from dynamo.tools.graph_operators import build_graph, div, potential

    if key not in adata.obsp:
        info(f"      potential skipped: obsp['{key}'] not found; have {list(adata.obsp.keys())}")
        return False
    if sign == "auto":
        sign = "pos" if "fp_transition" in key else "neg"
    info(f"[diffgeo] potential from obsp['{key}'] (div sign: {sign})")
    g = build_graph(adata.obsp[key])
    d = div(g)
    pot = potential(g, d if sign == "pos" else -d)
    adata.obs[f"potential_{key}"] = pot
    adata.obs[f"pseudotime_{key}"] = -pot if sign == "pos" else pot
    info(f"      wrote obs['potential_{key}'] and obs['pseudotime_{key}'] "
         f"(pseudotime = {'-' if sign == 'pos' else '+'}potential)")
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="dg.h5ad")
    p.add_argument("--basis", default="pca",
                   help="Vector-field basis to differentiate (default: pca)")
    p.add_argument("--quantities", nargs="+",
                   default=["jacobian", "acceleration", "curvature", "divergence"],
                   choices=["jacobian", "acceleration", "curvature", "divergence",
                            "curl", "speed", "potential"],
                   help="Geometric quantities to compute")
    p.add_argument("--transition-key", default="cosine_transition_matrix",
                   help="obsp transition matrix for --quantities potential "
                        "(default: cosine_transition_matrix; fp_transition_rate is "
                        "the other common choice)")
    p.add_argument("--potential-sign", choices=["auto", "neg", "pos"], default="auto",
                   help="Sign of div() fed to potential(). auto: pos for "
                        "fp_transition* keys, neg otherwise (default)")
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
        # potential is a graph quantity, not a derivative of the field in a basis
        if q == "potential":
            compute_potential(adata, args.transition_key, args.potential_sign)
            continue
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

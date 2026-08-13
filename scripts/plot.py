#!/usr/bin/env python3
"""
Generate a standard dynamo figure from a processed object. One --kind per call;
the figure is saved to --figdir with a predictable name.

Kinds:
  streamline   velocity streamlines on the embedding (needs cell_velocities)
  topography   streamlines + fixed points (needs vector_field.py --map-topography)
  phase        per-gene phase portraits (spliced vs unspliced / new vs total)
  umap         plain embedding colored by --color
  acceleration | curvature | divergence | speed | jacobian
               differential-geometry scalar fields (needs differential_geometry.py)
  kinetic_heatmap  gene expression kinetics ordered along the field

Examples:
    python plot.py vf.h5ad --kind streamline --color cell_type
    python plot.py dyn.h5ad --kind phase --genes GATA1 KLF1
    python plot.py dg.h5ad --kind acceleration --basis umap
"""

import argparse

from _common import configure_dynamo, die, info, load_adata, save_fig

SCALAR_FIELDS = {"acceleration", "curvature", "divergence", "speed", "jacobian"}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Processed .h5ad")
    p.add_argument("--kind", required=True,
                   choices=["streamline", "topography", "phase", "umap",
                            "kinetic_heatmap"] + sorted(SCALAR_FIELDS))
    p.add_argument("--basis", default="umap", help="Embedding basis (default: umap)")
    p.add_argument("--color", nargs="+", default=None, help="obs/var key(s) to color by")
    p.add_argument("--genes", nargs="+", default=None, help="Genes (phase / kinetic_heatmap)")
    p.add_argument("--figdir", default="figures", help="Figure directory")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)
    k = args.kind

    if k == "streamline":
        dyn.pl.streamline_plot(adata, basis=args.basis, color=args.color, save_show_or_return="return")
    elif k == "topography":
        dyn.pl.topography(adata, basis=args.basis, color=args.color, save_show_or_return="return")
    elif k == "phase":
        if not args.genes:
            die("--kind phase needs --genes")
        dyn.pl.phase_portraits(adata, genes=args.genes, save_show_or_return="return")
    elif k == "umap":
        dyn.pl.umap(adata, color=args.color, save_show_or_return="return")
    elif k == "kinetic_heatmap":
        if not args.genes:
            die("--kind kinetic_heatmap needs --genes")
        dyn.pl.kinetic_heatmap(adata, genes=args.genes, save_show_or_return="return")
    elif k in SCALAR_FIELDS:
        getattr(dyn.pl, k)(adata, basis=args.basis, save_show_or_return="return")
    else:  # unreachable given choices=, but explicit
        die(f"unsupported kind: {k}")

    save_fig(f"{args.figdir}/{k}.png")


if __name__ == "__main__":
    main()

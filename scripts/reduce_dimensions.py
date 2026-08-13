#!/usr/bin/env python3
"""
Compute a low-dimensional embedding and project RNA velocity onto it.

Two steps:
  1. dyn.tl.reduceDimension  -> UMAP (or t-SNE / PCA) embedding in obsm.
  2. dyn.tl.cell_velocities  -> projects the high-dim velocity onto that basis
                                so it can be visualized and used for the vector
                                field. ekey/vkey are auto-detected per modality
                                (spliced: M_s/velocity_S; labeling: M_t/velocity_T).

Run after dynamics.py.

Examples:
    python reduce_dimensions.py dyn.h5ad -o red.h5ad
    python reduce_dimensions.py dyn.h5ad -o red.h5ad --reduction umap --basis umap
    python reduce_dimensions.py dyn.h5ad -o red.h5ad --ekey M_t --vkey velocity_T
"""

import argparse

from _common import add_io_args, configure_dynamo, info, load_adata, save_adata


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="red.h5ad")
    p.add_argument("--reduction", default="umap", choices=["umap", "tsne", "psl", "sude"],
                   help="Non-linear reduction method (default: umap)")
    p.add_argument("--basis", default="umap",
                   help="Embedding basis to project velocity onto (default: umap)")
    p.add_argument("--n-pca-components", type=int, default=30,
                   help="PCs fed into the non-linear reduction (default: 30)")
    p.add_argument("--n-neighbors", type=int, default=30,
                   help="Neighbors for the embedding graph (default: 30)")
    p.add_argument("--ekey", default=None, help="Expression layer key (default: auto)")
    p.add_argument("--vkey", default=None, help="Velocity layer key (default: auto)")
    p.add_argument("--cores", type=int, default=1, help="Cores (t-SNE only)")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)

    info(f"[reduce] reduction={args.reduction}  n_pca={args.n_pca_components}  "
         f"n_neighbors={args.n_neighbors}")
    dyn.tl.reduceDimension(
        adata,
        reduction_method=args.reduction,
        n_pca_components=args.n_pca_components,
        n_neighbors=args.n_neighbors,
        cores=args.cores,
    )

    info(f"[cell_velocities] projecting velocity onto basis={args.basis}")
    dyn.tl.cell_velocities(adata, basis=args.basis, ekey=args.ekey, vkey=args.vkey)

    save_adata(adata, args.output)
    info("Next: vector_field.py to reconstruct the continuous vector field.")


if __name__ == "__main__":
    main()

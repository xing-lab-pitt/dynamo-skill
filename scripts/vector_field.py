#!/usr/bin/env python3
"""
Reconstruct a continuous vector field from projected cell velocities
(dyn.vf.VectorField, sparse VFC). This is the object that everything downstream
— differential geometry, fate prediction, least-action paths, perturbation —
depends on.

The field for a basis is stored in adata.uns['VecFld_<basis>'] (and, for the
default reconstruction, keyed for reuse by the downstream tools). With
--map-topography it also computes fixed points / potential during the fit.

Run after reduce_dimensions.py.

Examples:
    python vector_field.py red.h5ad -o vf.h5ad
    python vector_field.py red.h5ad -o vf.h5ad --basis umap --map-topography
    python vector_field.py red.h5ad -o vf.h5ad --basis pca   # for differential geometry
"""

import argparse

from _common import add_io_args, configure_dynamo, info, load_adata, save_adata


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="vf.h5ad")
    p.add_argument("--basis", default="umap",
                   help="Basis to reconstruct the field in (default: umap). "
                        "Use pca for high-dim differential geometry / perturbation.")
    p.add_argument("--map-topography", action="store_true",
                   help="Also compute fixed points / potential during reconstruction")
    p.add_argument("--pot-curl-div", action="store_true",
                   help="Also compute potential, curl and divergence (2D bases)")
    p.add_argument("--cores", type=int, default=1, help="Cores (default: 1)")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)

    # VectorField needs the velocity projected onto this basis (obsm['velocity_<basis>']).
    # reduce_dimensions.py projects one basis (usually umap); if the requested basis
    # was not projected yet (e.g. --basis pca for differential geometry), do it now.
    if f"velocity_{args.basis}" not in adata.obsm:
        info(f"[cell_velocities] projecting velocity onto basis={args.basis} (missing)")
        dyn.tl.cell_velocities(adata, basis=args.basis)

    info(f"[vector_field] basis={args.basis}  map_topography={args.map_topography}")
    dyn.vf.VectorField(
        adata,
        basis=args.basis,
        map_topography=args.map_topography,
        pot_curl_div=args.pot_curl_div,
        cores=args.cores,
    )

    save_adata(adata, args.output)
    info("Next: differential_geometry.py / fate.py / least_action.py / perturbation.py")


if __name__ == "__main__":
    main()

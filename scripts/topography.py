#!/usr/bin/env python3
"""
Inspect and curate the fixed points of a reconstructed vector field.

`dyn.vf.topography` finds fixed points by seeding n points and following the
field, so raw output over-detects: many "fixed points" are numerical artifacts
sitting off the data manifold or duplicated inside one attractor basin. Every
published topography figure therefore shows a curated handful — Qiu et al.
(Cell 2022) keep 6 of ~200. Curation is the normal path, not an edge case, but
vector_field.py --map-topography cannot do it: it neither exposes the seed
count nor lets you subset the result.

So: --list to see what was detected, --keep to choose. The kept indices are
recorded in uns[VecFld_<basis>]["kept_fixed_points"] so a figure can be traced
back to the selection that produced it.

ftype: -1 stable (attractor), 0 saddle, 1 unstable (repeller/source).

Examples:
    python topography.py vf.h5ad --basis umap --n 750 --list
    python topography.py vf.h5ad -o vf_fps.h5ad --basis umap --n 750 \\
        --keep 2 8 1 195 4 5 --color cell_type
"""

import argparse

import numpy as np

from _common import (add_io_args, configure_dynamo, die, has_vectorfield, info,
                     load_adata, save_adata, save_fig)

FTYPE_NAME = {-1: "stable", 0: "saddle", 1: "unstable"}


def fixed_points(adata, basis):
    """(Xss, ftype) from the field's uns entry, or die explaining what is missing."""
    key = f"VecFld_{basis}"
    if key not in adata.uns:
        die(f"uns['{key}'] missing — run vector_field.py --basis {basis} first")
    vf = adata.uns[key]
    if "Xss" not in vf or "ftype" not in vf:
        die(f"uns['{key}'] has no fixed points; run this script without --keep "
            "first (it calls dyn.vf.topography to detect them)")
    return np.atleast_2d(np.asarray(vf["Xss"])), np.atleast_1d(np.asarray(vf["ftype"]))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="topography.h5ad")
    p.add_argument("--basis", default="umap", help="Vector-field basis (default: umap)")
    p.add_argument("--n", type=int, default=None,
                   help="Seed points for dyn.vf.topography. Higher n finds more "
                        "candidates (and more artifacts). Omit to reuse points "
                        "already in the object.")
    p.add_argument("--list", action="store_true",
                   help="Print the detected fixed points and exit without writing")
    p.add_argument("--keep", nargs="+", type=int, default=None,
                   help="Indices (from --list) to retain; all others are dropped")
    p.add_argument("--color", nargs="+", default=None,
                   help="obs key(s) to color the topography plot by")
    p.add_argument("--markersize", type=int, default=500, help="Fixed-point marker size")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)
    if not has_vectorfield(adata, args.basis):
        die(f"no vector field for basis '{args.basis}' — run vector_field.py --basis {args.basis}")

    if args.n is not None:
        info(f"[topography] detecting fixed points with n={args.n}, basis={args.basis}")
        dyn.vf.topography(adata, n=args.n, basis=args.basis)

    Xss, ftype = fixed_points(adata, args.basis)
    info(f"[topography] {len(ftype)} fixed points detected")

    if args.list or args.keep is None:
        counts = {FTYPE_NAME.get(int(t), str(t)): int((ftype == t).sum())
                  for t in np.unique(ftype)}
        info(f"[topography] by type: {counts}")
        print(f"{'idx':>5s}  {'type':9s}  coordinates", flush=True)
        for i, (x, t) in enumerate(zip(Xss, ftype)):
            coords = ", ".join(f"{v:.3f}" for v in np.atleast_1d(x)[:4])
            print(f"{i:5d}  {FTYPE_NAME.get(int(t), str(t)):9s}  [{coords}]", flush=True)
        if args.keep is None:
            # Detection is the expensive half (minutes), curation is instant. Save
            # the detected-but-uncurated object so the --keep pass can reuse it
            # instead of re-running dyn.vf.topography.
            if args.n is not None and args.output:
                save_adata(adata, args.output)
                info(f"Detection cached. Now curate without re-detecting:\n"
                     f"       python topography.py {args.output} --basis {args.basis} "
                     f"--keep i j k --color cell_type")
            else:
                info("Pick the points that sit one-per-lineage and re-run with "
                     "--keep i j k ... to write a curated object.")
            return

    bad = [i for i in args.keep if i < 0 or i >= len(ftype)]
    if bad:
        die(f"--keep index out of range {bad}; there are {len(ftype)} fixed points (0-{len(ftype)-1})")

    vf = adata.uns[f"VecFld_{args.basis}"]
    vf["Xss"], vf["ftype"] = Xss[args.keep], ftype[args.keep]
    # Record the selection: a curated topography figure is otherwise impossible
    # to reproduce, since the indices depend on --n and on dynamo's seeding.
    vf["kept_fixed_points"] = list(map(int, args.keep))
    if args.n is not None:
        vf["kept_fixed_points_n"] = int(args.n)
    kept = [FTYPE_NAME.get(int(t), str(t)) for t in vf["ftype"]]
    info(f"[topography] kept {len(args.keep)} of {len(ftype)}: {args.keep} ({kept})")

    try:
        dyn.pl.topography(adata, markersize=args.markersize, basis=args.basis,
                          fps_basis=args.basis, color=args.color,
                          streamline_alpha=0.9, save_show_or_return="return")
        save_fig(f"{args.figdir}/topography_curated.png")
    except Exception as e:
        info(f"[topography] plot skipped ({type(e).__name__}: {e})")

    save_adata(adata, args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Predict cell fate by integrating cells forward/backward through the vector
field (dyn.pd.fate), and draw streamline + topography plots.

Pick the starting cells either explicitly (--init-cells barcode1 barcode2 ...)
or by an obs group (--init-group cell_type:HSC selects all cells whose
obs['cell_type'] == 'HSC'). Requires a vector field (vector_field.py first).

Examples:
    python fate.py vf.h5ad -o fate.h5ad --init-group cell_type:HSC --basis umap
    python fate.py vf.h5ad -o fate.h5ad --init-cells CELL_1 CELL_2 --direction forward
"""

import argparse

from _common import add_io_args, configure_dynamo, info, load_adata, save_adata, save_fig, select_cells


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="fate.h5ad")
    p.add_argument("--init-cells", nargs="+", default=None, help="Explicit starting cell names")
    p.add_argument("--init-group", default=None, help="'obs_column:value' selecting starting cells")
    p.add_argument("--basis", default="umap", help="Vector-field basis (default: umap)")
    p.add_argument("--direction", default="forward", choices=["forward", "backward", "both"],
                   help="Integration direction (default: forward)")
    p.add_argument("--interpolation-num", type=int, default=100,
                   help="Points along each predicted trajectory (default: 100)")
    p.add_argument("--color", default=None, help="obs/var key to color plots by")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)
    init_cells = select_cells(adata, args.init_cells, args.init_group, role="init")

    info(f"[fate] integrating basis={args.basis} direction={args.direction}")
    dyn.pd.fate(
        adata,
        init_cells=init_cells,
        basis=args.basis,
        direction=args.direction,
        interpolation_num=args.interpolation_num,
    )

    color = args.color or (args.init_group.split(":", 1)[0] if args.init_group else None)
    dyn.pl.streamline_plot(adata, basis=args.basis, color=color, save_show_or_return="return")
    save_fig(f"{args.figdir}/fate_streamline.png")
    try:
        dyn.pl.topography(adata, basis=args.basis, color=color, save_show_or_return="return")
        save_fig(f"{args.figdir}/fate_topography.png")
    except Exception as e:
        info(f"topography plot skipped ({e}); run vector_field.py --map-topography to enable")

    save_adata(adata, args.output)


if __name__ == "__main__":
    main()

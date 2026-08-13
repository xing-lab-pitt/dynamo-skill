#!/usr/bin/env python3
"""
Compute the least-action path (LAP) between two cell states — the most probable
transition / reprogramming route through the vector field (dyn.pd.least_action).

Select start and target cells by explicit names or by obs group
('column:value'). LAP is usually run in the PCA basis, so build the field with
`vector_field.py --basis pca` beforehand.

Examples:
    python least_action.py vf.h5ad -o lap.h5ad \\
        --init-group cell_type:HSC --target-group cell_type:Meg --basis pca
    python least_action.py vf.h5ad -o lap.h5ad \\
        --init-cells C1 C2 --target-cells C9 C10
"""

import argparse

from _common import add_io_args, configure_dynamo, info, load_adata, save_adata, save_fig, select_cells


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="lap.h5ad")
    p.add_argument("--init-cells", nargs="+", default=None, help="Explicit start cell names")
    p.add_argument("--init-group", default=None, help="'obs_column:value' for start cells")
    p.add_argument("--target-cells", nargs="+", default=None, help="Explicit target cell names")
    p.add_argument("--target-group", default=None, help="'obs_column:value' for target cells")
    p.add_argument("--basis", default="pca", help="Vector-field basis (default: pca)")
    p.add_argument("--n-points", type=int, default=25, help="Points along the path (default: 25)")
    p.add_argument("--n-cells", type=int, default=1,
                   help="Representative init/target cells to use (default: 1). LAP solves "
                        "one path per init cell, so a whole population is intractable; this "
                        "caps each endpoint to the first N selected cells.")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)
    init_cells = select_cells(adata, args.init_cells, args.init_group, role="init")
    target_cells = select_cells(adata, args.target_cells, args.target_group, role="target")

    # LAP is a per-init-cell optimization — cap endpoints to a few representative cells.
    if args.n_cells and len(init_cells) > args.n_cells:
        info(f"[least_action] capping init cells {len(init_cells)} -> {args.n_cells} (--n-cells)")
        init_cells = init_cells[: args.n_cells]
    if args.n_cells and len(target_cells) > args.n_cells:
        info(f"[least_action] capping target cells {len(target_cells)} -> {args.n_cells} (--n-cells)")
        target_cells = target_cells[: args.n_cells]

    info(f"[least_action] basis={args.basis}  n_points={args.n_points}")
    dyn.pd.least_action(
        adata,
        init_cells=init_cells,
        target_cells=target_cells,
        basis=args.basis,
        n_points=args.n_points,
    )

    # LAP stores its gene list as a pandas Index in uns, which the h5ad writer
    # cannot serialize nested in a dict; coerce it to a plain list.
    lap = adata.uns.get(f"LAP_{args.basis}")
    if isinstance(lap, dict) and hasattr(lap.get("genes"), "tolist"):
        lap["genes"] = lap["genes"].tolist()

    try:
        dyn.pl.least_action(adata, basis=args.basis, save_show_or_return="return")
        save_fig(f"{args.figdir}/least_action_path.png")
    except Exception as e:
        info(f"LAP plot skipped ({e})")

    save_adata(adata, args.output)
    info("Tip: pd.GeneLeastActionPath + pl.kinetic_heatmap show gene dynamics along the path.")


if __name__ == "__main__":
    main()

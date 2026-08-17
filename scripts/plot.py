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
               differential-geometry scalar fields (needs differential_geometry.py).
               --kind jacobian takes --genes as regulator/effector pairs.
  jacobian_heatmap  regulator x effector Jacobian matrix (needs --genes)
  scalar_by_group   distribution of a per-cell scalar split by an obs column —
               the way to compare e.g. speed across cell types, which painting
               the scalar on the embedding cannot show
  lap_paths    least-action paths overlaid on the field, colored by action
               (needs lap_matrix.py --save-paths)
  kinetic_heatmap  gene expression kinetics ordered along the field

Examples:
    python plot.py vf.h5ad --kind streamline --color cell_type
    python plot.py dyn.h5ad --kind phase --genes GATA1 KLF1
    python plot.py dg.h5ad --kind acceleration --basis umap
    python plot.py dg.h5ad --kind scalar_by_group --scalar speed_pca --group cell_type
    python plot.py dg.h5ad --kind jacobian --genes SPI1 GATA1
    python plot.py vf.h5ad --kind lap_paths --paths results/lap/laps.pkl \\
        --pairs "HSC->Meg" "HSC->Ery" --out-name lap_development
"""

import argparse
import json
import os
import pickle

from _common import configure_dynamo, die, info, load_adata, save_fig

SCALAR_FIELDS = {"acceleration", "curvature", "divergence", "speed", "jacobian"}


def plot_scalar_by_group(adata, scalar, group, palette=None):
    """Violin + box of a per-cell scalar, split by an obs column, sorted by median."""
    import matplotlib.pyplot as plt
    import numpy as np

    for col, what in ((scalar, "--scalar"), (group, "--group")):
        if col not in adata.obs.columns:
            die(f"obs['{col}'] not found ({what}). Available: "
                f"{', '.join(list(adata.obs.columns)[:25])}...")

    vals = adata.obs[scalar].astype(float)
    keys = adata.obs[group].astype(str)
    groups = sorted(keys.unique(), key=lambda g: float(np.median(vals[keys == g])))
    data = [vals[keys == g].values for g in groups]

    fig, ax = plt.subplots(figsize=(0.75 * len(groups) + 3, 4))
    parts = ax.violinplot(data, showextrema=False)
    if palette:
        for body, g in zip(parts["bodies"], groups):
            if g in palette:
                body.set_facecolor(palette[g])
                body.set_alpha(0.85)
    ax.boxplot(data, widths=0.15, showfliers=False,
               medianprops=dict(color="k"), whiskerprops=dict(color="k"))
    ax.set_xticks(range(1, len(groups) + 1), groups, rotation=45, ha="right")
    ax.set_ylabel(scalar)
    ax.set_xlabel(group)
    ax.set_title(f"{scalar} by {group} (sorted by median)")
    fig.tight_layout()

    order = ", ".join(f"{g}={np.median(vals[keys == g]):.4g}" for g in reversed(groups))
    info(f"[plot] {scalar} median, highest first: {order}")


def plot_lap_paths(dyn, adata, paths_file, pairs, basis, color):
    """Overlay saved least-action paths on the streamline field, colored by action.

    Path geometry comes from lap_matrix.py --save-paths; this only draws it, so
    the expensive optimization is paid once and can be re-plotted in any grouping.
    """
    import matplotlib.pyplot as plt
    from dynamo.plot.utils import map2color

    if not os.path.exists(paths_file):
        die(f"--paths {paths_file} not found; produce it with "
            "`lap_matrix.py --save-paths {paths_file}`")
    with open(paths_file, "rb") as fh:
        saved = pickle.load(fh)
    if not isinstance(saved, dict) or not saved:
        die(f"{paths_file} holds no paths")

    missing = [p for p in pairs if p not in saved]
    if missing:
        die(f"pairs not in {paths_file}: {missing}\n"
            f"       available: {', '.join(sorted(saved))}")

    # streamline_plot wants a single color key here, but --color is nargs="+"
    key = color[0] if isinstance(color, (list, tuple)) and color else (color or "cell_type")
    fig, ax = plt.subplots(figsize=(6, 4))
    axes = dyn.pl.streamline_plot(adata, basis=basis, color=key, frontier=True,
                                  save_show_or_return="return", ax=ax)
    ax = axes[0] if isinstance(axes, (list, tuple)) else axes

    drawn = 0
    for pair in pairs:
        rec = saved[pair]
        for prediction, action in zip(rec["prediction"], rec["action"]):
            xy = prediction[:, [0, 1]]
            ax.scatter(*xy.T, c=map2color(action), zorder=3, s=18)
            ax.plot(*xy.T, c="k", zorder=2, lw=1)
            drawn += 1
    info(f"[plot] drew {drawn} path(s) for {len(pairs)} pair(s)")
    fig.tight_layout()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Processed .h5ad")
    p.add_argument("--kind", required=True,
                   choices=["streamline", "topography", "phase", "umap",
                            "kinetic_heatmap", "jacobian_heatmap",
                            "scalar_by_group", "lap_paths"] + sorted(SCALAR_FIELDS))
    p.add_argument("--basis", default="umap", help="Embedding basis to draw on (default: umap)")
    p.add_argument("--compute-basis", default=None,
                   help="Basis the scalar/Jacobian was COMPUTED in, when it differs "
                        "from the embedding (differential geometry is normally run in "
                        "pca and viewed on umap). Default: same as --basis.")
    p.add_argument("--color", nargs="+", default=None, help="obs/var key(s) to color by")
    p.add_argument("--genes", nargs="+", default=None,
                   help="Genes (phase / kinetic_heatmap; regulators+effectors for "
                        "jacobian / jacobian_heatmap)")
    p.add_argument("--cell-idx", nargs="+", type=int, default=None,
                   help="Numeric cell indices for --kind jacobian_heatmap. Default: "
                        "average over all cells (one matrix) — without this it draws "
                        "one panel per cell.")
    p.add_argument("--average", action="store_true",
                   help="Average the Jacobian over the selected --cell-idx cells")
    p.add_argument("--scalar", default=None,
                   help="obs key holding the per-cell scalar (--kind scalar_by_group)")
    p.add_argument("--group", default=None,
                   help="obs key to split by (--kind scalar_by_group)")
    p.add_argument("--palette", default=None,
                   help="JSON file or inline JSON mapping group -> color "
                        "(--kind scalar_by_group)")
    p.add_argument("--paths", default=None,
                   help="Pickle from lap_matrix.py --save-paths (--kind lap_paths)")
    p.add_argument("--pairs", nargs="+", default=None,
                   help="Path keys to draw, e.g. 'HSC->Meg' (--kind lap_paths)")
    p.add_argument("--pairs-file", default=None,
                   help="File with one path key per line, instead of --pairs")
    p.add_argument("--out-name", default=None,
                   help="Figure basename (default: the --kind). Use it when one kind "
                        "produces several figures, e.g. lap_development vs lap_reprogram.")
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
    elif k == "jacobian_heatmap":
        if not args.genes:
            die("--kind jacobian_heatmap needs --genes (regulators/effectors)")
        # dyn.pl.jacobian_heatmap draws ONE PANEL PER CELL unless average=True, so
        # passing a whole dataset silently asks for thousands of panels. Default to
        # the population-averaged matrix; --cell-idx opts into specific cells.
        if args.cell_idx:
            cells, average = list(args.cell_idx), args.average
        else:
            cells, average = list(range(adata.n_obs)), True
            info(f"[plot] averaging the Jacobian over all {adata.n_obs} cells "
                 "(pass --cell-idx for per-cell panels)")
        dyn.pl.jacobian_heatmap(adata, cell_idx=cells, average=average,
                                regulators=args.genes, effectors=args.genes,
                                basis=args.basis, save_show_or_return="return")
    elif k == "scalar_by_group":
        if not args.scalar or not args.group:
            die("--kind scalar_by_group needs --scalar and --group")
        palette = None
        if args.palette:
            raw = (open(args.palette).read() if os.path.exists(args.palette)
                   else args.palette)
            palette = json.loads(raw)
        plot_scalar_by_group(adata, args.scalar, args.group, palette)
    elif k == "lap_paths":
        pairs = args.pairs
        if args.pairs_file:
            pairs = [ln.strip() for ln in open(args.pairs_file) if ln.strip()]
        if not args.paths or not pairs:
            die("--kind lap_paths needs --paths and --pairs/--pairs-file")
        plot_lap_paths(dyn, adata, args.paths, pairs, args.basis, args.color)
    elif k in SCALAR_FIELDS:
        cbasis = args.compute_basis or args.basis
        if k == "jacobian":
            # dyn.pl.jacobian already separates the embedding (basis) from where the
            # Jacobian was computed (j_basis), and needs --genes for a named pair.
            if not args.genes:
                die("--kind jacobian needs --genes (regulators/effectors)")
            dyn.pl.jacobian(adata, regulators=args.genes, effectors=args.genes,
                            basis=args.basis, j_basis=cbasis,
                            save_show_or_return="return")
        elif cbasis == args.basis:
            getattr(dyn.pl, k)(adata, basis=args.basis, save_show_or_return="return")
        else:
            # dyn.pl.speed/divergence/acceleration/curvature use one `basis` for BOTH
            # the obs key and the embedding, so they cannot show a pca-computed scalar
            # on a umap embedding — the normal case, since differential geometry is
            # run in pca. Paint the scalar on the requested embedding instead.
            key = f"{k}_{cbasis}"
            if key not in adata.obs.columns:
                die(f"obs['{key}'] not found — run differential_geometry.py "
                    f"--basis {cbasis} --quantities {k} first")
            info(f"[plot] {k} computed in '{cbasis}', drawn on '{args.basis}' "
                 f"(coloring by obs['{key}'])")
            dyn.pl.scatters(adata, basis=args.basis, color=key, frontier=True,
                            save_show_or_return="return")
    else:  # unreachable given choices=, but explicit
        die(f"unsupported kind: {k}")

    save_fig(f"{args.figdir}/{args.out_name or k}.png")


if __name__ == "__main__":
    main()

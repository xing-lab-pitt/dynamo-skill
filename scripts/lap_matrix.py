#!/usr/bin/env python3
"""
Least-action-path (LAP) action matrix across cell states.

least_action.py answers "what is the most probable route from A to B". This
script answers the follow-up that dynamo's LAP analysis is really for: how the
action compares across EVERY ordered pair of states. In Qiu et al. (Cell 2022)
that comparison is the result — developmental transitions cost less action than
de-differentiation or transdifferentiation, so the action matrix ranks how hard
each conversion is.

Solves one LAP per ordered pair of the selected groups and reports the final
action of each, as a matrix plus a source->target vs target->source comparison.

Endpoint choice: one representative cell per group, the cell CLOSEST TO ITS
GROUP'S CENTROID in the field's basis. (least_action.py's --n-cells takes the
first N cells of a group, which makes the answer depend on row order; for an
all-pairs matrix that would be reported as signal.) Use --endpoint first to
reproduce the old behaviour.

Run in the pca basis: the action is defined in the high-dimensional state space.
Build the field with `vector_field.py --basis pca` beforehand.

Cost: one optimization per ordered pair — N groups is N*(N-1) solves. Expect a
few minutes per pair on a few thousand cells; 6 groups (30 pairs) is hours, not
minutes. Start with 3-4 groups.

Examples:
    python lap_matrix.py vf_pca.h5ad -o lap.h5ad --group-col cell_type \\
        --source HSC --targets Meg Ery Bas Mon Neu --out-dir results/lap
    python lap_matrix.py vf_pca.h5ad -o lap.h5ad --group-col cell_type \\
        --groups HSC Meg Ery
"""

import argparse
import os

import numpy as np
import pandas as pd

from _common import (add_io_args, configure_dynamo, die, info, load_adata,
                     save_adata, save_fig)


def representative_cell(adata, col, group, basis, how):
    """One cell name standing in for a group (centroid-closest, or first)."""
    mask = (adata.obs[col].astype(str) == group).values
    if not mask.sum():
        return None
    names = adata.obs_names[mask]
    if how == "first":
        return names[0]
    key = f"X_{basis}"
    if key not in adata.obsm:
        die(f"obsm['{key}'] missing — needed for centroid endpoint selection")
    X = np.asarray(adata.obsm[key])[mask]
    return names[int(np.argmin(np.linalg.norm(X - X.mean(0), axis=1)))]


def pick_adj_key(adata, explicit=None):
    if explicit:
        if explicit not in adata.obsp:
            die(f"obsp['{explicit}'] not found; have {list(adata.obsp.keys())}")
        return explicit
    for k in ("pearson_transition_matrix", "cosine_transition_matrix"):
        if k in adata.obsp:
            return k
    die(f"no transition matrix in obsp; have {list(adata.obsp.keys())}. "
        "Run reduce_dimensions.py (cell_velocities) first.")


def final_action(adata, lap_obj, key):
    """Final action of a LAP: from uns if dynamo stored it, else from the object."""
    u = adata.uns.get(key, {})
    if isinstance(u, dict) and "action" in u:
        a = np.atleast_1d(np.asarray(u["action"], dtype=float)).ravel()
        if a.size:
            return float(a[-1])
    try:
        a = np.atleast_1d(np.asarray(lap_obj.action_t(), dtype=float)).ravel()
        return float(a[-1])
    except Exception:
        return float("nan")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="lap_matrix.h5ad")
    p.add_argument("--group-col", default="cell_type",
                   help="obs column holding the cell states (default: cell_type)")
    p.add_argument("--groups", nargs="+", default=None,
                   help="States to include. Default: every category in --group-col")
    p.add_argument("--source", default=None,
                   help="Origin state for the differentiation-vs-reverse table "
                        "(e.g. the progenitor). Added to --groups if missing.")
    p.add_argument("--targets", nargs="+", default=None,
                   help="Terminal states, used with --source. Shorthand for "
                        "--groups <source> <targets...>")
    p.add_argument("--basis", default="pca", help="Vector-field basis (default: pca)")
    p.add_argument("--n-points", type=int, default=25, help="Points along each path")
    p.add_argument("--endpoint", choices=["centroid", "first"], default="centroid",
                   help="Representative cell per group (default: centroid)")
    p.add_argument("--adj-key", default=None, help="obsp transition-matrix key")
    p.add_argument("--min-cells", type=int, default=50,
                   help="Warn when a group has fewer cells than this (default: 50)")
    p.add_argument("--out-dir", default="results/lap", help="CSV output directory")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)

    if args.group_col not in adata.obs.columns:
        die(f"obs column '{args.group_col}' not found")
    col = adata.obs[args.group_col].astype(str)

    if args.targets:
        groups = ([args.source] if args.source else []) + list(args.targets)
    elif args.groups:
        groups = list(args.groups)
        if args.source and args.source not in groups:
            groups = [args.source] + groups
    else:
        groups = list(pd.Index(col.unique()))
    groups = [g for g in dict.fromkeys(groups) if (col == g).sum() > 0]
    if len(groups) < 2:
        die(f"need >= 2 groups present in '{args.group_col}', got {groups}")

    counts = {g: int((col == g).sum()) for g in groups}
    info(f"[lap_matrix] basis={args.basis} groups={groups}")
    info(f"[lap_matrix] cells per group: {counts}")
    thin = [g for g, n in counts.items() if n < args.min_cells]
    if thin:
        info(f"[lap_matrix] WARNING low-n groups, paths weakly constrained: {thin}")

    adj_key = pick_adj_key(adata, args.adj_key)
    reps = {g: representative_cell(adata, args.group_col, g, args.basis, args.endpoint)
            for g in groups}
    info(f"[lap_matrix] endpoint selection: {args.endpoint}")
    info(f"[lap_matrix] solving {len(groups) * (len(groups) - 1)} ordered pairs "
         "(one optimization each) ...")

    M = pd.DataFrame(np.nan, index=groups, columns=groups, dtype=float)
    failures = []
    for src in groups:
        for tgt in groups:
            if src == tgt:
                continue
            key = f"LAP_{src}_to_{tgt}"
            try:
                lap = dyn.pd.least_action(
                    adata, init_cells=[reps[src]], target_cells=[reps[tgt]],
                    basis=args.basis, adj_key=adj_key, n_points=args.n_points,
                    add_key=key,
                )
                if isinstance(lap, list):
                    lap = lap[0]
                M.loc[src, tgt] = final_action(adata, lap, key)
            except Exception as e:  # one bad pair must not lose the whole matrix
                failures.append((src, tgt, type(e).__name__, str(e)[:120]))

    os.makedirs(args.out_dir, exist_ok=True)
    M.to_csv(os.path.join(args.out_dir, "lap_action_matrix.csv"))
    n_ok = int(M.notna().sum().sum())
    info(f"[lap_matrix] solved {n_ok}/{len(groups) * (len(groups) - 1)} pairs")
    for f in failures:
        info(f"[lap_matrix] FAILED {f[0]}->{f[1]}: {f[2]}: {f[3]}")

    # ---- matrix heatmap -----------------------------------------------------
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1.1 * len(groups) + 3, 1.0 * len(groups) + 2.5))
    im = ax.imshow(M.values, cmap="viridis")
    ax.set_xticks(range(len(groups)), groups, rotation=45, ha="right")
    ax.set_yticks(range(len(groups)), groups)
    ax.set_xlabel("target state")
    ax.set_ylabel("source state")
    ax.set_title("LAP action (lower = more probable transition)")
    for i in range(len(groups)):
        for j in range(len(groups)):
            v = M.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.3g}", ha="center", va="center", color="w", fontsize=7)
    fig.colorbar(im, ax=ax, label="action")
    fig.tight_layout()
    save_fig(f"{args.figdir}/lap_action_matrix.png")

    # ---- forward vs reverse from --source ----------------------------------
    if args.source and args.source in groups:
        others = [g for g in groups if g != args.source]
        cmp_df = pd.DataFrame(
            {f"{args.source}->x": M.loc[args.source, others].astype(float).values,
             f"x->{args.source}": M.loc[others, args.source].astype(float).values},
            index=others,
        )
        cmp_df["ratio_reverse_over_forward"] = cmp_df.iloc[:, 1] / cmp_df.iloc[:, 0]
        cmp_df.to_csv(os.path.join(args.out_dir, "lap_forward_vs_reverse.csv"))
        info("[lap_matrix] forward vs reverse action:")
        print(cmp_df.to_string(float_format=lambda x: f"{x:.5g}"), flush=True)

        x = np.arange(len(others))
        fig, ax = plt.subplots(figsize=(1.2 * len(others) + 3, 4))
        ax.bar(x - 0.2, cmp_df.iloc[:, 0], 0.4, label=f"{args.source} → x")
        ax.bar(x + 0.2, cmp_df.iloc[:, 1], 0.4, label=f"x → {args.source}")
        ax.set_xticks(x, others)
        ax.set_ylabel("LAP action")
        ax.set_title(f"Transitions out of and back into {args.source}")
        ax.legend()
        fig.tight_layout()
        save_fig(f"{args.figdir}/lap_action_barplot.png")

    # dyn.pd.least_action leaves pandas Index objects in uns; anndata cannot
    # serialize them, so coerce before writing (same fix as least_action.py).
    for k, v in list(adata.uns.items()):
        if isinstance(v, dict):
            for kk, vv in list(v.items()):
                if isinstance(vv, pd.Index):
                    v[kk] = list(vv)

    save_adata(adata, args.output)
    if n_ok == 0:
        die("no LAP solved — check the vector field and transition matrix")
    info("Tip: pd.GeneLeastActionPath + pl.kinetic_heatmap show genes along one path.")


if __name__ == "__main__":
    main()

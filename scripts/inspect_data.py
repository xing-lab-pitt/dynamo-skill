#!/usr/bin/env python3
"""
Inspect an AnnData file for a dynamo workflow and print a structured summary.

Reports shape, detected modality (splicing vs metabolic labeling), obs/var
columns, layers, obsm/uns keys, and — importantly — which dynamo pipeline steps
have already run (velocity, dimension reduction, cell velocities projection,
vector field, differential geometry). Use this before analysis to decide which
stages still need to run.

Examples:
    python inspect_data.py data.h5ad
    python inspect_data.py processed.h5ad --max-cols 60
"""

import argparse

from _common import configure_dynamo, detect_modality, load_adata


def _named(mapping):
    return [k for k in mapping.keys() if isinstance(k, str)]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="Input .h5ad (or other anndata-readable file)")
    p.add_argument("--max-cols", type=int, default=40, help="Max obs/var columns to list")
    args = p.parse_args()

    configure_dynamo()
    adata = load_adata(args.input)
    layers = _named(adata.layers)

    print("=" * 70)
    print(f"AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
    print("=" * 70)

    modality = detect_modality(adata)
    print(f"\nDetected modality: {modality}")
    print(f"layers: {layers}")

    print(f"\nobs columns ({len(adata.obs.columns)}):")
    for c in adata.obs.columns[:args.max_cols]:
        col = adata.obs[c]
        if str(col.dtype) in ("category", "object"):
            nuniq = col.nunique()
            extra = f"  {nuniq} categories" + (f": {list(col.unique()[:8])}" if nuniq <= 8 else "")
        else:
            try:
                extra = f"  range=[{col.min():.2f}, {col.max():.2f}]"
            except (TypeError, ValueError):
                extra = ""
        print(f"  - {c} ({col.dtype}){extra}")

    print(f"\nvar columns ({len(adata.var.columns)}): {list(adata.var.columns[:args.max_cols])}")
    print(f"\nobsm: {_named(adata.obsm)}")
    print(f"uns:  {_named(adata.uns)}")

    # Which dynamo steps have run?
    done = []
    if "X_pca" in adata.obsm:
        done.append("PCA")
    if any(k.startswith("velocity_") for k in layers) or "velocity_S" in layers or "velocity_T" in layers:
        done.append("velocity (dynamics)")
    embeddings = [k for k in _named(adata.obsm) if k.startswith("X_") and k != "X_pca"]
    if embeddings:
        done.append("dim-reduction " + str(embeddings))
    if any(k.startswith("velocity_") and k.endswith(("umap", "pca", "tsne")) for k in _named(adata.obsm)):
        done.append("cell_velocities (projected)")
    vecflds = [k for k in _named(adata.uns) if k.startswith("VecFld")]
    if vecflds:
        done.append("vector field " + str(vecflds))
    if any(k in adata.uns for k in ("jacobian_pca", "acceleration", "curvature", "divergence", "curl")):
        done.append("differential geometry")
    print("\nDynamo steps already present: " + (", ".join(done) if done else "none"))

    # Next-step hint
    if "velocity (dynamics)" not in " ".join(done):
        print("\nNext: preprocess.py -> dynamics.py")
    elif not embeddings:
        print("\nNext: reduce_dimensions.py")
    elif not vecflds:
        print("\nNext: vector_field.py")
    else:
        print("\nReady for downstream: differential_geometry.py / fate.py / least_action.py / perturbation.py")


if __name__ == "__main__":
    main()

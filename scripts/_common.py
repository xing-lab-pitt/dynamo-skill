#!/usr/bin/env python3
"""
Shared helpers for the dynamo script toolkit.

Every CLI script in this directory imports from this module so that data
loading, saving, figure configuration, modality detection, and logging behave
consistently. This file is NOT a CLI itself; import it:

    from _common import load_adata, save_adata, configure_dynamo, info, detect_modality
"""

import os
import sys


# ---- pandas >= 2 compatibility ----------------------------------------------
# dynamo was written against pandas < 2 and turns size-factor Series into column
# vectors with ``series[:, None]`` in ~130 places (e.g. tools/moments.py:724
# ``sfs, tot_sfs = sfs[:, None], sfs[:, None]``). pandas 2 rejects that:
#     ValueError: Multi-dimensional indexing (e.g. `obj[:, None]`) is no longer
#     supported. Convert to a numpy array before indexing instead.
# On pandas 2 this makes dynamics() fail outright for metabolic-labeling data,
# so every script in this toolkit needs the fix — hence applying it here, at
# import, rather than in each script.
#
# Safe by construction: it only intercepts the exact ``[:, None]`` key shape,
# which raises on pandas 2 today, so no working code path can change meaning.
# On pandas < 2 nothing is patched. Idempotent.

_NDIM_SHIM_FLAG = "_dynamo_skill_ndim_shim"


def _is_ndim_key(key):
    return (
        isinstance(key, tuple)
        and len(key) == 2
        and isinstance(key[0], slice)
        and key[1] is None
    )


def apply_pandas_ndim_shim():
    """Restore pre-2.0 ``series[:, None]`` behaviour. Returns #classes patched."""
    try:
        import numpy as np
        import pandas as pd
    except ImportError:  # pragma: no cover - pandas/numpy always present with dynamo
        return 0

    if tuple(int(p) for p in pd.__version__.split(".")[:1]) < (2,):
        return 0

    n = 0
    for cls in (pd.Series, pd.Index):
        if getattr(cls, _NDIM_SHIM_FLAG, False):
            continue
        orig = cls.__getitem__

        def _make(orig=orig):
            def __getitem__(self, key):
                if _is_ndim_key(key):
                    return np.asarray(self)[key[0], None]
                return orig(self, key)
            return __getitem__

        cls.__getitem__ = _make()
        setattr(cls, _NDIM_SHIM_FLAG, True)
        n += 1
    return n


apply_pandas_ndim_shim()
# -----------------------------------------------------------------------------


def info(msg):
    """Print a progress message with a marker."""
    print(f"[dynamo] {msg}", flush=True)


def die(msg, code=1):
    """Print an error and exit."""
    print(f"Error: {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


def _import_dynamo():
    try:
        import dynamo as dyn  # noqa: F401
        return dyn
    except ImportError:
        die("dynamo not installed. Activate the venv noted in SKILL.md, "
            "or install with: pip install dynamo-release")


def configure_dynamo(figdir="figures", dpi=120, background="white"):
    """Apply consistent dynamo/matplotlib settings and return the dynamo module.

    dynamo has no global figure directory (unlike scanpy). The scripts save
    figures explicitly with ``save_fig`` below, so this only sets the
    non-interactive backend, figure style, and creates ``figdir``.
    """
    import matplotlib
    matplotlib.use("Agg")  # headless: never try to open a window
    dyn = _import_dynamo()
    try:
        dyn.configuration.set_figure_params("dynamo", background=background, dpi=dpi)
    except Exception:  # keep going even if the style call changes across versions
        pass
    os.makedirs(figdir, exist_ok=True)
    return dyn


def save_fig(path, dpi=150):
    """Save the current matplotlib figure to ``path`` and close it.

    dynamo plotting functions draw onto the current figure when called with
    ``save_show_or_return="return"``; this captures that figure to disk with a
    predictable filename instead of using each function's bespoke save logic.
    """
    import matplotlib.pyplot as plt
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    info(f"Wrote figure {path}")


def load_adata(path):
    """Load an AnnData object, dispatching on the file extension.

      * ``.h5ad``          -> dyn.read_h5ad
      * ``.loom``          -> dyn.read_loom
      * anything else      -> dyn.read (delegates to anndata)
    """
    dyn = _import_dynamo()
    if not os.path.exists(path):
        die(f"input not found: {path}")
    lower = path.lower()
    if lower.endswith(".h5ad"):
        return dyn.read_h5ad(path)
    if lower.endswith(".loom"):
        return dyn.read_loom(path)
    return dyn.read(path)


def save_adata(adata, path):
    """Write an AnnData object to .h5ad, creating parent dirs as needed."""
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    adata.write_h5ad(path)
    info(f"Wrote {path}  ({adata.n_obs} cells x {adata.n_vars} genes)")


# ---- dynamo-specific helpers -------------------------------------------------

# Layers that identify each data modality. Splicing-based velocity needs
# spliced/unspliced counts; metabolic-labeling velocity needs new (labeled) and
# total counts (scNT-seq, scEU-seq, scSLAM-seq, ...).
_SPLICING_LAYERS = ("spliced", "unspliced")
_LABELING_LAYERS = ("new", "total")


def _named_layers(adata):
    return [k for k in adata.layers.keys() if isinstance(k, str)]


def detect_modality(adata):
    """Return 'labeling', 'splicing', or 'unknown' from the layers present.

    Labeling takes precedence when both label and splice layers exist, because
    the labeling kinetics model is the more informative one for such data.
    """
    layers = set(_named_layers(adata))
    has_label = all(k in layers for k in _LABELING_LAYERS)
    has_splice = all(k in layers for k in _SPLICING_LAYERS)
    if has_label:
        return "labeling"
    if has_splice:
        return "splicing"
    return "unknown"


def resolve_modality(adata, requested):
    """Turn a --modality flag ('auto'/'splicing'/'labeling') into a concrete one."""
    if requested and requested != "auto":
        return requested
    modality = detect_modality(adata)
    if modality == "unknown":
        die("could not detect modality: need spliced/unspliced (splicing) or "
            "new/total (labeling) layers. Pass --modality explicitly and add the "
            "required layers first.")
    info(f"Auto-detected modality: {modality}")
    return modality


def select_cells(adata, explicit=None, group=None, role="init"):
    """Resolve a set of cell names from an explicit list or an 'obs_col:value' group.

    Used by the fate / least-action scripts to pick starting and target cells.
    """
    if explicit:
        return list(explicit)
    if group:
        if ":" not in group:
            die(f"--{role}-group must be 'column:value' (e.g. cell_type:HSC)")
        col, val = group.split(":", 1)
        if col not in adata.obs.columns:
            die(f"obs column '{col}' not found (for --{role}-group)")
        cells = list(adata.obs_names[adata.obs[col].astype(str) == val])
        if not cells:
            die(f"no cells with obs['{col}'] == '{val}' (for --{role}-group)")
        info(f"[{role}] {len(cells)} cells from {col}=={val}")
        return cells
    die(f"provide --{role}-cells or --{role}-group")


def representative_cell(adata, col, group, basis="pca", how="centroid"):
    """One cell name standing in for a group: centroid-closest, or first.

    Used wherever a whole group must collapse to a single endpoint — LAP solves
    one path per init cell, so a whole population is intractable. "first"
    depends on row order, which turns an arbitrary ordering into apparent
    signal, so "centroid" is the default. least_action.py and lap_matrix.py both
    call this, so their endpoints agree by construction.

    Returns None when the group is absent; callers decide whether that is fatal.
    """
    import numpy as np

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


def has_vectorfield(adata, basis=None):
    """Whether a reconstructed vector field is present (optionally for a basis)."""
    key = "VecFld" if basis is None else f"VecFld_{basis}"
    return key in adata.uns


def add_io_args(parser, default_output=None):
    """Attach the standard input/output/figdir arguments to an argparse parser."""
    parser.add_argument("input", help="Input .h5ad (or .loom / other anndata-readable file)")
    parser.add_argument("-o", "--output", default=default_output,
                        help="Output .h5ad path" +
                             (f" (default: {default_output})" if default_output else ""))
    parser.add_argument("--figdir", default="figures",
                        help="Directory for saved figures (default: figures)")
    return parser

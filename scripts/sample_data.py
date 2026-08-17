#!/usr/bin/env python3
"""
Fetch one of dynamo's bundled sample datasets into a local .h5ad.

Most dynamo tutorials open by loading a published, already-processed object
(`dyn.sample_data.hematopoiesis()` and friends) rather than re-deriving one.
Those objects ship the paper's embedding, velocity, reconstructed vector field
and Jacobians, so reproducing a tutorial figure means plotting them, not
rebuilding them — and a local re-derivation can easily disagree with the
published result. This script is the toolkit's entry point for that.

Downloads are cached in --data-dir, so a repeat call is free.

Note `human_tfs` is a transcription-factor table, not an AnnData; write it with
a .csv/.tsv --output.

Examples:
    python sample_data.py --list
    python sample_data.py --dataset hematopoiesis -o data/hematopoiesis.h5ad
    python sample_data.py --dataset human_tfs -o data/human_tfs.csv
"""

import argparse
import inspect
import os

from _common import configure_dynamo, die, info, save_adata

# Callables in dynamo.sample_data that fetch a *dataset*; these three are the
# generic download plumbing those functions are built from, not datasets.
_HELPERS = {"download_data", "download_data_requests", "get_adata"}


def list_datasets(sd):
    """[(name, default filename or '', url or '')] for every bundled dataset."""
    out = []
    for name in sorted(dir(sd)):
        if name.startswith("_") or name in _HELPERS:
            continue
        fn = getattr(sd, name)
        if not inspect.isfunction(fn) or getattr(fn, "__module__", "") != sd.__name__:
            continue
        defaults = {}
        try:
            for p, v in inspect.signature(fn).parameters.items():
                if v.default is not inspect.Parameter.empty:
                    defaults[p] = v.default
        except (TypeError, ValueError):
            pass
        out.append((name, str(defaults.get("filename", "") or ""),
                    str(defaults.get("url", "") or "")))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list", action="store_true",
                   help="List the bundled datasets and whether each is cached, then exit")
    p.add_argument("--dataset", default=None, help="Dataset name (see --list)")
    p.add_argument("-o", "--output", default=None,
                   help="Output path (.h5ad; .csv/.tsv for table datasets)")
    p.add_argument("--data-dir", default=".",
                   help="Working directory for the download cache (default: '.'). "
                        "dynamo's loaders always cache into a 'data/' subdirectory "
                        "of it, so files land in <data-dir>/data/.")
    args = p.parse_args()

    dyn = configure_dynamo()
    sd = dyn.sample_data
    datasets = list_datasets(sd)

    cache_dir = os.path.join(args.data_dir, "data")
    if args.list or not args.dataset:
        info(f"[sample_data] {len(datasets)} datasets  (cache: {os.path.abspath(cache_dir)})")
        print(f"{'name':32s} {'cached':7s} {'filename':34s} source", flush=True)
        for name, fname, url in datasets:
            cached = bool(fname) and os.path.exists(os.path.join(cache_dir, fname))
            host = url.split("/")[2] if "//" in url else ""
            print(f"{name:32s} {'yes' if cached else '-':7s} {fname:34s} {host}", flush=True)
        if not args.dataset:
            if not args.list:
                die("pass --dataset NAME (or --list)")
            return

    names = [d[0] for d in datasets]
    if args.dataset not in names:
        die(f"unknown dataset '{args.dataset}'. Available: {', '.join(names)}")
    if not args.output:
        die("pass -o/--output")

    os.makedirs(args.data_dir, exist_ok=True)
    # dynamo's loaders cache relative to the working directory, so run the call
    # from --data-dir instead of threading a path through every loader signature.
    cwd = os.getcwd()
    out_abs = os.path.abspath(args.output)
    try:
        os.chdir(args.data_dir)
        info(f"[sample_data] fetching '{args.dataset}' (cached downloads are reused) ...")
        try:
            obj = getattr(sd, args.dataset)()
        except Exception as e:
            url = dict((n, u) for n, _, u in datasets).get(args.dataset, "")
            die(f"could not fetch '{args.dataset}': {type(e).__name__}: {e}\n"
                f"       source: {url or 'unknown'}\n"
                f"       If the host is unreachable, download it by hand into "
                f"{os.path.abspath('.')} and re-run.")
    finally:
        os.chdir(cwd)

    os.makedirs(os.path.dirname(out_abs) or ".", exist_ok=True)
    if hasattr(obj, "obs") and hasattr(obj, "var"):
        info(f"[sample_data] {args.dataset}: {obj.n_obs} cells x {obj.n_vars} genes")
        save_adata(obj, out_abs)
    elif hasattr(obj, "to_csv"):
        sep = "\t" if out_abs.endswith((".tsv", ".txt")) else ","
        obj.to_csv(out_abs, sep=sep, index=False)
        info(f"[sample_data] {args.dataset}: table {obj.shape} -> {out_abs}")
    else:
        die(f"'{args.dataset}' returned {type(obj).__name__}, which this script "
            "cannot write; load it in Python instead.")


if __name__ == "__main__":
    main()

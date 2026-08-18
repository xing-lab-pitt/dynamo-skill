#!/usr/bin/env python3
"""
Preprocess an scRNA-seq AnnData for dynamo: QC, gene selection, size-factor
normalization, PCA, and (optionally) cell-cycle scoring, via dyn.pp.Preprocessor.

Works for both modalities:
  * splicing  -> spliced/unspliced layers, conventional recipe.
  * labeling  -> new/total layers; pass --tkey and --experiment-type so the
                 preprocessor standardizes the labeling design before normalizing.

The Preprocessor recipe bundles the individual dyn.pp.* steps (calc_sz_factor,
select_genes, normalize, pca). Use --recipe to switch the gene-selection /
normalization strategy.

Examples:
    python preprocess.py raw.h5ad -o pp.h5ad
    python preprocess.py raw.h5ad -o pp.h5ad --recipe seurat --n-top-genes 3000
    python preprocess.py labeled.h5ad -o pp.h5ad --modality labeling \\
        --tkey time --experiment-type kin
"""

import argparse

from _common import (add_io_args, configure_dynamo, info, load_adata, resolve_modality,
                     save_adata, symbol_converter)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    add_io_args(p, default_output="pp.h5ad")
    p.add_argument("--recipe", default="monocle",
                   choices=["monocle", "seurat", "sctransform", "pearson_residuals",
                            "monocle_pearson_residuals"],
                   help="Preprocessing recipe (default: monocle)")
    p.add_argument("--modality", default="auto", choices=["auto", "splicing", "labeling"],
                   help="Data modality (default: auto-detect from layers)")
    p.add_argument("--n-top-genes", type=int, default=None,
                   help="Number of feature genes to select (default: 2000, or the "
                        "length of --force-gene-list* when one is given)")
    p.add_argument("--force-gene-list", default=None,
                   help="File with one gene per line to use INSTEAD of gene selection "
                        "(papers usually pin an author-curated list)")
    p.add_argument("--force-gene-list-uns", default=None,
                   help="uns key holding that gene list, e.g. genes_to_use")
    p.add_argument("--ensembl-release", type=int, default=None,
                   help="Ensembl release used to turn Ensembl IDs into gene symbols. "
                        "dynamo's own default is 77 (2014), whose symbols predate "
                        "current HGNC names; pin a release so the vintage is explicit "
                        "and matches whatever gene list you compare against")
    p.add_argument("--tkey", default=None,
                   help="obs column with labeling time (labeling data only)")
    p.add_argument("--experiment-type", default=None,
                   choices=["conventional", "kin", "deg", "one-shot", "mix_std_stm"],
                   help="Labeling experiment type (labeling data only)")
    p.add_argument("--no-cell-cycle", action="store_true",
                   help="Disable cell-cycle scoring (enabled by default)")
    args = p.parse_args()

    dyn = configure_dynamo(figdir=args.figdir)
    adata = load_adata(args.input)
    modality = resolve_modality(adata, args.modality)

    if modality == "labeling" and args.tkey is None:
        info("WARNING: labeling data usually needs --tkey (labeling time in obs); "
             "continuing without it.")

    if args.force_gene_list and args.force_gene_list_uns:
        p.error("use --force-gene-list or --force-gene-list-uns, not both")
    force_genes = None
    if args.force_gene_list:
        with open(args.force_gene_list) as fh:
            force_genes = [ln.strip() for ln in fh if ln.strip()]
    elif args.force_gene_list_uns:
        if args.force_gene_list_uns not in adata.uns:
            p.error(f"uns['{args.force_gene_list_uns}'] not found; "
                    f"available: {sorted(adata.uns)}")
        force_genes = list(adata.uns[args.force_gene_list_uns])

    # A forced list replaces gene selection, so it also sets how many genes to keep.
    n_top_genes = args.n_top_genes
    if n_top_genes is None:
        n_top_genes = len(force_genes) if force_genes else 2000

    kwargs = {"force_gene_list": force_genes} if force_genes else {}
    if args.ensembl_release is not None:
        kwargs["convert_gene_name_function"] = symbol_converter(args.ensembl_release)
    preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=not args.no_cell_cycle,
                                       **kwargs)

    # Tune the feature-gene count for whichever recipe is selected.
    if args.recipe == "monocle":
        preprocessor.config_monocle_recipe(adata, n_top_genes=n_top_genes)
    else:
        # Other recipes read n_top_genes from their select_genes kwargs.
        preprocessor.select_genes_kwargs["n_top_genes"] = n_top_genes

    info(f"[preprocess] recipe={args.recipe}  modality={modality}  "
         f"n_top_genes={n_top_genes}"
         + (f"  force_gene_list={len(force_genes)} genes" if force_genes else ""))
    preprocessor.preprocess_adata(
        adata,
        recipe=args.recipe,
        tkey=args.tkey,
        experiment_type=args.experiment_type,
    )

    save_adata(adata, args.output)
    info("Next: dynamics.py to estimate kinetics + RNA velocity.")


if __name__ == "__main__":
    main()

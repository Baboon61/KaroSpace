"""
Command-line interface for KaroSpace.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional


def _parse_section_rotations_arg(raw: str) -> Optional[Dict[str, float]]:
    text = str(raw or "").strip()
    if not text:
        return None

    rotations: Dict[str, float] = {}
    for item in text.split(","):
        token = item.strip()
        if not token:
            continue
        if ":" not in token:
            raise ValueError(
                "each section rotation must use section_id:angle format"
            )
        section_id, angle_text = token.split(":", 1)
        section_id = section_id.strip()
        angle_text = angle_text.strip()
        if not section_id:
            raise ValueError("section rotation entries must include a section_id")
        try:
            rotations[section_id] = float(angle_text)
        except ValueError as exc:
            raise ValueError(
                f"invalid angle for section {section_id!r}: {angle_text!r}"
            ) from exc

    return rotations or None


def _run_export_cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate HTML viewer for Xenium spatial transcriptomics data"
    )
    io_args = parser.add_argument_group("Input/output")
    dataset_args = parser.add_argument_group("Dataset loading and coordinates")
    metadata_args = parser.add_argument_group("Metadata and labels")
    viewer_args = parser.add_argument_group("Viewer layout")
    gene_args = parser.add_argument_group("Gene content and storage")
    pseudobulk_args = parser.add_argument_group("Pseudobulk DE")
    neighborhood_args = parser.add_argument_group("Neighborhoods and interactions")
    overlay_args = parser.add_argument_group("Images, deconvolution, and utilities")

    io_args.add_argument(
        "input",
        type=str,
        help="Path to input .h5ad file"
    )
    io_args.add_argument(
        "-o", "--output",
        type=str,
        default="karospace.html",
        help="Output HTML file path (default: karospace.html)"
    )
    viewer_args.add_argument(
        "--annotation",
        type=str,
        default="leiden",
        dest="annotation",
        help="Initial cell annotation column or gene (default: leiden)"
    )
    viewer_args.add_argument(
        "--additional-annotations",
        type=str,
        default="",
        dest="additional_annotations",
        help="Comma-separated extra obs columns to embed as selectable annotations "
             "(e.g. a second clustering). Needed to compare annotations in the River plot."
    )
    gene_args.add_argument(
        "--genes",
        type=str,
        default="",
        help="Comma-separated genes to preload for expression visualization. Significant DE genes are embedded automatically."
    )
    dataset_args.add_argument(
        "-g", "--groupby",
        type=str,
        default="sample_id",
        help="Column to group sections by (default: sample_id)"
    )
    dataset_args.add_argument(
        "--group-order",
        type=str,
        default="",
        help="Comma-separated section/group IDs to control section ordering."
    )
    dataset_args.add_argument(
        "--spatial-key",
        type=str,
        default="spatial",
        help="Key in obsm containing spatial coordinates (default: spatial)"
    )
    dataset_args.add_argument(
        "--spatial-x",
        type=str,
        default=None,
        help=(
            "Obs/metadata column to use as X coordinates. Must be used with "
            "--spatial-y; creates adata.obsm[spatial_key] before export."
        ),
    )
    dataset_args.add_argument(
        "--spatial-y",
        type=str,
        default=None,
        help=(
            "Obs/metadata column to use as Y coordinates. Must be used with "
            "--spatial-x; creates adata.obsm[spatial_key] before export."
        ),
    )
    metadata_args.add_argument(
        "--metadata-labels",
        type=str,
        default="",
        help=(
            "JSON object mapping metadata/obs column keys to display names in the viewer "
            '(example: {"sample_id":"Sample","last_score":"Disease score"}).'
        ),
    )
    metadata_args.add_argument(
        "--metadata-columns",
        type=str,
        default="",
        help=(
            "Comma-separated obs columns to use as section metadata and filter chips "
            "(e.g. strain,region,Batch,Slide). Empty uses loader defaults."
        ),
    )
    metadata_args.add_argument(
        "--metadata-value-order",
        type=str,
        default="",
        help='JSON object mapping metadata columns to ordered value lists (example: {"strain":["WT","KO"]}).',
    )
    metadata_args.add_argument(
        "--metadata-max-columns",
        type=int,
        default=None,
        help="Limit the number of metadata columns used, preserving order."
    )

    viewer_args.add_argument(
        "--min-panel-size",
        type=int,
        default=150,
        help="Minimum panel width in pixels (default: 150). Grid auto-adjusts columns."
    )
    viewer_args.add_argument(
        "--spot-size",
        type=str,
        default="auto",
        help="Default spot size. Use a positive number or 'auto' (default: auto)."
    )
    viewer_args.add_argument(
        "--downsample",
        type=int,
        default=None,
        help="Downsample to N cells per section (for large datasets)"
    )
    viewer_args.add_argument(
        "--title",
        type=str,
        default="KaroSpace",
        help="Page title"
    )
    viewer_args.add_argument(
        "--outlineby",
        dest="outline_by",
        type=str,
        default="course",
        help=(
            "Metadata column used to paint panel outlines. Use 'None' to disable outlines. "
            "When the column is embedded as metadata/annotation, outlines reuse that palette. (default: course)"
        )
    )
    viewer_args.add_argument(
        "--outline-by",
        dest="outline_by",
        type=str,
        help=argparse.SUPPRESS,
    )
    viewer_args.add_argument(
        "--viewer-info-html",
        type=str,
        default=None,
        help="HTML string shown in the viewer Info tab."
    )
    viewer_args.add_argument(
        "--viewer-info-html-file",
        type=str,
        default=None,
        help="Path to an HTML fragment file shown in the viewer Info tab."
    )
    gene_args.add_argument(
        "--gene-encoding",
        choices=["auto", "dense", "sparse"],
        default="auto",
        help="Gene vector encoding. 'sparse' stores only non-zero indices/values (smaller HTML for zero-inflated data). (default: auto)"
    )
    gene_args.add_argument(
        "--gene-value-encoding",
        choices=["uint16", "uint8"],
        default="uint16",
        help="Sidecar/package gene value encoding for binary shards. (default: uint16)"
    )
    gene_args.add_argument(
        "--gene-storage",
        choices=["embedded", "sidecar"],
        default="embedded",
        help="Store genes in the HTML (`embedded`) or write non-embedded genes to an auxiliary JSON sidecar (`sidecar`). (default: embedded)"
    )
    gene_args.add_argument(
        "--gene-aux-path",
        type=str,
        default=None,
        help="Optional output path for the gene sidecar JSON when --gene-storage sidecar."
    )
    gene_args.add_argument(
        "--gene-sidecar-shard-size",
        type=int,
        default=256,
        help="Number of genes/features per sidecar shard. (default: 256)"
    )
    gene_args.add_argument(
        "--modalities",
        type=str,
        default=None,
        help=(
            "Comma-separated list of modalities to export (e.g. 'rna,protein'). "
            "Defaults to all detected. Non-default modalities require --gene-storage sidecar."
        ),
    )
    gene_args.add_argument(
        "--gene-sparse-zero-threshold",
        type=float,
        default=0.8,
        help="Only used with --gene-encoding auto. Use sparse encoding when zero fraction >= threshold. (default: 0.8)"
    )
    neighborhood_args.add_argument(
        "--neighbor-permutations",
        type=str,
        default="auto",
        help="Neighbor enrichment permutation count. Use 0 to disable, or 'auto' (default) which disables for very large datasets."
    )
    neighborhood_args.add_argument(
        "--neighbor-stats-groupby",
        type=str,
        default="auto",
        help="Comma-separated obs columns to compute neighbor composition stats for. Use 'auto' (default) to match the initial annotation; empty disables."
    )
    pseudobulk_args.add_argument(
        "--pseudobulk",
        type=str,
        default="auto",
        help=(
            "Category pseudobulk DE mode. Use 'auto' to analyze the initial --annotation "
            "and --pseudobulk-additional-annotations, or 'None' to disable. (default: auto)"
        )
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-additional-annotations",
        type=str,
        default="",
        help=(
            "Comma-separated additional annotation columns to analyze when pseudobulk or "
            "interaction markers are enabled. The initial --annotation is included automatically."
        )
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-counts-layer",
        type=str,
        default="counts",
        help="AnnData layer containing raw counts for pseudobulk DE. Use 'None' for adata.X. (default: counts)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-min-replicates",
        type=int,
        default=2,
        help="Minimum paired replicates required for each pseudobulk contrast. (default: 2)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-min-pct-expressed",
        type=float,
        default=0.0,
        help="Minimum fraction of cells expressing a gene required in both compared groups. Values >1 are interpreted as percentages. (default: 0)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-p-adjust-method",
        choices=["fdr_bh", "bonferroni", "holm", "none"],
        default="fdr_bh",
        help="Multiple-testing correction method for pseudobulk p-values. (default: fdr_bh)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-padj-cutoff",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for volcano highlighting and DE table inclusion. (default: 0.05)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-log2fc-cutoff",
        type=float,
        default=0.5,
        help="Absolute log2FC cutoff for volcano highlighting and DE table inclusion. (default: 0.5)"
    )
    pseudobulk_args.add_argument(
        "--pseudobulk-deseq2-fit-type",
        choices=["parametric", "mean"],
        default="parametric",
        help="PyDESeq2 dispersion trend fit type. Use 'mean' to avoid parametric trend fallback warnings. (default: parametric)"
    )
    neighborhood_args.add_argument(
        "--neighbor-stats-seed",
        type=int,
        default=0,
        help="Random seed for neighbor enrichment permutations. (default: 0)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers",
        type=str,
        default="auto",
        help=(
            "Contact-conditioned pseudobulk marker mode. Use 'auto' to analyze the initial --annotation "
            "and --pseudobulk-additional-annotations, or 'None' to disable. (default: auto)"
        )
    )
    neighborhood_args.add_argument(
        "--interaction-markers-top-targets",
        type=int,
        default=8,
        help="Number of target cell types to evaluate per source for contact-conditioned markers. (default: 8)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers-top-genes",
        type=int,
        default=20,
        help="Number of top DE genes to keep per source-target interaction. (default: 20)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers-min-cells",
        type=int,
        default=30,
        help="Minimum cells required per replicate contact+ and contact- pseudobulk sample. (default: 30)"
    )
    neighborhood_args.add_argument(
        "--interaction-markers-min-neighbors",
        type=int,
        default=1,
        help="Minimum target neighbors required to classify source cells as contact+. (default: 1)"
    )
    overlay_args.add_argument(
        "--section-rotations",
        type=str,
        default="",
        help="Comma-separated section_id:angle pairs for initial per-section rotations with exact degree values (example: S1:37.5,S2:-90)."
    )
    gene_args.add_argument(
        "--gene-correlation-top-n",
        type=int,
        default=10,
        help="Number of top correlated genes to show per embedded gene in the discovery panel. Use 0 to disable. (default: 10)"
    )
    gene_args.add_argument(
        "--cluster-means-n-genes",
        type=int,
        default=500,
        help="Maximum embedded pseudobulk-DE genes to expose in category mean summaries. Use 0 to disable. (default: 500)"
    )
    gene_args.add_argument(
        "--spatial-variable-genes-n",
        type=int,
        default=200,
        help="Number of top variable genes to score with Moran's I spatial autocorrelation. Requires spatial graph in obsp. Use 0 to disable. (default: 200)"
    )
    viewer_args.add_argument(
        "--scalebar-unit",
        type=str,
        default="μm",
        help="Unit label for the scalebar (default: μm). Assumes spatial coordinates are in this unit."
    )
    overlay_args.add_argument(
        "--deconvolutions",
        type=str,
        default="",
        help='JSON object mapping deconvolution labels to obs/obsm keys.'
    )
    overlay_args.add_argument(
        "--section-images",
        type=str,
        default="",
        help="JSON object mapping section IDs to image paths or image specs."
    )
    overlay_args.add_argument(
        "--section-images-max-px",
        type=int,
        default=4096,
        help="Maximum image dimension when embedding section images. (default: 4096)"
    )

    args = parser.parse_args(argv)

    # Check input file
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if not input_path.suffix == ".h5ad":
        print(f"Warning: Expected .h5ad file, got: {input_path.suffix}", file=sys.stderr)

    # Import here to avoid slow startup for --help
    from .data_loader import load_spatial_data
    from .exporter import export_to_html

    neighbor_perms: Optional[int]
    if str(args.neighbor_permutations).lower() == "auto":
        neighbor_perms = None
    else:
        try:
            neighbor_perms = int(args.neighbor_permutations)
        except ValueError:
            print("Error: --neighbor-permutations must be an integer or 'auto'", file=sys.stderr)
            sys.exit(2)

    spot_token = str(args.spot_size).strip()
    if spot_token.lower() in {"auto", "adaptive", "density"}:
        spot_size_value = "auto"
    else:
        try:
            spot_size_value = float(spot_token)
        except ValueError:
            print("Error: --spot-size must be a positive number or 'auto'", file=sys.stderr)
            sys.exit(2)
        if spot_size_value <= 0:
            print("Error: --spot-size must be a positive number or 'auto'", file=sys.stderr)
            sys.exit(2)

    def _parse_csv(value: str):
        cleaned = [v.strip() for v in str(value).split(",") if v.strip()]
        return cleaned or None

    def _parse_json_object(value: str, option_name: str):
        text = str(value or "").strip()
        if not text:
            return None
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError as exc:
            print(f"Error: {option_name} must be valid JSON: {exc}", file=sys.stderr)
            sys.exit(2)
        if not isinstance(parsed, dict):
            print(f"Error: {option_name} must be a JSON object/dictionary", file=sys.stderr)
            sys.exit(2)
        return parsed

    def _parse_optional_layer(value: str):
        text = str(value or "").strip()
        if not text or text.lower() in {"none", "null"}:
            return None
        return text

    def _parse_auto_or_none(value: str, option_name: str) -> Optional[str]:
        text = str(value or "").strip().lower()
        if text in {"", "auto"}:
            return "auto"
        if text in {"none", "null"}:
            return None
        print(f"Error: {option_name} must be 'auto' or 'None'", file=sys.stderr)
        sys.exit(2)

    if str(args.neighbor_stats_groupby).lower() == "auto":
        neighbor_stats_groupby = [args.annotation]
    else:
        neighbor_stats_groupby = _parse_csv(args.neighbor_stats_groupby)
    pseudobulk_mode = _parse_auto_or_none(args.pseudobulk, "--pseudobulk")
    interaction_markers_mode = _parse_auto_or_none(args.interaction_markers, "--interaction-markers")
    pseudobulk_additional_annotations = _parse_csv(args.pseudobulk_additional_annotations)
    additional_annotations = _parse_csv(args.additional_annotations)
    genes = _parse_csv(args.genes)
    group_order = _parse_csv(args.group_order)
    metadata_columns = _parse_csv(args.metadata_columns)
    metadata_labels_raw = _parse_json_object(args.metadata_labels, "--metadata-labels")
    metadata_labels = {
        str(key): str(value)
        for key, value in (metadata_labels_raw or {}).items()
        if str(key).strip() and value is not None
    } or None
    metadata_value_order_raw = _parse_json_object(args.metadata_value_order, "--metadata-value-order")
    metadata_value_order = None
    if metadata_value_order_raw:
        metadata_value_order = {}
        for key, values in metadata_value_order_raw.items():
            if not isinstance(values, list):
                print("Error: --metadata-value-order values must be lists", file=sys.stderr)
                sys.exit(2)
            metadata_value_order[str(key)] = [str(v) for v in values]
    deconvolutions_raw = _parse_json_object(args.deconvolutions, "--deconvolutions")
    deconvolutions = {
        str(key): str(value)
        for key, value in (deconvolutions_raw or {}).items()
        if str(key).strip() and value is not None
    } or None
    section_images = _parse_json_object(args.section_images, "--section-images")
    viewer_info_html = args.viewer_info_html
    if args.viewer_info_html_file:
        try:
            viewer_info_html = Path(args.viewer_info_html_file).expanduser().read_text(encoding="utf-8")
        except OSError as exc:
            print(f"Error: could not read --viewer-info-html-file: {exc}", file=sys.stderr)
            sys.exit(2)
    try:
        section_rotations = _parse_section_rotations_arg(args.section_rotations)
    except ValueError as exc:
        print(f"Error: --section-rotations {exc}", file=sys.stderr)
        sys.exit(2)

    spatial_columns = None
    if bool(args.spatial_x) != bool(args.spatial_y):
        print("Error: --spatial-x and --spatial-y must be provided together", file=sys.stderr)
        sys.exit(2)
    if args.spatial_x and args.spatial_y:
        spatial_columns = (args.spatial_x, args.spatial_y)
    outline_token = str(args.outline_by or "").strip()
    outline_by = None if outline_token.lower() in {"", "none", "null"} else outline_token

    # Load and export
    print(f"Loading data from: {args.input}")
    load_kwargs = {
        "groupby": args.groupby,
        "spatial_key": args.spatial_key,
    }
    if spatial_columns is not None:
        load_kwargs["spatial_columns"] = spatial_columns
    if group_order is not None:
        load_kwargs["group_order"] = group_order
    if metadata_columns is not None:
        load_kwargs["metadata_columns"] = metadata_columns
    if metadata_value_order is not None:
        load_kwargs["metadata_value_order"] = metadata_value_order
    if args.metadata_max_columns is not None:
        load_kwargs["metadata_max_columns"] = args.metadata_max_columns
    dataset = load_spatial_data(args.input, **load_kwargs)

    modalities_arg: Optional[List[str]] = None
    if args.modalities:
        modalities_arg = [m.strip() for m in args.modalities.split(",") if m.strip()]

    print(f"Exporting to HTML...")
    output_path = export_to_html(
        dataset,
        output_path=args.output,
        annotation=args.annotation,
        additional_annotations=additional_annotations,
        genes=genes,
        title=args.title,
        modalities=modalities_arg,
        min_panel_size=args.min_panel_size,
        spot_size=spot_size_value,
        downsample=args.downsample,
        outline_by=outline_by,
        metadata_labels=metadata_labels,
        viewer_info_html=viewer_info_html,
        gene_encoding=args.gene_encoding,
        gene_value_encoding=args.gene_value_encoding,
        gene_storage=args.gene_storage,
        gene_aux_path=args.gene_aux_path,
        gene_sidecar_shard_size=args.gene_sidecar_shard_size,
        gene_sparse_zero_threshold=args.gene_sparse_zero_threshold,
        neighbor_stats_permutations=neighbor_perms,
        neighbor_stats_groupby=neighbor_stats_groupby,
        neighbor_stats_seed=args.neighbor_stats_seed,
        interaction_markers_top_targets=args.interaction_markers_top_targets,
        interaction_markers_top_genes=args.interaction_markers_top_genes,
        interaction_markers_min_cells=args.interaction_markers_min_cells,
        interaction_markers_min_neighbors=args.interaction_markers_min_neighbors,
        pseudobulk=pseudobulk_mode,
        pseudobulk_additional_annotations=pseudobulk_additional_annotations,
        pseudobulk_counts_layer=_parse_optional_layer(args.pseudobulk_counts_layer),
        pseudobulk_min_replicates=args.pseudobulk_min_replicates,
        pseudobulk_min_pct_expressed=args.pseudobulk_min_pct_expressed,
        pseudobulk_p_adjust_method=args.pseudobulk_p_adjust_method,
        pseudobulk_padj_cutoff=args.pseudobulk_padj_cutoff,
        pseudobulk_log2fc_cutoff=args.pseudobulk_log2fc_cutoff,
        pseudobulk_deseq2_fit_type=args.pseudobulk_deseq2_fit_type,
        interaction_markers=interaction_markers_mode,
        section_rotations=section_rotations,
        deconvolutions=deconvolutions,
        gene_correlation_top_n=args.gene_correlation_top_n,
        cluster_means_n_genes=args.cluster_means_n_genes,
        spatial_variable_genes_n=args.spatial_variable_genes_n,
        scalebar_unit=args.scalebar_unit,
        section_images=section_images,
        section_images_max_px=args.section_images_max_px,
    )

    if args.gene_storage == "sidecar":
        output_obj = Path(output_path).expanduser()
        print(
            "Done! Sidecar gene loading requires HTTP(S). "
            f"Serve the output directory with: python -m http.server --directory {output_obj.parent}"
        )
        print(f"Then open http://localhost:8000/{output_obj.name}")
    else:
        print(f"Done! Open {output_path} in a browser to view.")


def _run_package_sidecar_cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Package an existing KaroSpace sidecar viewer into a .karospace archive"
    )
    parser.add_argument(
        "html",
        type=str,
        help="Path to an existing sidecar HTML viewer",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output .karospace file path (default: <html stem>.karospace)",
    )
    parser.add_argument(
        "--gene-aux-path",
        type=str,
        default=None,
        help="Optional actual path to the sidecar gene manifest JSON if it differs from the path referenced in the HTML.",
    )
    parser.add_argument(
        "--gene-shard-dir",
        type=str,
        default=None,
        help="Optional actual path to the sidecar shard directory if it differs from the manifest stem directory.",
    )
    parser.add_argument(
        "--loader-output",
        type=str,
        default=None,
        help="Optional output path for the companion .loader.html file.",
    )

    args = parser.parse_args(argv)

    html_path = Path(args.html)
    if not html_path.exists():
        print(f"Error: Sidecar HTML not found: {args.html}", file=sys.stderr)
        sys.exit(1)

    from .exporter import package_sidecar_viewer

    print(f"Packaging existing sidecar viewer: {args.html}")
    package_path = package_sidecar_viewer(
        html_path,
        output_path=args.output,
        gene_manifest_path=args.gene_aux_path,
        gene_shard_dir=args.gene_shard_dir,
        loader_output_path=args.loader_output,
    )
    print(f"Done! Share {package_path} together with its .loader.html opener.")


def _run_ome_convert_cli(argv=None):
    parser = argparse.ArgumentParser(
        prog="karospace ome-convert",
        description=(
            "Optional helper: convert stitched TIFF(s) to pyramidal tiled "
            "OME-TIFF for import into Xenium Explorer. Independent of the HTML "
            "export — use it only if you need OME-TIFFs."
        ),
    )
    parser.add_argument("inputs", nargs="+", help="Input .tif file(s)")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: alongside each input)",
    )
    parser.add_argument(
        "--pyramid-levels",
        type=int,
        default=4,
        help="Number of pyramid levels (default: 4)",
    )
    args = parser.parse_args(argv)

    from .omeconvert import convert_to_ome

    for input_path in args.inputs:
        output_path = None
        if args.output_dir:
            p = Path(input_path)
            output_path = str(Path(args.output_dir) / (p.stem + ".ome.tif"))
        convert_to_ome(input_path, output_path, args.pyramid_levels)


def main():
    argv = list(sys.argv[1:])
    if argv and argv[0] in {"package-sidecar", "package"}:
        _run_package_sidecar_cli(argv[1:])
        return
    if argv and argv[0] in {"ome-convert", "omeconvert"}:
        _run_ome_convert_cli(argv[1:])
        return
    _run_export_cli(argv)


if __name__ == "__main__":
    main()

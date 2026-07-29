# KaroSpace

**KaroSpace** is an interactive HTML viewer for exploring spatial transcriptomics data. It generates standalone HTML files from h5ad files that can be shared and viewed in any web browser — no server or Python installation required.

Originally developed at Karolinska Institutet for visualizing Xenium spatial transcriptomics data across multiple tissue sections.

## Live Demo
- [KaroSpace website](https://karospace.se/)
- **Pancreas viewer**: [Open hosted demo](https://christoffermattssonlangseth.github.io/KaroSpace/pancreas.html)

## Features

- **Grid + modal exploration** — Browse many sections in a responsive grid, then zoom and pan any section in detail
- **Per-section rotation** — Set exact initial section angles at export time and adjust them interactively
- **Linked UMAP + section selection** — Magic Wand lasso works in UMAP and modal view with synced highlights
- **Selection summaries** — Selected-cell totals and per-type counts with expandable scrollable lists; minimize the panel to a compact header to keep the view clear
- **Polygon annotation workflow** — Save lasso selections as persistent annotations, reorder labels, and export JSON for downstream `adata` integration
- **Region-to-region DE** — Compare drawn annotations directly in the viewer, export JSON/CSV reports, and search top hits
- **Split compare slider** — Compare two variables side-by-side in the modal (`Cell type` or `Gene`, including `All categories`), draggable directly on the canvas
- **Legend controls + spotlight** — Toggle/hide categories and spotlight one class across grid and UMAP
- **Flexible coloring + gene discovery** — Switch between annotation columns, fuzzy-search genes, and reuse recent or saved gene panels
- **Insights panel** — `Summary`, `Compare`, `Genes`, `Neighborhood`, and `Regions` tabs with pseudobulk category DE, neighbor composition, interaction markers, and region comparison
- **Annotation river plot** — `Insights → Compare → River` draws a Sankey of how two annotations correspond (e.g. `leiden_1` ↔ `leiden_2`); click a node to recolor and spotlight it, or export the crosstab as CSV
- **Numeric category ordering** — Numeric cluster labels sort naturally (`2` before `10`) in legends, dropdowns, and plots, with `adata.uns` palettes kept aligned
- **Modal selection workflow** — Lasso in the sample view, save selections as annotations, open a focused subview, browse `Genes in selection`, and use `Space` + drag to pan while Select is active
- **Shareable packages** — Export as `.karospace` bundles (ZIP + viewer HTML) for offline sharing; open via the hosted loader at [karospace.se/open](https://karospace.se/open) or a local `loader.html`. See [KAROSPACE_PACKAGE_FORMAT_SPEC.md](KAROSPACE_PACKAGE_FORMAT_SPEC.md)
- **Compact sidecar options** — JSON and binary shard formats, sparse-first encoding, and optional `uint16`/`uint8` quantization for large datasets
- **Metadata-aware browsing** — Filter sections by metadata and outline by course or another column
- **Neighbor graph tools** — Graph overlay and hover rings (1–3 hops) when `adata.obsp` contains a spatial graph
- **Quality-of-life controls** — Hideable toolbar, screenshots, theme toggle, keyboard shortcuts, and adjustable spot size
- **Standalone export** — One self-contained HTML file, no backend required

## Browser Considerations

KaroSpace is canvas-heavy. Chrome/Chromium is generally fastest; Safari can be noticeably slower on large datasets. The viewer caps canvas DPR at `1.0` in Safari by default to reduce pixel work on Retina displays.

For large datasets:
- Use `downsample` to limit cells per section
- Lower `min_panel_size` to reduce pixels drawn per thumbnail
- Keep the neighbor graph toggle off unless needed

## Installation

### From source

```bash
git clone https://github.com/christoffermattssonlangseth/karospace.git
cd karospace
pip install -e .
```

### Dependencies

- Python >= 3.9
- scanpy >= 1.9.0
- anndata >= 0.8.0
- numpy >= 1.20.0
- pandas >= 1.3.0
- scipy >= 1.7.0

## Quick Start

### Desktop GUI (KaroSpaceBuilder)

Prebuilt executables are available from the
[`KaroSpaceBuilder` Releases](https://github.com/christoffermattssonlangseth/KaroSpaceBuilder/releases) page:
- Apple Silicon: `KaroSpaceBuilder-macos-arm64.zip`
- Windows: `KaroSpaceBuilder-windows.zip`
- Linux: `KaroSpaceBuilder-linux.zip`

Download, unzip, and run — no Python required.

To install from source:

```bash
git clone https://github.com/christoffermattssonlangseth/KaroSpaceBuilder
cd KaroSpaceBuilder
python -m pip install "git+https://github.com/christoffermattssonlangseth/KaroSpace.git"
python -m pip install -e .
KaroSpaceBuilder
```

If KaroSpaceBuilder is already installed, launch with:

```bash
karospacebuilder   # or: karospace-gui
```

**Workflow:**
1. Launch the app and pick a preset (`Default`, `Pancreas`, or `Lightweight`)
2. In `Basic`, set input `.h5ad` and output path, then click `Inspect`
3. In `Annotations & Genes`, build cell-annotation columns and gene lists with the searchable editors
4. Open `Advanced` only if needed, then export

Output is written to `<output_dir>/index.html`.

### Python API

```python
from karospace import load_spatial_data, export_to_html

dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",  # Column identifying each section
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)

export_to_html(
    dataset,
    output_path="viewer.html",
    annotation="cell_type",           # Initial cell-annotation column
    title="KaroSpace",
    min_panel_size=150,          # Min panel width (responsive autoscaling)
    spot_size="auto",            # Adaptive by section density (or set a fixed number)
    downsample=30000,            # Max cells per section
    additional_annotations=[          # Extra columns for annotation dropdown
        "leiden",
        "condition",
    ],
    genes=[                      # Pre-load genes for expression view
        "Cd4",
        "Cd8a",
        "Gfap",
    ],
    gene_encoding="auto",        # "auto" | "dense" | "sparse"
    gene_storage="embedded",     # "embedded" | "sidecar"
    gene_aux_path=None,          # Optional manifest path; defaults to viewer.genes.json
    gene_sparse_zero_threshold=0.8,
    neighbor_stats_groupby=["cell_type"],
    neighbor_stats_permutations=20,
    pseudobulk="auto",           # Use None to disable category pseudobulk DE
    pseudobulk_additional_annotations=["cell_type"],
    pseudobulk_counts_layer="counts",
    pseudobulk_min_pct_expressed=0.0,
    pseudobulk_p_adjust_method="fdr_bh",
    pseudobulk_padj_cutoff=0.05,
    pseudobulk_log2fc_cutoff=0.5,
    pseudobulk_deseq2_fit_type="parametric",
    interaction_markers="auto",  # Use None to disable contact-conditioned marker DE
    section_rotations={
        "sample_a": 37.5,
        "sample_b": -90,
    },
)
```

### Command Line

```bash
karospace your_data.h5ad -o viewer.html --annotation leiden
```

#### CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output HTML file path | `karospace.html` |
| `--annotation` | Initial cell-annotation column | `leiden` |
| `--additional-annotations` | Comma-separated extra obs columns to embed as selectable annotations (needed to compare two annotations in the River plot) | empty |
| `--genes` | Comma-separated genes to preload; significant pseudobulk DE genes are embedded automatically | empty |
| `--metadata-labels` | JSON object mapping metadata/obs column keys to display labels in the viewer UI | empty |
| `--metadata-columns` | Comma-separated obs columns to use as section metadata and filter chips | loader defaults |
| `--metadata-value-order` | JSON object mapping metadata columns to ordered value lists | empty |
| `--metadata-max-columns` | Limit metadata columns used, preserving order | empty |
| `-g, --groupby` | Column to group sections by | `sample_id` |
| `--group-order` | Comma-separated section/group IDs to control section order | empty |
| `--spatial-key` | Key in `adata.obsm` containing spatial coordinates, or target key created from `--spatial-x/--spatial-y` | `spatial` |
| `--spatial-x` | Obs/metadata column to use as X coordinates; requires `--spatial-y` | empty |
| `--spatial-y` | Obs/metadata column to use as Y coordinates; requires `--spatial-x` | empty |
| `--min-panel-size` | Minimum panel width in pixels | `150` |
| `--spot-size` | Cell/spot size (`auto` or positive number) | `auto` |
| `--downsample` | Max cells per section | None |
| `--title` | Page title | `KaroSpace` |
| `--outlineby` | Metadata column used to paint panel outlines; use `None` to disable. If the same column is in metadata, outlines reuse the metadata/annotation palette | `course` |
| `--viewer-info-html` | HTML string shown in the viewer Info tab | default info |
| `--viewer-info-html-file` | Path to an HTML fragment shown in the viewer Info tab | empty |
| `--gene-encoding` | Gene vector encoding (`auto`, `dense`, `sparse`) | `auto` |
| `--gene-value-encoding` | Sidecar/package gene value encoding for binary shards (`uint16`, `uint8`) | `uint16` |
| `--gene-storage` | Gene storage mode (`embedded`, `sidecar`) | `embedded` |
| `--gene-aux-path` | Path for the gene sidecar manifest JSON | auto |
| `--gene-sidecar-shard-size` | Genes/features per sidecar shard | `256` |
| `--gene-sparse-zero-threshold` | Zero fraction threshold for `auto` sparse encoding | `0.8` |
| `--modalities` | Comma-separated modalities to export | all detected |
| `--neighbor-permutations` | Permutations for neighbor enrichment z-scores | `auto` |
| `--neighbor-stats-groupby` | Obs columns for neighbor composition stats (`auto` or comma-separated) | `auto` |
| `--neighbor-stats-seed` | Random seed for neighbor enrichment permutations | `0` |
| `--interaction-markers` | Contact-conditioned pseudobulk marker mode (`auto`, `None`) | `auto` |
| `--interaction-markers-top-targets` | Target cell types evaluated per source for contact-conditioned markers | `8` |
| `--interaction-markers-top-genes` | Top DE genes kept per source-target interaction | `20` |
| `--interaction-markers-min-cells` | Minimum cells per replicate contact+ and contact- pseudobulk sample | `30` |
| `--interaction-markers-min-neighbors` | Minimum target neighbors to classify contact+ source cells | `1` |
| `--pseudobulk` | Category pseudobulk DE mode (`auto`, `None`) | `auto` |
| `--pseudobulk-additional-annotations` | Additional annotation columns to analyze when pseudobulk or interaction markers are enabled. The initial `--annotation` is included automatically | empty |
| `--pseudobulk-counts-layer` | Raw-count AnnData layer for pseudobulk aggregation; use `None` for `adata.X` | `counts` |
| `--pseudobulk-min-replicates` | Minimum paired replicates required for each contrast | `2` |
| `--pseudobulk-min-pct-expressed` | Minimum fraction of cells expressing a gene required in both compared groups; values >1 are interpreted as percentages | `0` |
| `--pseudobulk-p-adjust-method` | Multiple-testing correction method (`fdr_bh`, `bonferroni`, `holm`, `none`) | `fdr_bh` |
| `--pseudobulk-padj-cutoff` | Adjusted p-value cutoff for volcano highlighting and DE table inclusion | `0.05` |
| `--pseudobulk-log2fc-cutoff` | Absolute log2FC cutoff for volcano highlighting and DE table inclusion | `0.5` |
| `--pseudobulk-deseq2-fit-type` | PyDESeq2 dispersion trend fit type; use `mean` to avoid parametric trend fallback warnings | `parametric` |
| `--section-rotations` | Comma-separated `section_id:angle` pairs | empty |
| `--gene-correlation-top-n` | Correlated genes shown per embedded gene in discovery panel | `10` |
| `--cluster-means-n-genes` | Maximum embedded pseudobulk-DE genes used for category mean summaries; use `0` to disable | `500` |
| `--spatial-variable-genes-n` | Top variable genes scored with Moran's I; use `0` to disable | `200` |
| `--deconvolutions` | JSON object mapping deconvolution labels to obs/obsm keys | empty |
| `--section-images` | JSON object mapping section IDs to image paths/specs | empty |
| `--section-images-max-px` | Maximum image dimension when embedding section images | `4096` |
| `--scalebar-unit` | Unit label for the scalebar | `μm` |

## Data Requirements

- **`adata.obsm['spatial']`** — 2D coordinates for each cell (x, y)
- **`adata.obs[groupby]`** — Column identifying which section each cell belongs to
- **Categorical or numeric columns in `adata.obs`** — For assigning cell annotations and visualizing cells

If coordinates are stored as separate obs columns instead of an `obsm` matrix,
pass them on the CLI:

```bash
karospace your_data.h5ad -o viewer.html --spatial-x x_centroid --spatial-y y_centroid
```

This creates `adata.obsm["spatial"]` during loading. Use `--spatial-key` to pick
a different target key.

### Optional metadata

- `course` — Experimental phase (e.g., `"naive"`, `"peak_I"`); sections are outlined by this column when present
- `region`, `condition`, `timepoint` — Used for filter chips

Control display order of metadata values and section ordering via `metadata_value_order`:

```python
dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)
```

### Optional category palettes

If `adata.uns["{col}_colors"]` exists (scanpy convention — list of hex aligned to `adata.obs[col].cat.categories`), KaroSpace uses it for that column everywhere (legend, spots, neighbor views, samples panel). Length mismatch or missing key falls back to the default palette.

### Optional neighborhood graph

If `adata.obsp` contains a neighbor graph (`spatial_connectivities`, `connectivities`, `neighbors`, or `neighbor_graph`), KaroSpace exposes graph overlay and neighbor-hover controls.

Contact-conditioned interaction markers are computed automatically for the initial `color` column and for `pseudobulk_additional_annotations` unless `interaction_markers=None` / `--interaction-markers None` is used. KaroSpace classifies source cells as contact+ or contact- within each `groupby` replicate, aggregates raw counts by replicate/contact status, and fits pseudobulk DESeq2 with a paired replicate design.

Pseudobulk category-vs-category DE is precomputed automatically for the initial `color` column unless `pseudobulk=None` / `--pseudobulk None` is used, and shown in `Insights → Compare → Cell DE`. It aggregates raw counts by `groupby` replicate and category. Significant genes passing both adjusted p-value and log2FC thresholds are embedded automatically, exposed in `Insights → Genes → DE Genes`, and reused for category means/correlations. Use `pseudobulk_additional_annotations=[...]` or `--pseudobulk-additional-annotations ...` to compute category DE and interaction pseudobulk DE for extra annotation columns.

## Examples

See [`examples/`](examples/) for complete dataset-specific export scripts.

## Viewer Controls

### Grid View
- **Click a section** — Open detailed modal view
- **Annotation dropdown** — Switch between annotation columns
- **Gene input** — Fuzzy search with keyboard navigation, recent genes, saved panels, and pseudobulk DE gene suggestions
- **Size slider** — Adjust spot size
- **Filter chips** — Filter sections by metadata
- **Legend items** — Toggle categories; spotlight one across grid and UMAP
- **Insights button** — Toggle the insights panel (`Summary`, `Compare`, `Genes`, `Neighborhood`, `Regions`)
- **Screenshot / Theme** — Download snapshot or toggle dark/light mode

### Modal View
- **Scroll** — Zoom in/out
- **Drag** — Pan; `Space` + drag to pan while Select is active
- **Split button** — Open the A/B comparison panel; drag the divider line directly on the canvas
- **Magic Wand** — Draw lasso selection; opens `Genes in selection` ranked gene panel
- **Annotate** — Save lasso selections as persistent annotations; export as JSON with cell indices
- **Pick type** — Click a cell to spotlight its category
- **Focused view** — Open the current selection as a filtered subview
- **Hide tools** — Collapse the toolbar for an unobstructed canvas
- **Escape or click outside** — Close modal

### UMAP View
- **UMAP button** — Toggle panel; dock to any corner
- **Magic Wand** — Lasso selection synced to spatial views
- **Scroll / drag** — Zoom and pan the embedding

## Deployment and Sharing

- **Single-file HTML** — Default `gene_storage="embedded"` writes one standalone file; works offline
- **Sidecar mode** — `gene_storage="sidecar"` writes HTML + `<name>.genes.json` + `<name>.genes/` shards; requires HTTP(S) serving
- **`.karospace` packages** — Self-contained ZIP bundles combining the viewer HTML and all sidecar assets; open via the hosted loader or a local `loader.html`. All processing happens in the browser — data never leaves the user's machine
- **Static hosting** — Both embedded and sidecar modes work on GitHub Pages, S3, or any lab intranet
- **Sharing** — Send the HTML directly (embedded), or share HTML + manifest + shard directory via a hosted location (sidecar)

### Package an existing sidecar into `.karospace`

```bash
# Short form (auto-detects sidecar paths from the HTML)
karospace package-sidecar BALO.html --output BALO.karospace

# Explicit form
karospace package-sidecar BALO.html \
  --output BALO.karospace \
  --gene-aux-path BALO.genes.json \
  --gene-shard-dir BALO.genes \
  --loader-output BALO.loader.html
```

This wraps existing sidecar assets into a `.karospace` archive without recomputing analytics or rewriting the viewer HTML.

### Integrate polygon annotations back into AnnData

```python
import scanpy as sc
from karospace import integrate_polygon_annotations

adata = sc.read_h5ad("your_data.h5ad")
integrate_polygon_annotations(
    adata,
    "karospace-annotations-2026-02-12T12-00-00-000Z.json",
    label_key="lesion_labels",
    count_key="lesion_label_count",
    uns_key="lesion_polygons",
)
adata.write_h5ad("your_data_with_polygons.h5ad")
```

## License

MIT License

## Author

Christoffer Mattsson Langseth — Karolinska Institutet

<p align="center">
  <img src="assets/logo.png" alt="KaroSpace Logo" width="800">
</p>

**KaroSpace** is an interactive HTML viewer for exploring spatial transcriptomics data. It generates standalone HTML files from h5ad files that can be shared and viewed in any web browser without requiring a server or Python installation.

Originally developed at Karolinska Institutet for visualizing Xenium spatial transcriptomics data across multiple tissue sections.

## Live Demo

- **Pancreas viewer (GitHub Pages)**: [Open hosted demo](https://christoffermattssonlangseth.github.io/KaroSpace/pancreas.html)

## Features

- **Multi-section grid view** - Display dozens or hundreds of tissue sections in a responsive grid layout
- **Interactive zoom and pan** - Click any section to open a detailed view with mouse wheel zoom and drag-to-pan
- **UMAP view with Magic Wand selection** - Toggle UMAP panel, draw to select cells, highlights sync across views
- **Modal section Magic Wand** - Draw lasso selections directly inside the zoomed section view
- **Selection summaries** - Selected-cell totals + per-type counts with expandable category list and scrollable panel
- **Polygon annotation export** - Draw multiple persistent polygons per section and export as JSON for downstream `adata` integration
- **Modal split compare slider** - In zoomed section view, compare variable A vs B with a left/right split
- **Flexible split inputs** - A/B can be cell type/cluster columns or genes; cell-type sides support one category or `All categories`
- **Category toggling** - Click legend items to show/hide specific cell types or clusters; hidden cells appear as grey
- **Linked spotlight mode** - Spotlight a legend category across grid + UMAP for faster visual comparison
- **Gene expression visualization** - Pre-load genes of interest and switch between them with a viridis colormap
- **Multiple color columns** - Switch between different annotation columns (e.g., cell types, clusters, conditions)
- **Insights panel** - Search color columns, view categorical stats, and marker genes by color
- **Neighbor stats + enrichment** - Neighbor composition and permutation z-scores for categorical cell types
- **Cell-cell interaction browser** - Pick a source cell type to rank neighboring targets with enrichment, type markers, and contact-conditioned markers
- **Metadata filtering** - Filter sections by metadata like experimental condition, timepoint, or region
- **Cell tooltips** - Hover over cells to see their type or expression value
- **Course-based borders** - Section panels are outlined with colors indicating their experimental course/condition
- **Neighborhood graph overlay** - Toggle adjacency edges when a neighbor graph is present in `adata.obsp`
- **Neighbor rings on hover** - Highlight 1–3 hop neighbors around a hovered cell (when neighbors are available)
- **Screenshot export** - Download a full-page image of the current view
- **Dark/light mode** - Toggle between themes with preference saved to browser localStorage
- **Adjustable spot size** - Control cell/spot size with slider and +/- steps in grid/UMAP/modal views
- **Standalone HTML** - Generated files are self-contained with embedded data and JavaScript

## Browser Considerations

KaroSpace is canvas-heavy. Chrome/Chromium browsers are generally fastest; Safari can be noticeably slower on large datasets.
To mitigate this, the viewer caps canvas device pixel ratio (DPR) in Safari at `1.0` by default, reducing pixel work on Retina displays.

Performance tips for large datasets:
- Reduce `downsample` to limit cells per section.
- Lower `min_panel_size` to reduce the number of pixels drawn.
- Keep the neighbor graph toggle off unless needed.

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

### Python API

```python
from karospace import load_spatial_data, export_to_html

# Load your h5ad file
dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",  # Column identifying each section
    # Optional: custom ordering for metadata values (affects filter chips + outlines)
    # If group_order isn't set, the first key here is also used to order sections.
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)

# Export to HTML
export_to_html(
    dataset,
    output_path="viewer.html",
    color="cell_type",           # Initial color column
    title="KaroSpace",
    min_panel_size=150,          # Min panel width (responsive autoscaling)
    spot_size=2.0,               # Cell/spot size
    downsample=30000,            # Max cells per section (for large datasets)
    theme="light",               # "light" or "dark"
    additional_colors=[          # Extra columns for color dropdown
        "leiden",
        "condition",
    ],
    genes=[                      # Pre-load genes for expression view
        "Cd4",
        "Cd8a",
        "Gfap",
    ],
    # For zero-inflated expression matrices, store gene vectors sparsely (smaller HTML)
    gene_encoding="auto",        # "auto" | "dense" | "sparse"
    gene_sparse_zero_threshold=0.8,
    pack_arrays=True,            # Pack coords/colors/UMAP as base64 typed arrays (smaller + faster load)
    pack_arrays_min_len=1024,
    use_hvgs=True,               # Use adata.var['highly_variable'] when available
    hvg_limit=20,                # Max number of HVGs to include
    marker_genes_groupby=None,   # Marker genes off by default (enable with a categorical column)
    marker_genes_top_n=30,       # Top N markers per group
    neighbor_stats_groupby=[         # Compute neighbor composition stats for these categorical obs columns
        "cell_type",
    ],
    neighbor_stats_permutations=20,   # Permutations for neighbor enrichment z-scores (0 disables)
    neighbor_stats_seed=0,       # Random seed for permutations
    interaction_markers_groupby=None,           # Interaction markers off by default
    interaction_markers_top_targets=8,          # Targets evaluated per source type
    interaction_markers_top_genes=20,           # Genes shown per source-target interaction
    interaction_markers_min_cells=30,           # Minimum cells in contact+ and contact- groups
    interaction_markers_min_neighbors=1,        # Source cell needs >= this many target neighbors to be contact+
    interaction_markers_method="wilcoxon",      # DE method for contact+ vs contact-
    interaction_markers_layer="normalized",     # Layer used for DE (falls back to adata.X if missing)
)
```

### Command Line

```bash
karospace your_data.h5ad -o viewer.html --color leiden --cols 6
```

#### CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output HTML file path | `karospace.html` |
| `-c, --color` | Initial color column | `leiden` |
| `-g, --groupby` | Column to group sections by | `sample_id` |
| `--min-panel-size` | Minimum panel width in pixels (responsive autoscaling) | `150` |
| `--spot-size` | Cell/spot size | `2.0` |
| `--downsample` | Max cells per section | None |
| `--theme` | Color theme (`light` or `dark`) | `light` |
| `--title` | Page title | `KaroSpace` |
| `--gene-encoding` | Gene vector encoding (`auto`, `dense`, `sparse`) | `auto` |
| `--gene-sparse-zero-threshold` | Zero fraction threshold for `auto` sparse encoding | `0.8` |
| `--no-pack-arrays` | Disable base64 packing of large per-section arrays | off |
| `--pack-arrays-min-len` | Only pack arrays when section cell count ≥ this value | `1024` |
| `--neighbor-permutations` | Permutations for neighbor enrichment z-scores (`auto`, `0`, `100`, …) | `auto` |
| `--neighbor-stats-groupby` | Comma-separated obs columns to compute neighbor composition stats for (`auto`, empty disables) | `auto` |
| `--marker-genes-groupby` | Comma-separated obs columns to compute marker genes for (empty disables) | empty |
| `--interaction-markers-groupby` | Comma-separated obs columns to compute contact-conditioned markers for (empty disables) | empty |

## Data Requirements

Your h5ad file should have:

- **`adata.obsm['spatial']`** - 2D coordinates for each cell (x, y)
- **`adata.obs[groupby]`** - Column identifying which section each cell belongs to
- **Categorical or numeric columns in `adata.obs`** - For coloring cells (e.g., cell types, clusters)

### Optional metadata columns

For filtering functionality, include these columns in `adata.obs`:
- `course` - Experimental course/phase (e.g., "peak_I", "peak_II", "naive")
- `region` - Tissue region
- `condition` - Experimental condition
- `timepoint` - Time point

Sections will be outlined with colors based on their `course` metadata if present.

### Optional metadata ordering

You can control the display order of metadata values (filter chips + outline legend)
and section ordering by passing `metadata_value_order` to `load_spatial_data`.
If `group_order` is not provided, sections are ordered by the first key in
`metadata_value_order` (unknown values appear after the custom list).

```python
dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",
    metadata_value_order={
        "course": ["naive", "peak_I", "peak_II", "peak_III"],
    },
)
```

### Optional neighborhood graph

If `adata.obsp` contains a neighbor graph (e.g., `spatial_connectivities`, `connectivities`,
`neighbors`, or `neighbor_graph`), KaroSpace will expose graph/neighbor-hover controls.
Neighbor composition stats are **enabled by default** for the initial color. To customize,
pass `neighbor_stats_groupby=[...]` to `export_to_html` (or `--neighbor-stats-groupby`
in the CLI).

To enable contact-conditioned interaction markers in the interaction browser, pass
`interaction_markers_groupby=[...]` to `export_to_html` (typically the same categorical
column used for coloring, e.g. `cell_type`).

Interaction markers are defined per source-target pair as:
- **contact+**: source cells with at least `interaction_markers_min_neighbors` neighbors of the target type
- **contact-**: source cells of the same source type with zero neighbors of that target type

### Optional cell polygons

If you have per-cell polygons, store them in `adata.uns["polygons"]` using a flat vertex
buffer with offsets so each cell can have a variable number of vertices:

```python
# n_cells = adata.n_obs
# vertices is a flat (M, 2) array of x/y polygon points for all cells
# offsets is length n_cells + 1, with vertices for cell i in
# vertices[offsets[i]:offsets[i+1]]

adata.uns["polygons"] = {
    "vertices": vertices,  # shape (M, 2), float32/float64
    "offsets": offsets,    # shape (n_cells + 1,), int64
}
```

You can keep `adata.obsm["spatial"]` as the cell centroid coordinates for fallback rendering.

## Example

See [example.py](example.py) for a complete working example.

```python
from karospace import load_spatial_data, export_to_html

dataset = load_spatial_data(
    "your_data.h5ad",
    groupby="sample_id",
)

print(f"Loaded {dataset.n_sections} sections with {dataset.n_cells:,} cells")

export_to_html(
    dataset,
    output_path="viewer.html",
    color="anno_L2",
    title="KaroSpace",
    min_panel_size=120,          # Min panel width (responsive autoscaling)
    spot_size=1.5,
    downsample=30000,
    additional_colors=['anno_L3', 'anno_L2', 'anno_L1', 'leiden'],
    genes=["Cd4", "Cd8a", "Gfap", "Mbp"],
    use_hvgs=False,
    hvg_limit=20,
    marker_genes_groupby=["anno_L2"],
    marker_genes_top_n=30,
)
```

## Viewer Controls

### Grid View
- **Click a section** - Open detailed modal view
- **Color dropdown** - Switch between annotation columns
- **Gene input** - Type a gene name to view expression (must be pre-loaded)
- **Size slider** - Adjust spot size (drag or +/- buttons)
- **Filter chips** - Click to filter sections by metadata
- **Legend items** - Click to toggle categories on/off
- **Spotlight button (legend)** - Enable linked spotlight; hover/click legend categories to dim others across grid + UMAP
- **Legend button** - Show/hide the legend panel
- **Insights button** - Toggle the insights panel (colors + stats + marker genes)
- **Insights usage** - Use “Stats” to aggregate categorical colors by metadata, “Neighbors” for neighbor composition + z-scores + interaction browser (Contact markers + Type markers), and “Marker genes” for top markers
- **Graph button** - Toggle neighborhood graph overlay (if available)
- **Neighbors button** - Toggle neighbor rings on hover (if available)
- **Hop selector** - Choose which neighbor hop(s) to display (if available)
- **Screenshot button** - Download a snapshot of the current view
- **Theme button** - Toggle dark/light mode
- **Footer badge** - “KaroSpace” label with a GitHub link

### Modal View (Detailed Section)
- **Mouse wheel** - Zoom in/out
- **Click and drag** - Pan around
- **Zoom buttons** - +/- zoom controls
- **Reset button** - Return to default zoom/pan
- **Split button** - Open the A/B comparison panel in modal view
- **A/B variable selectors** - Choose per side: `Cell type` or `Gene`
- **Cell type mode** - Pick annotation column + either one category or `All categories`
- **Gene mode** - Pick a pre-loaded gene from the gene list
- **Split slider** - Controls left/right display share (e.g. 10% = left 10% uses A, right 90% uses B)
- **Magic Wand** - Draw a lasso directly on the section to select cells of interest
- **Annotate button** - Draw persistent polygons in the section (multiple polygons supported)
- **Annotation panel** - Rename polygons, select polygon cells, delete polygons, or clear section/all
- **Export JSON** - Download all polygon annotations including vertices + mapped cell indices
- **Screenshot button** - Download a PNG of the current sample (modal) view
- **Clear selection** - Remove current selected cells
- **Graph button** - Toggle neighborhood graph overlay (if available)
- **Neighbors button** - Toggle neighbor rings on hover (if available)
- **Hop selector** - Choose which neighbor hop(s) to display (if available)
- **Size slider** - Adjust spot size for this view
- **Selection summary** - Shows selected-cell count and per-type counts for the active annotation
- **Selection summary expand** - Click `+N more categories` to expand; click again to collapse
- **Selection summary scroll** - Scroll inside the summary panel to inspect long category lists
- **Hover over cells** - See cell type or expression value
- **Escape or click outside** - Close modal

### UMAP View (if available)
If your h5ad file contains UMAP coordinates (`adata.obsm['X_umap']`), a UMAP toggle button appears:
- **UMAP button** - Toggle the UMAP panel on/off
- **Dock button (TR/TL/BR/BL)** - Cycle the panel corner placement
- **Panel +/-** - Resize the UMAP panel
- **Magic Wand** - Activate lasso selection mode
- **Draw selection** - Click and drag to draw a selection area
- **Clear** - Clear the current selection
- **Size slider** - Adjust point size in the UMAP view
- **Mouse wheel** - Zoom the UMAP view
- **Selection summary** - Shows selected-cell count and per-type counts synced with the section view
- **Selection summary expand** - Click `+N more categories` to expand; click again to collapse
- **Selection summary scroll** - Scroll inside the summary panel to inspect long category lists

## Implementation, Deployment, and Use

### Implementation
- **Single-file HTML**: `export_to_html` writes one standalone HTML file that embeds all data (JSON) and viewer logic (CSS/JS) directly in the document.
- **No backend required**: All interactions (filtering, coloring, legends, zoom/pan) run entirely in the browser.
- **Data pipeline**: `load_spatial_data` reads `.h5ad`, builds section metadata, and serializes coordinates, colors, and optional gene vectors into JSON.

### Deployment
- **Local use**: Open the generated `.html` in any modern browser (Chrome, Firefox, Safari). No server needed.
- **Static hosting**: You can host the HTML as a static asset (e.g., GitHub Pages, S3, lab intranet). Since it’s self-contained, there are no runtime dependencies.
- **GitHub Pages (this repo)**: Use **Source = GitHub Actions** (static HTML artifact deploy). Do **not** use Jekyll/deploy-from-branch mode for this workflow.
- **Pancreas deployment workflow**: `.github/workflows/pages-balo.yml` publishes `pancreas.html` and creates an `index.html` redirect.
- **Pancreas test URL**: `https://christoffermattssonlangseth.github.io/KaroSpace/pancreas.html` (available after a successful Pages workflow run).
- **File size note**: Large datasets create large HTML files. Consider `downsample` and limiting `genes`/`additional_colors` to keep the file manageable.

### Use
- **Sharing**: Send the HTML file directly to collaborators; it will work offline once downloaded.
- **Viewing**: Just double-click or drag into a browser. All controls work immediately.
- **Updates**: Re-run `export_to_html` to refresh the file when annotations or metadata change.

Selected cells are highlighted with a yellow/gold outline in both UMAP and spatial views.

### Integrate Polygon Annotations Back Into AnnData

Polygon exports use global cell indices when available, so you can map annotations directly back to the original `adata`.

```python
import scanpy as sc
from karospace import integrate_polygon_annotations

adata = sc.read_h5ad("your_data.h5ad")
integrate_polygon_annotations(
    adata,
    "karospace-annotations-2026-02-12T12-00-00-000Z.json",
    label_key="lesion_labels",           # per-cell joined labels
    count_key="lesion_label_count",      # number of polygons covering each cell
    uns_key="lesion_polygons",           # full polygon metadata
)

adata.write_h5ad("your_data_with_polygons.h5ad")
```

## Performance Tips

- Use `downsample` parameter for datasets with >50,000 cells per section
- Limit `genes` list to only essential genes (each adds to file size)
- If you enable `use_hvgs`, the viewer preloads up to 20 HVGs to limit file size
- Consider splitting very large datasets into multiple viewers

## License

MIT License

## Author

Christoffer Mattsson Langseth - Karolinska Institutet

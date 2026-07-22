# KaroSpace Feature Summary

This document summarizes the current capabilities in this repository, focused on what the generated HTML viewer can do and how it is configured.

## 1. Product Surfaces

- Python API for loading `.h5ad` and exporting a standalone HTML viewer.
- Command-line interface (`karospace`) for scriptable exports.
- Desktop GUI (`karospace.gui`) for non-code export configuration.
- Self-contained HTML output for sharing and interactive exploration in a browser.

## 2. Data and Export Inputs

- Input: AnnData (`.h5ad`) with spatial coordinates in `adata.obsm`.
- Section grouping by any obs column (`groupby`).
- Optional metadata columns for filter chips and panel metadata.
- Optional metadata value ordering for deterministic section/filter order.
- Optional UMAP support from `adata.obsm["X_umap"]`.
- Optional neighborhood graph support from `adata.obsp` keys (for graph overlays and neighbor analytics).

## 3. HTML Export Characteristics

- Single-file HTML output with embedded data payload.
- Optional sidecar gene storage (`gene_storage="sidecar"`) with manifest + shard directory for lazy loading.
- Per-section arrays are packed into base64 typed arrays by default for smaller/faster files.
- Gene vectors support `dense`, `sparse`, or `auto` encoding.
- Gene loading can be manual (`genes=[...]`) and significant pseudobulk DE genes are embedded automatically.
- Gene discovery can be enriched with DE-gene and pseudobulk/category correlation suggestions (`gene_correlation_top_n`).
- Adaptive spot size mode (`spot_size="auto"`), or fixed numeric spot size.
- Optional exact initial per-section rotations via `section_rotations={section_id: angle}`.

## 4. Grid Viewer Features

- Responsive multi-section grid rendering.
- Per-section rotation controls with shared section state between grid and modal views.
- Metadata filter chips with reset.
- Dynamic layout behavior:
  - Single true section dataset gets centered single-section layout.
  - Filtered-to-one-section view gets capped-width layout (prevents full-page stretch).
- Category legend with hide/show toggles and spotlight behavior.
- Light viewer styling is built in.
- Screenshot export for current grid view.
- Category/gene color switching from the top controls.
- Gene discovery panel with fuzzy search, keyboard navigation, recent genes, saved panels, and pseudobulk-DE-driven suggestions.
- Performance-aware incremental rendering and offscreen skipping.

## 5. Modal (Section Detail) Features

- Click any section to open detailed modal view.
- Pan/zoom controls with mouse and buttons.
- Rotation controls that preserve exact stored per-section angles.
- Graph and neighbor hover controls when neighbor graph data exists.
- Magic Wand lasso selection and linked selection summaries.
- Polygon annotation workflow:
  - Draw persistent polygon annotations.
  - Export annotations as JSON.
  - Clear section or clear all annotations.
- Compact grouped sample-view toolbar:
  - Hide/show toggle for the full bottom toolbar.
  - Draggable dock layout.
  - Contextual control visibility for selection/type/graph actions.
- Modal screenshot export.

## 6. Split (Variable Slider) Features

- A/B split rendering inside modal:
  - Each side can be `cell` (categorical column + category) or `gene`.
  - Split slider controls left/right blend ratio.
- Per-side gene scale controls (manual min/max and auto percentile).
- Split gene views lazy-load sidecar genes when needed, so comparison mode can recover without manually reselecting the gene.
- Split legend in modal reflects active A/B variables.
- Pseudobulk DE genes are integrated into split metadata when side is `cell`:
  - If side category is `All categories`, DE-gene lists are shown across categories.
  - If a specific category is selected, only that category's DE genes are shown.

## 7. UMAP Features

- Optional UMAP panel with toggle button when UMAP exists.
- Dock positions (corner cycling), panel sizing controls.
- UMAP pan/zoom and point-size control.
- Linked lasso selection with grid/modal synchronization.

## 8. Insights Panel Features

### 8.1 Stats Tab

- Aggregate summaries by selected metadata column.
- Expand/collapse grouped summaries.
- Cell-type trend panel with search.

### 8.2 Neighbors Tab

- Neighbor composition stats by categorical color.
- Optional permutation-based enrichment z-scores.
- Target search/filter.
- Interaction browser:
  - Source/target interaction summaries.
  - Contact-conditioned DE genes (if precomputed).
  - Type DE genes per target.

### 8.3 Genes Tab

- Dotplot:
  - Group by categorical color.
  - Optional metadata-based aggregation filter.
  - Gene list input (comma-separated).
  - Dot size as fraction expressing; dot color as mean expression.
  - Dotplot gene input supports datalist-based suggestions and Tab token autocomplete.
- DE Genes:
  - Pseudobulk DE genes per category for current categorical color.
  - Search/filter DE-gene output.
- Means:
  - Pseudobulk-derived category mean expression for embedded DE genes.
- Discovery:
  - Fuzzy-searchable top-bar gene picker.
  - Recent genes and saved gene panels stored per viewer.
  - DE-gene and correlated-gene suggestions when available.
- Compare:
  - Pairwise pseudobulk category-vs-category DE for precomputed categorical columns.
  - Select source category and reference category within an annotation.
  - Ranked result table with log fold change, adjusted p-value, score, and percent expressing.
  - Clicking a DE gene activates that gene in the viewer.

## 9. Annotation Integration Back to AnnData

- Utility function: `integrate_polygon_annotations(...)`.
- Imports exported polygon JSON and maps annotation labels back to AnnData obs/uns.
- Supports global cell index mapping when present.

## 10. Interfaces and Configuration

- API entry points:
  - `load_spatial_data(...)`
  - `export_to_html(...)`
  - `integrate_polygon_annotations(...)`
- CLI supports core export settings:
  - color, groupby, panel size, spot size, downsample
  - gene encoding options
  - sidecar gene storage options
  - neighbor stats controls
  - pseudobulk DE controls
  - interaction pseudobulk controls
  - section rotations
- GUI supports:
  - Searchable list editors for additional colors and genes
  - Advanced controls for gene storage, pseudobulk DE, neighbor stats, and interaction markers
  - Inspect/validate + export workflow with logs

## 11. Reliability and UX Improvements Included

- Loading warning no longer persists indefinitely:
  - startup warning is removable and auto-dismissed
  - loader cleanup occurs even if init throws
- Main and modal spot-size sliders are constrained to max `5`.
- Dotplot genes autocomplete support added.
- Compare-tab gene clicks now keep the DE view working without a `refreshInsights` runtime error.
- Sample-view toolbar tightened and made draggable without shrinking text.
- Example script (`examples/CODEX.py`) now forces local package import and prints correct output filename.

## 12. Practical Notes

- Generated HTML is static and self-contained. Re-export is required to pick up new code/features.
- Old exported files will not gain new behavior automatically.
- Large datasets may produce very large HTML files; packing, sparse encoding, and downsampling are the primary size/performance levers.

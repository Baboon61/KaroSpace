"""
Data loading utilities for spatial transcriptomics data.

Handles loading h5ad files with scanpy and extracting spatial coordinates,
gene expression, and metadata for visualization.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import scipy.sparse as sp
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import json


@dataclass
class SectionData:
    """Data for a single tissue section."""
    section_id: str
    coordinates: np.ndarray  # (n_cells, 2) array of x, y coordinates
    metadata: Dict[str, str] = field(default_factory=dict)

    @property
    def n_cells(self) -> int:
        return self.coordinates.shape[0]

    @property
    def bounds(self) -> Tuple[float, float, float, float]:
        """Return (xmin, xmax, ymin, ymax)."""
        return (
            float(self.coordinates[:, 0].min()),
            float(self.coordinates[:, 0].max()),
            float(self.coordinates[:, 1].min()),
            float(self.coordinates[:, 1].max()),
        )


@dataclass
class SpatialDataset:
    """Container for spatial transcriptomics dataset."""
    adata: sc.AnnData
    sections: List[SectionData]
    groupby: str
    obs_columns: List[str]
    var_names: List[str]
    metadata_columns: List[str]
    metadata_value_order: Optional[Dict[str, List[str]]] = None

    @property
    def n_sections(self) -> int:
        return len(self.sections)

    @property
    def n_cells(self) -> int:
        return self.adata.n_obs

    @property
    def has_umap(self) -> bool:
        """Check if UMAP coordinates are available."""
        return "X_umap" in self.adata.obsm

    def get_color_data(
        self,
        color: str,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None
    ) -> Tuple[np.ndarray, bool, Optional[List[str]]]:
        """
        Get color values for all cells.

        Parameters
        ----------
        color : str
            Column in obs or gene name
        vmin, vmax : float, optional
            Min/max for continuous data

        Returns
        -------
        values : np.ndarray
            Numeric values for each cell
        is_continuous : bool
            Whether data is continuous
        categories : list or None
            Category names if categorical, else None
        """
        if color in self.adata.obs.columns:
            col = self.adata.obs[color]
            if pd.api.types.is_categorical_dtype(col):
                categories = list(col.cat.categories)
                values = col.cat.codes.to_numpy().astype(float)
                # Handle NaN codes (-1)
                values[values < 0] = np.nan
                return values, False, categories
            elif pd.api.types.is_numeric_dtype(col):
                values = col.to_numpy(dtype=float)
                return values, True, None
            else:
                # Convert to categorical
                cat = col.astype("category")
                categories = list(cat.cat.categories)
                values = cat.cat.codes.to_numpy().astype(float)
                values[values < 0] = np.nan
                return values, False, categories
        elif color in self.adata.var_names:
            # Gene expression (prefer normalized layer when available)
            gene_idx = self.adata.var_names.get_loc(color)
            expr_layer = None
            if "normalized" in self.adata.layers:
                expr_layer = self.adata.layers["normalized"]
            x = expr_layer[:, gene_idx] if expr_layer is not None else self.adata.X[:, gene_idx]
            if issparse(x):
                values = np.asarray(x.toarray()).ravel()
            else:
                values = np.asarray(x).ravel()
            return values, True, None
        else:
            raise KeyError(f"{color!r} not found in obs columns or var_names")

    def get_section_indices(self) -> Dict[str, np.ndarray]:
        """Get cell indices for each section."""
        indices = {}
        gvals = self.adata.obs[self.groupby].astype(str).to_numpy()
        for section in self.sections:
            indices[section.section_id] = np.flatnonzero(gvals == section.section_id)
        return indices

    def get_metadata_filters(self) -> Dict[str, List[str]]:
        """Get unique values for filterable metadata columns."""
        filters = {}
        for col in self.metadata_columns:
            if col in self.adata.obs.columns:
                unique_vals = list(self.adata.obs[col].dropna().astype(str).unique())
                custom_order = None
                if self.metadata_value_order and col in self.metadata_value_order:
                    custom_order = [str(v) for v in self.metadata_value_order[col]]
                if custom_order:
                    custom_set = set(custom_order)
                    ordered = [v for v in custom_order if v in unique_vals]
                    remaining = [v for v in unique_vals if v not in custom_set]
                    if col == "last_day":
                        def _sort_key(v):
                            try:
                                return (0, float(v))
                            except ValueError:
                                return (1, v)
                        remaining = sorted(remaining, key=_sort_key)
                    else:
                        remaining = sorted(remaining)
                    filters[col] = ordered + remaining
                elif col == "last_day":
                    def _sort_key(v):
                        try:
                            return (0, float(v))
                        except ValueError:
                            return (1, v)
                    filters[col] = sorted(unique_vals, key=_sort_key)
                else:
                    filters[col] = sorted(unique_vals)
        return filters

    def to_json_data(
        self,
        color: str,
        downsample: Optional[int] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        additional_colors: Optional[List[str]] = None,
        genes: Optional[List[str]] = None,
        marker_genes_groupby: Optional[List[str]] = None,
        marker_genes_top_n: int = 30,
        neighbor_stats_groupby: Optional[List[str]] = None,
        neighbor_stats_permutations: int = 0,
        neighbor_stats_seed: int = 0,
    ) -> Dict:
        """
        Export dataset to JSON-serializable format for the HTML viewer.

        Parameters
        ----------
        color : str
            Initial color column or gene
        downsample : int, optional
            If set, randomly downsample to this many cells per section
        vmin, vmax : float, optional
            Min/max for continuous color scale
        additional_colors : list, optional
            Additional obs columns to include for color switching
        genes : list, optional
            Gene names to include for expression visualization
        marker_genes_groupby : list, optional
            Obs columns to compute marker genes for (categorical only)
        marker_genes_top_n : int
            Number of top marker genes to keep per group
        neighbor_stats_groupby : list, optional
            Obs columns to compute neighbor composition stats for (categorical only)
        neighbor_stats_permutations : int
            Number of permutations for neighbor enrichment z-scores (0 disables)
        neighbor_stats_seed : int
            Random seed used for neighbor permutations

        Returns
        -------
        dict
            JSON-serializable data structure
        """
        coords = np.asarray(self.adata.obsm["spatial"])[:, :2]
        section_indices = self.get_section_indices()

        # Get UMAP coordinates if available
        umap_coords = None
        umap_bounds = None
        if self.has_umap:
            umap_coords = np.asarray(self.adata.obsm["X_umap"])[:, :2]
            # Compute global UMAP bounds for consistent scaling across all sections
            umap_bounds = {
                "xmin": float(umap_coords[:, 0].min()),
                "xmax": float(umap_coords[:, 0].max()),
                "ymin": float(umap_coords[:, 1].min()),
                "ymax": float(umap_coords[:, 1].max()),
            }

        # Get neighborhood graph if available
        neighbor_graph = None
        neighbor_graph_key = None
        for key in ("spatial_connectivities", "connectivities", "neighbors", "neighbor_graph"):
            if key in self.adata.obsp:
                neighbor_graph = self.adata.obsp[key]
                neighbor_graph_key = key
                break
        if neighbor_graph is not None:
            if not issparse(neighbor_graph):
                neighbor_graph = sp.csr_matrix(neighbor_graph)
            else:
                neighbor_graph = neighbor_graph.tocsr()

        # Get initial color data
        values, is_continuous, categories = self.get_color_data(color, vmin, vmax)

        # Compute global bounds for initial color
        if is_continuous:
            finite_mask = np.isfinite(values)
            if finite_mask.any():
                global_vmin = float(np.nanmin(values[finite_mask])) if vmin is None else vmin
                global_vmax = float(np.nanmax(values[finite_mask])) if vmax is None else vmax
            else:
                global_vmin, global_vmax = 0.0, 1.0
        else:
            global_vmin, global_vmax = None, None

        # Build list of all colors to export
        all_colors = [color]
        if additional_colors:
            all_colors.extend([c for c in additional_colors if c != color and c in self.obs_columns])

        # Pre-compute all color data
        color_data = {}
        for col in all_colors:
            try:
                vals, is_cont, cats = self.get_color_data(col)
                if is_cont:
                    finite = np.isfinite(vals)
                    col_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                    col_vmax = float(np.nanmax(vals[finite])) if finite.any() else 1.0
                else:
                    col_vmin, col_vmax = None, None
                color_data[col] = {
                    "values": vals,
                    "is_continuous": is_cont,
                    "categories": cats,
                    "vmin": col_vmin,
                    "vmax": col_vmax,
                }
            except Exception as e:
                print(f"  Warning: Could not load color '{col}': {e}")

        # Pre-compute gene expression data
        gene_data = {}
        if genes:
            for gene in genes:
                if gene in self.adata.var_names:
                    try:
                        vals, _, _ = self.get_color_data(gene)
                        finite = np.isfinite(vals)
                        gene_vmin = float(np.nanmin(vals[finite])) if finite.any() else 0.0
                        gene_vmax = float(np.nanmax(vals[finite])) if finite.any() else 1.0
                        gene_data[gene] = {
                            "values": vals,
                            "vmin": gene_vmin,
                            "vmax": gene_vmax,
                        }
                    except Exception as e:
                        print(f"  Warning: Could not load gene '{gene}': {e}")

        # Get metadata filters
        metadata_filters = self.get_metadata_filters()

        # Compute marker genes for requested groupby columns
        marker_genes = {}
        if marker_genes_groupby:
            for groupby in marker_genes_groupby:
                if groupby not in self.adata.obs.columns:
                    print(f"  Warning: marker_genes groupby '{groupby}' not found in obs.")
                    continue
                col = self.adata.obs[groupby]
                if not pd.api.types.is_categorical_dtype(col):
                    self.adata.obs[groupby] = col.astype("category")
                key_added = f"rank_genes_groups_{groupby}"
                alt_key_added = f"rank_genes_groups__{groupby}"
                existing_key = None
                if key_added in self.adata.uns:
                    existing_key = key_added
                elif alt_key_added in self.adata.uns:
                    existing_key = alt_key_added

                if existing_key is None:
                    try:
                        sc.tl.rank_genes_groups(
                            self.adata,
                            groupby=groupby,
                            reference="rest",
                            method="t-test",
                            pts=False,
                            key_added=key_added,
                        )
                    except Exception as e:
                        print(f"  Warning: Could not compute marker genes for '{groupby}': {e}")
                        continue

                rg = self.adata.uns.get(existing_key or key_added)
                if not rg:
                    print(f"  Warning: marker genes not found for '{groupby}'.")
                    continue

                names = rg.get("names")
                if names is None:
                    print(f"  Warning: marker genes missing names for '{groupby}'.")
                    continue

                if isinstance(names, pd.DataFrame):
                    marker_genes[groupby] = {
                        col_name: names[col_name].astype(str).tolist()[:marker_genes_top_n]
                        for col_name in names.columns
                    }
                elif isinstance(names, np.ndarray) and names.dtype.names:
                    marker_genes[groupby] = {
                        group: [str(x) for x in names[group][:marker_genes_top_n]]
                        for group in names.dtype.names
                    }
                else:
                    print(f"  Warning: Unrecognized marker gene format for '{groupby}'.")

        # Compute neighbor composition stats
        neighbor_stats = {}
        if neighbor_graph is not None and neighbor_stats_groupby:
            for groupby in neighbor_stats_groupby:
                if groupby not in self.adata.obs.columns:
                    print(f"  Warning: neighbor stats groupby '{groupby}' not found in obs.")
                    continue
                col = self.adata.obs[groupby]
                if pd.api.types.is_numeric_dtype(col):
                    print(f"  Warning: neighbor stats '{groupby}' is numeric; skipping.")
                    continue
                if not pd.api.types.is_categorical_dtype(col):
                    col = col.astype("category")
                categories = list(col.cat.categories)
                codes = col.cat.codes.to_numpy()
                valid_mask = codes >= 0
                if not valid_mask.any():
                    print(f"  Warning: neighbor stats '{groupby}' has no valid categories.")
                    continue

                if valid_mask.all():
                    graph = neighbor_graph
                    labels = codes
                else:
                    valid_idx = np.flatnonzero(valid_mask)
                    graph = neighbor_graph[valid_idx][:, valid_idx]
                    labels = codes[valid_mask]

                n_cells = np.bincount(labels, minlength=len(categories)).astype(int)
                if graph is None or graph.shape[0] == 0:
                    print(f"  Warning: neighbor stats '{groupby}' has empty graph.")
                    continue

                onehot = sp.csr_matrix(
                    (np.ones(len(labels), dtype=float), (np.arange(len(labels)), labels)),
                    shape=(len(labels), len(categories)),
                )
                counts = onehot.T.dot(graph).dot(onehot)
                if issparse(counts):
                    counts = counts.toarray()
                counts = np.asarray(counts, dtype=float)
                row_sums = counts.sum(axis=1)
                mean_degree = np.zeros(len(categories), dtype=float)
                valid_cells = n_cells > 0
                mean_degree[valid_cells] = row_sums[valid_cells] / n_cells[valid_cells]

                entry = {
                    "categories": categories,
                    "counts": counts.tolist(),
                    "n_cells": n_cells.tolist(),
                    "mean_degree": mean_degree.tolist(),
                }
                if neighbor_stats_permutations and neighbor_stats_permutations > 0:
                    rng = np.random.default_rng(int(neighbor_stats_seed))
                    perm_mean = np.zeros_like(counts, dtype=float)
                    perm_m2 = np.zeros_like(counts, dtype=float)
                    for i in range(int(neighbor_stats_permutations)):
                        perm_labels = rng.permutation(labels)
                        perm_onehot = sp.csr_matrix(
                            (np.ones(len(perm_labels), dtype=float), (np.arange(len(perm_labels)), perm_labels)),
                            shape=(len(perm_labels), len(categories)),
                        )
                        perm_counts = perm_onehot.T.dot(graph).dot(perm_onehot)
                        if issparse(perm_counts):
                            perm_counts = perm_counts.toarray()
                        perm_counts = np.asarray(perm_counts, dtype=float)
                        delta = perm_counts - perm_mean
                        perm_mean += delta / (i + 1)
                        perm_m2 += delta * (perm_counts - perm_mean)
                    if neighbor_stats_permutations > 1:
                        perm_var = perm_m2 / (neighbor_stats_permutations - 1)
                    else:
                        perm_var = np.zeros_like(counts, dtype=float)
                    perm_std = np.sqrt(perm_var)
                    zscore = np.zeros_like(counts, dtype=float)
                    valid_std = perm_std > 0
                    zscore[valid_std] = (counts[valid_std] - perm_mean[valid_std]) / perm_std[valid_std]
                    entry["perm_n"] = int(neighbor_stats_permutations)
                    entry["zscore"] = zscore.tolist()
                neighbor_stats[groupby] = entry

        # Build section data with all color layers
        sections_data = []
        for section in self.sections:
            idx = section_indices[section.section_id]

            if downsample and len(idx) > downsample:
                rng = np.random.default_rng(42)
                idx = rng.choice(idx, size=downsample, replace=False)
                idx = np.sort(idx)

            section_coords = coords[idx]

            # Get UMAP coordinates for this section if available
            section_umap = None
            if umap_coords is not None:
                section_umap = umap_coords[idx]

            # Build color values for this section
            section_colors = {}
            for col, cdata in color_data.items():
                section_vals = cdata["values"][idx]
                # Convert numpy types to native Python types for JSON serialization
                section_colors[col] = [
                    float(v) if np.isfinite(v) else None for v in section_vals
                ]

            # Build gene expression values for this section
            section_genes = {}
            for gene, gdata in gene_data.items():
                section_vals = gdata["values"][idx]
                section_genes[gene] = [
                    float(v) if np.isfinite(v) else None for v in section_vals
                ]

            section_entry = {
                "id": section.section_id,
                "metadata": section.metadata,
                "n_cells": int(len(idx)),
                "x": section_coords[:, 0].tolist(),
                "y": section_coords[:, 1].tolist(),
                "colors": section_colors,
                "genes": section_genes,
                "bounds": {
                    "xmin": float(section_coords[:, 0].min()) if len(idx) > 0 else 0,
                    "xmax": float(section_coords[:, 0].max()) if len(idx) > 0 else 0,
                    "ymin": float(section_coords[:, 1].min()) if len(idx) > 0 else 0,
                    "ymax": float(section_coords[:, 1].max()) if len(idx) > 0 else 0,
                }
            }

            # Add UMAP coordinates if available
            if section_umap is not None:
                section_entry["umap_x"] = section_umap[:, 0].tolist()
                section_entry["umap_y"] = section_umap[:, 1].tolist()

            if neighbor_graph is not None:
                subgraph = neighbor_graph[idx][:, idx]
                if issparse(subgraph) and subgraph.nnz > 0:
                    upper = sp.triu(subgraph, k=1).tocoo()
                    section_entry["edges"] = list(zip(upper.row.tolist(), upper.col.tolist()))
                else:
                    section_entry["edges"] = []

            sections_data.append(section_entry)

        # Build color metadata
        colors_meta = {}
        for col, cdata in color_data.items():
            colors_meta[col] = {
                "is_continuous": cdata["is_continuous"],
                "categories": cdata["categories"],
                "vmin": cdata["vmin"],
                "vmax": cdata["vmax"],
            }

        # Build gene metadata
        genes_meta = {}
        for gene, gdata in gene_data.items():
            genes_meta[gene] = {
                "vmin": gdata["vmin"],
                "vmax": gdata["vmax"],
            }

        return {
            "initial_color": color,
            "colors_meta": colors_meta,
            "genes_meta": genes_meta,
            "metadata_filters": metadata_filters,
            "n_sections": len(sections_data),
            "total_cells": sum(s["n_cells"] for s in sections_data),
            "sections": sections_data,
            "available_colors": list(color_data.keys()),
            "available_genes": list(gene_data.keys()),
            "marker_genes": marker_genes,
            "has_umap": umap_coords is not None,
            "umap_bounds": umap_bounds,
            "has_neighbors": neighbor_graph is not None,
            "neighbors_key": neighbor_graph_key,
            "neighbor_stats": neighbor_stats,
        }


def load_spatial_data(
    path: str,
    groupby: str = "sample_id",
    spatial_key: str = "spatial",
    group_order: Optional[List[str]] = None,
    metadata_columns: Optional[List[str]] = None,
    metadata_value_order: Optional[Dict[str, List[str]]] = None,
    metadata_max_columns: Optional[int] = None,
) -> SpatialDataset:
    """
    Load spatial transcriptomics data from h5ad file.

    Parameters
    ----------
    path : str
        Path to .h5ad file
    groupby : str
        Column in obs to group sections by
    spatial_key : str
        Key in obsm containing spatial coordinates
    group_order : list, optional
        Custom order for sections
    metadata_columns : list, optional
        Obs columns to use for section metadata and filter chips
    metadata_value_order : dict, optional
        Custom ordering for metadata values per column (e.g. {"course": ["A", "B"]})
        If group_order is not provided, the first key in this dict is used to order sections
        by that metadata column (unknowns last, then section_id sort).
    metadata_max_columns : int, optional
        Limit the number of metadata columns used (order preserved)

    Returns
    -------
    SpatialDataset
        Loaded dataset ready for visualization
    """
    print(f"Loading {path}...")
    adata = sc.read_h5ad(path)
    print(f"  Loaded {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")

    if groupby not in adata.obs.columns:
        raise ValueError(f"Groupby column '{groupby}' not found in adata.obs")

    # Determine section order
    gser = adata.obs[groupby]
    gser_str = gser.astype(str)
    if group_order is not None:
        section_ids = [str(g) for g in group_order if str(g) in gser_str.unique()]
    else:
        order_by_meta = None
        if metadata_value_order:
            order_by_meta = next(iter(metadata_value_order.keys()), None)
        if order_by_meta and order_by_meta in adata.obs.columns:
            desired_order = [str(v) for v in metadata_value_order.get(order_by_meta, [])]
            desired_index = {v: i for i, v in enumerate(desired_order)}
            section_ids = []
            for sid in gser_str.unique():
                mask = gser_str == str(sid)
                vals = adata.obs.loc[mask, order_by_meta].dropna().astype(str).unique()
                meta_value = vals[0] if len(vals) == 1 else "mixed"
                section_ids.append((str(sid), meta_value))
            def _order_key(item):
                sid, meta_value = item
                if meta_value in desired_index:
                    return (0, desired_index[meta_value], sid)
                return (1, meta_value, sid)
            section_ids = [sid for sid, _ in sorted(section_ids, key=_order_key)]
        elif pd.api.types.is_categorical_dtype(gser) and gser.cat.ordered:
            section_ids = [str(c) for c in gser.cat.categories if str(c) in gser_str.unique()]
        else:
            section_ids = sorted(gser_str.unique())

    print(f"  Found {len(section_ids)} sections")

    # Determine metadata columns
    if metadata_columns is None:
        metadata_columns = ["course", "region", "condition", "timepoint", "last_score", "last_day"]
    if metadata_max_columns is not None:
        if metadata_max_columns < 0:
            raise ValueError("metadata_max_columns must be >= 0")
        metadata_columns = metadata_columns[:metadata_max_columns]

    # Build section data
    coords = np.asarray(adata.obsm[spatial_key])[:, :2]
    gvals = gser.astype(str).to_numpy()

    sections = []
    for sid in section_ids:
        mask = gvals == sid
        section_coords = coords[mask]

        # Extract metadata
        metadata = {}
        for meta_col in metadata_columns:
            if meta_col in adata.obs.columns:
                vals = adata.obs.loc[mask, meta_col].dropna().astype(str).unique()
                if len(vals) == 1:
                    metadata[meta_col] = vals[0]
                elif len(vals) > 1:
                    metadata[meta_col] = "mixed"

        sections.append(SectionData(
            section_id=sid,
            coordinates=section_coords,
            metadata=metadata,
        ))

    # Get available columns for coloring
    obs_columns = [
        col for col in adata.obs.columns
        if pd.api.types.is_categorical_dtype(adata.obs[col])
        or pd.api.types.is_numeric_dtype(adata.obs[col])
    ]

    return SpatialDataset(
        adata=adata,
        sections=sections,
        groupby=groupby,
        obs_columns=obs_columns,
        var_names=list(adata.var_names),
        metadata_columns=metadata_columns,
        metadata_value_order=metadata_value_order,
    )

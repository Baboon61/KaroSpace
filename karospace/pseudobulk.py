"""Pseudobulk differential expression utilities."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from pandas.api.types import CategoricalDtype


def _as_count_matrix(adata, counts_layer: Optional[str]) -> Tuple[Any, str, Optional[str]]:
    if counts_layer:
        if counts_layer in adata.layers:
            return adata.layers[counts_layer], str(counts_layer), None
        return adata.X, "X", f"counts layer '{counts_layer}' not found; using adata.X"
    return adata.X, "X", None


def _to_dense_counts(matrix) -> np.ndarray:
    dense = matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)
    dense = np.asarray(dense, dtype=np.float64)
    dense[~np.isfinite(dense)] = 0
    dense[dense < 0] = 0
    return np.rint(dense).astype(np.int64, copy=False)


def _positive_fraction(matrix, mask: np.ndarray, gene_indices: Sequence[int]) -> List[Optional[float]]:
    if not mask.any():
        return [None for _ in gene_indices]
    subset = matrix[mask]
    if sp.issparse(subset):
        subset = subset[:, list(gene_indices)]
        counts = np.asarray((subset > 0).sum(axis=0)).ravel()
    else:
        subset = np.asarray(subset)[:, list(gene_indices)]
        counts = np.count_nonzero(subset > 0, axis=0)
    denom = int(mask.sum())
    return [float(v) / denom for v in counts]


def _adjust_pvalues(pvalues: np.ndarray, method: str = "fdr_bh") -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    adjusted = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return adjusted
    idx = np.flatnonzero(finite)
    vals = np.clip(p[idx], 0, 1)
    method_norm = str(method or "fdr_bh").strip().lower().replace("-", "_")
    if method_norm in {"none", "raw", "pvalue", "pvalues"}:
        adjusted[idx] = vals
        adjusted[(adjusted == 0) & np.isfinite(adjusted)] = np.nextafter(0.0, 1.0)
        return adjusted
    if method_norm in {"bonferroni", "bonf"}:
        adjusted[idx] = np.clip(vals * len(vals), 0, 1)
        adjusted[(adjusted == 0) & np.isfinite(adjusted)] = np.nextafter(0.0, 1.0)
        return adjusted

    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    if method_norm in {"holm", "holm_bonferroni"}:
        adj = (n - np.arange(n, dtype=float)) * ranked
        adj = np.maximum.accumulate(adj)
    else:
        adj = ranked * n / (np.arange(n, dtype=float) + 1.0)
        adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    adjusted[idx[order]] = adj
    adjusted[(adjusted == 0) & np.isfinite(adjusted)] = np.nextafter(0.0, 1.0)
    return adjusted


def _normalize_pct_threshold(value: float) -> float:
    threshold = float(value or 0.0)
    if threshold > 1.0:
        threshold = threshold / 100.0
    return min(max(threshold, 0.0), 1.0)


def _json_float(value: float, digits: int = 6) -> Optional[float]:
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(val):
        return None
    return float(round(val, digits))


def _compute_pseudobulk_sample_diagnostics(
    pair_counts: np.ndarray,
    pair_meta: pd.DataFrame,
    *,
    max_genes: int = 1000,
) -> Dict[str, Any]:
    counts = np.array(pair_counts, dtype=float, copy=True)
    counts[~np.isfinite(counts)] = 0
    counts[counts < 0] = 0
    if counts.ndim != 2 or counts.shape[0] < 2 or counts.shape[1] < 1:
        return {}

    library_sizes = counts.sum(axis=1)
    safe_libs = np.where(library_sizes > 0, library_sizes, 1.0)
    log_cpm = np.log1p((counts / safe_libs[:, None]) * 1_000_000.0)

    variances = np.var(log_cpm, axis=0)
    finite_var = np.isfinite(variances) & (variances > 0)
    if not finite_var.any():
        feature_idx = np.arange(log_cpm.shape[1])
    else:
        feature_idx = np.flatnonzero(finite_var)
    if feature_idx.size > int(max_genes):
        order = np.argsort(variances[feature_idx])[::-1][: int(max_genes)]
        feature_idx = feature_idx[order]

    matrix = log_cpm[:, feature_idx]
    matrix = matrix - np.nanmean(matrix, axis=0, keepdims=True)
    matrix[~np.isfinite(matrix)] = 0

    pca_points: List[List[Optional[float]]] = []
    pca_variance: List[Optional[float]] = [None, None]
    try:
        u, s, _vt = np.linalg.svd(matrix, full_matrices=False)
        scores = u[:, :2] * s[:2]
        if scores.shape[1] == 1:
            scores = np.column_stack([scores[:, 0], np.zeros(scores.shape[0])])
        eig = s ** 2
        total = float(eig.sum())
        if total > 0:
            pca_variance = [_json_float(eig[i] / total, 6) if i < eig.size else None for i in range(2)]
        pca_points = [[_json_float(x, 6), _json_float(y, 6)] for x, y in scores[:, :2]]
    except Exception:
        pca_points = [[None, None] for _ in range(counts.shape[0])]

    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist, squareform

        condensed = pdist(matrix, metric="euclidean")
        dist = squareform(condensed)
        if condensed.size:
            linkage_matrix = linkage(condensed, method="average")
            order = leaves_list(linkage_matrix).astype(int).tolist()
        else:
            order = list(range(counts.shape[0]))
    except Exception:
        diff = matrix[:, None, :] - matrix[None, :, :]
        dist = np.sqrt(np.sum(diff * diff, axis=2))
        order = list(range(counts.shape[0]))

    ordered_dist = dist[np.ix_(order, order)]
    replicate_values = pair_meta["_pb_replicate"].astype(str).tolist()
    group_values = pair_meta["_pb_group"].astype(str).tolist()
    labels = [f"{replicate_values[i]} | {group_values[i]}" for i in range(len(replicate_values))]
    n_cells = pair_meta["n_cells"].to_numpy(dtype=float) if "n_cells" in pair_meta.columns else np.zeros(len(labels))

    return {
        "labels": labels,
        "replicates": replicate_values,
        "groups": group_values,
        "n_cells": [int(v) if np.isfinite(v) else None for v in n_cells],
        "library_size": [int(v) if np.isfinite(v) else None for v in library_sizes],
        "pca": pca_points,
        "pca_variance": pca_variance,
        "distance_order": [int(i) for i in order],
        "distance_labels": [labels[i] for i in order],
        "distance_groups": [group_values[i] for i in order],
        "distance_matrix": [
            [_json_float(v, 5) for v in row]
            for row in ordered_dist
        ],
        "distance_metric": "euclidean_log1p_cpm",
        "pca_features": int(feature_idx.size),
    }


def _compute_category_gene_means_from_aggregate(
    aggregate: np.ndarray,
    pb_meta: pd.DataFrame,
    categories: Sequence[str],
    gene_names: Sequence[str],
) -> Dict[str, Any]:
    """Summarize category-level means from replicate-level pseudobulk counts."""
    genes = [str(g) for g in gene_names]
    category_means: Dict[str, List[Optional[float]]] = {}
    category_cells: Dict[str, int] = {}
    aggregate = np.asarray(aggregate, dtype=float)
    n_cells = pb_meta["n_cells"].to_numpy(dtype=float) if "n_cells" in pb_meta else np.zeros(aggregate.shape[0])
    n_cells[~np.isfinite(n_cells)] = 0
    n_cells[n_cells < 0] = 0
    valid_rows = n_cells > 0
    per_sample_means = np.zeros_like(aggregate, dtype=float)
    if aggregate.size:
        np.divide(
            aggregate,
            n_cells[:, None],
            out=per_sample_means,
            where=valid_rows[:, None],
        )
    group_values = pb_meta["_pb_group"].astype(str).to_numpy()

    for category in [str(c) for c in categories]:
        mask = (group_values == category) & valid_rows
        cells = float(np.sum(n_cells[mask])) if mask.any() else 0.0
        category_cells[category] = int(cells)
        if mask.any():
            means = per_sample_means[mask].mean(axis=0)
            category_means[category] = [_json_float(v, 6) for v in means]
        else:
            category_means[category] = [0.0 for _ in genes]

    if "_pb_replicate" in pb_meta:
        replicate_values = pb_meta["_pb_replicate"].astype(str).to_numpy()
        replicate_means = []
        for replicate in sorted(set(replicate_values)):
            mask = (replicate_values == replicate) & valid_rows
            cells = float(np.sum(n_cells[mask])) if mask.any() else 0.0
            if cells > 0:
                replicate_means.append(np.asarray(aggregate[mask].sum(axis=0), dtype=float) / cells)
        if replicate_means:
            background = [_json_float(v, 6) for v in np.vstack(replicate_means).mean(axis=0)]
        else:
            background = [0.0 for _ in genes]
    elif valid_rows.any():
        background = [_json_float(v, 6) for v in per_sample_means[valid_rows].mean(axis=0)]
    else:
        background = [0.0 for _ in genes]

    return {
        "genes": genes,
        "categories": [str(c) for c in categories],
        "means": category_means,
        "background": background,
        "n_cells": category_cells,
        "source": "pseudobulk_aggregate",
    }


def _log2fc_from_pseudobulk_means(
    counts: np.ndarray,
    metadata: pd.DataFrame,
    source: str,
    reference: str,
) -> np.ndarray:
    """Compute source/reference log2FC from replicate-level aggregate per-cell means."""
    dense = np.asarray(counts, dtype=float)
    dense[~np.isfinite(dense)] = 0
    dense[dense < 0] = 0
    groups = metadata["_pb_group"].astype(str).to_numpy()
    if "n_cells" in metadata.columns:
        n_cells = metadata["n_cells"].to_numpy(dtype=float)
        n_cells[~np.isfinite(n_cells)] = 0
        n_cells[n_cells < 0] = 0
    else:
        n_cells = np.ones(dense.shape[0], dtype=float)

    source_mask = groups == str(source)
    reference_mask = groups == str(reference)
    n_vars = int(dense.shape[1])
    valid_rows = n_cells > 0
    per_sample_means = np.zeros_like(dense, dtype=float)
    if dense.size:
        np.divide(
            dense,
            n_cells[:, None],
            out=per_sample_means,
            where=valid_rows[:, None],
        )
    source_rows = source_mask & valid_rows
    reference_rows = reference_mask & valid_rows
    source_mean = (
        per_sample_means[source_rows].mean(axis=0)
        if source_rows.any()
        else np.zeros(n_vars, dtype=float)
    )
    reference_mean = (
        per_sample_means[reference_rows].mean(axis=0)
        if reference_rows.any()
        else np.zeros(n_vars, dtype=float)
    )
    tiny = np.nextafter(0.0, 1.0)
    log2fc = np.log2(np.maximum(source_mean, 0.0) + tiny) - np.log2(
        np.maximum(reference_mean, 0.0) + tiny
    )
    log2fc[~np.isfinite(log2fc)] = 0.0
    return np.asarray(log2fc, dtype=float)


def _fit_deseq2_pair(
    counts: np.ndarray,
    metadata: pd.DataFrame,
    source: str,
    reference: str,
    fit_type: str = "parametric",
) -> pd.DataFrame:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as exc:  # pragma: no cover - depends on optional runtime import
        raise RuntimeError(
            "PyDESeq2 is required for pseudobulk cluster DE. Install KaroSpace with "
            "pydeseq2 support or run `pip install pydeseq2`."
        ) from exc

    gene_names = [str(g) for g in metadata.attrs["gene_names"]]
    counts_df = pd.DataFrame(counts, index=metadata.index, columns=gene_names)
    meta = metadata[["_pb_replicate", "_pb_group"]].copy()
    meta["_pb_replicate"] = pd.Categorical(meta["_pb_replicate"])
    meta["_pb_group"] = pd.Categorical(
        meta["_pb_group"],
        categories=[str(reference), str(source)],
    )

    try:
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=meta,
            design="~ _pb_replicate + _pb_group",
            fit_type=fit_type,
            quiet=True,
        )
    except TypeError:
        dds = DeseqDataSet(
            counts=counts_df,
            clinical=meta,
            design_factors=["_pb_replicate", "_pb_group"],
            fit_type=fit_type,
            refit_cooks=True,
            n_cpus=1,
        )
    dds.deseq2()

    try:
        stat_res = DeseqStats(
            dds,
            contrast=["_pb_group", str(source), str(reference)],
            quiet=True,
        )
    except TypeError:
        stat_res = DeseqStats(
            dds,
            contrast=["_pb_group", str(source), str(reference)],
            n_cpus=1,
        )
    stat_res.summary()
    result = stat_res.results_df.copy()
    gene_names = [str(g) for g in metadata.attrs["gene_names"]]
    mean_log2fc = _log2fc_from_pseudobulk_means(counts, metadata, source, reference)
    if len(mean_log2fc) == len(gene_names):
        result["log2FoldChange"] = pd.Series(mean_log2fc, index=gene_names).reindex(result.index)
    return result


def _matrix_mean_var(matrix, mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    n = int(np.count_nonzero(mask))
    n_vars = int(matrix.shape[1])
    if n <= 0:
        return np.zeros(n_vars, dtype=float), np.zeros(n_vars, dtype=float), n
    subset = matrix[mask]
    if sp.issparse(subset):
        mean = np.asarray(subset.mean(axis=0)).ravel().astype(float, copy=False)
        mean_sq = np.asarray(subset.power(2).mean(axis=0)).ravel().astype(float, copy=False)
    else:
        arr = np.asarray(subset, dtype=float)
        arr[~np.isfinite(arr)] = 0
        mean = arr.mean(axis=0)
        mean_sq = np.square(arr).mean(axis=0)
    var = mean_sq - np.square(mean)
    var[~np.isfinite(var)] = 0
    var[var < 0] = 0
    if n > 1:
        var = var * (float(n) / float(n - 1))
    else:
        var = np.zeros_like(var)
    return mean, var, n


def _fit_welch_pair(
    expression_matrix,
    source_mask: np.ndarray,
    reference_mask: np.ndarray,
    var_names: Sequence[str],
) -> pd.DataFrame:
    source_mean, source_var, n_source = _matrix_mean_var(expression_matrix, source_mask)
    reference_mean, reference_var, n_reference = _matrix_mean_var(expression_matrix, reference_mask)
    denom_sq = (source_var / max(n_source, 1)) + (reference_var / max(n_reference, 1))
    denom = np.sqrt(denom_sq)
    stat = np.zeros_like(source_mean, dtype=float)
    valid = np.isfinite(denom) & (denom > 0)
    stat[valid] = (source_mean[valid] - reference_mean[valid]) / denom[valid]

    df = np.ones_like(source_mean, dtype=float)
    if n_source > 1 and n_reference > 1:
        a = source_var / n_source
        b = reference_var / n_reference
        denom_df = (np.square(a) / (n_source - 1)) + (np.square(b) / (n_reference - 1))
        valid_df = denom_df > 0
        df[valid_df] = np.square(a[valid_df] + b[valid_df]) / denom_df[valid_df]

    pvalue = np.ones_like(source_mean, dtype=float)
    valid_p = valid & np.isfinite(df) & (df > 0)
    pvalue[valid_p] = 2.0 * stats.t.sf(np.abs(stat[valid_p]), df[valid_p])
    pvalue[~np.isfinite(pvalue)] = 1.0

    tiny = np.nextafter(0.0, 1.0)
    log2fc = np.log2(np.maximum(source_mean, 0.0) + tiny) - np.log2(
        np.maximum(reference_mean, 0.0) + tiny
    )
    log2fc[~np.isfinite(log2fc)] = 0.0
    base_mean = (source_mean + reference_mean) / 2.0

    return pd.DataFrame(
        {
            "baseMean": base_mean,
            "log2FoldChange": log2fc,
            "stat": stat,
            "pvalue": np.clip(pvalue, 0, 1),
        },
        index=[str(g) for g in var_names],
    )


def _empty_result(
    reason: str,
    *,
    n_source: int,
    n_reference: int,
    min_cells: int,
    min_replicates: int,
    details: Optional[str] = None,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "available": False,
        "reason": reason,
        "method": "pseudobulk-deseq2",
        "genes": [],
        "logfoldchanges": [],
        "pvals": [],
        "pvals_adj": [],
        "scores": [],
        "pct_source": [],
        "pct_reference": [],
        "base_mean": [],
        "n_source": int(n_source),
        "n_reference": int(n_reference),
        "min_cells_required": int(min_cells),
        "min_replicates_required": int(min_replicates),
    }
    if details:
        result["details"] = details
    return result


def _empty_interaction_result(
    reason: str,
    *,
    n_contact: int,
    n_non_contact: int,
    min_cells: int,
    min_replicates: int,
    details: Optional[str] = None,
    pct_contact: float = 0.0,
    mean_target_neighbors_contact: float = 0.0,
    mean_target_neighbors_non_contact: float = 0.0,
    target_edge_count: float = 0.0,
    target_zscore: Optional[float] = None,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "available": False,
        "reason": reason,
        "method": "pseudobulk-deseq2-contact",
        "genes": [],
        "logfoldchanges": [],
        "log2foldchanges": [],
        "pvals": [],
        "pvals_adj": [],
        "scores": [],
        "pct_source": [],
        "pct_reference": [],
        "base_mean": [],
        "n_contact": int(n_contact),
        "n_non_contact": int(n_non_contact),
        "min_cells_required": int(min_cells),
        "min_replicates_required": int(min_replicates),
        "pct_contact": float(pct_contact),
        "mean_target_neighbors_contact": float(mean_target_neighbors_contact),
        "mean_target_neighbors_non_contact": float(mean_target_neighbors_non_contact),
        "target_edge_count": float(target_edge_count),
        "target_zscore": target_zscore,
    }
    if details:
        result["details"] = details
    return result


def _format_result(
    df: pd.DataFrame,
    *,
    source_mask: np.ndarray,
    reference_mask: np.ndarray,
    expression_matrix,
    var_names: Sequence[str],
    top_n: int,
    n_source: int,
    n_reference: int,
    n_replicates: int,
    counts_layer_used: str,
    warning: Optional[str],
    p_adjust_method: str,
    min_pct_expressed: float,
    padj_cutoff: float,
    log2fc_cutoff: float,
    sample_diagnostics: Optional[Dict[str, Any]] = None,
    method: str = "pseudobulk-deseq2",
) -> Dict[str, Any]:
    p_adjust_method = str(p_adjust_method or "fdr_bh").strip().lower().replace("-", "_")
    min_pct_expressed = _normalize_pct_threshold(min_pct_expressed)
    padj_cutoff = min(max(float(padj_cutoff), 0.0), 1.0)
    log2fc_cutoff = max(float(log2fc_cutoff), 0.0)
    base_payload = {
        "method": str(method or "pseudobulk-deseq2"),
        "p_adjust_method": p_adjust_method,
        "min_pct_expressed": min_pct_expressed,
        "padj_cutoff": padj_cutoff,
        "log2fc_cutoff": log2fc_cutoff,
        "table_top_n": int(top_n),
        **({"pseudobulk_samples": sample_diagnostics} if sample_diagnostics else {}),
        "n_source": int(n_source),
        "n_reference": int(n_reference),
        "n_replicates": int(n_replicates),
        "counts_layer": counts_layer_used,
        **({"warning": warning} if warning else {}),
    }
    if df is None or df.empty:
        result = _empty_result("no_results", n_source=n_source, n_reference=n_reference, min_cells=0, min_replicates=0)
        result.update(base_payload)
        return result

    work = df.copy()
    if "log2FoldChange" not in work.columns:
        work["log2FoldChange"] = np.nan
    if "pvalue" not in work.columns:
        work["pvalue"] = np.nan
    work["padj"] = _adjust_pvalues(work["pvalue"].to_numpy(dtype=float), method=p_adjust_method)
    if "stat" not in work.columns:
        work["stat"] = np.nan
    if "baseMean" not in work.columns:
        work["baseMean"] = np.nan

    work["_gene"] = [str(idx) for idx in work.index]
    work = work[np.isfinite(work["log2FoldChange"].to_numpy(dtype=float))]
    if work.empty:
        return {
            "available": True,
            "genes": [],
            "log2foldchanges": [],
            "logfoldchanges": [],
            "pvals": [],
            "pvals_adj": [],
            "scores": [],
            "pct_source": [],
            "pct_reference": [],
            "base_mean": [],
            **base_payload,
        }

    gene_to_idx = {str(g): i for i, g in enumerate(var_names)}
    gene_indices = [gene_to_idx.get(g) for g in work["_gene"]]
    valid_positions = [i for i, idx in enumerate(gene_indices) if idx is not None]
    if len(valid_positions) != len(gene_indices):
        work = work.iloc[valid_positions]
        gene_indices = [gene_indices[i] for i in valid_positions]

    pct_source = _positive_fraction(expression_matrix, source_mask, gene_indices)
    pct_reference = _positive_fraction(expression_matrix, reference_mask, gene_indices)
    work["_pct_source"] = pct_source
    work["_pct_reference"] = pct_reference
    if min_pct_expressed > 0:
        pct_source_arr = np.asarray([v if v is not None else 0.0 for v in pct_source], dtype=float)
        pct_reference_arr = np.asarray([v if v is not None else 0.0 for v in pct_reference], dtype=float)
        keep_pct = (pct_source_arr >= min_pct_expressed) | (pct_reference_arr >= min_pct_expressed)
        work = work.iloc[np.flatnonzero(keep_pct)].copy()
        pct_source = [pct_source[i] for i in np.flatnonzero(keep_pct)]
        pct_reference = [pct_reference[i] for i in np.flatnonzero(keep_pct)]

    work["_padj_sort"] = work["padj"].fillna(np.inf)
    work["_pvalue_sort"] = work["pvalue"].fillna(np.inf)
    work["_abs_lfc"] = np.abs(work["log2FoldChange"].to_numpy(dtype=float))
    work = work.sort_values(
        ["_padj_sort", "_pvalue_sort", "_abs_lfc", "_gene"],
        ascending=[True, True, False, True],
    )
    pct_source = [work["_pct_source"].iloc[i] for i in range(len(work))]
    pct_reference = [work["_pct_reference"].iloc[i] for i in range(len(work))]

    def _finite_or_none(value):
        try:
            val = float(value)
        except (TypeError, ValueError):
            return None
        return val if np.isfinite(val) else None

    result = {
        "available": True,
        "genes": work["_gene"].astype(str).tolist(),
        "log2foldchanges": [_finite_or_none(v) for v in work["log2FoldChange"]],
        "logfoldchanges": [_finite_or_none(v) for v in work["log2FoldChange"]],
        "pvals": [_finite_or_none(v) for v in work["pvalue"]],
        "pvals_adj": [_finite_or_none(v) for v in work["padj"]],
        "scores": [_finite_or_none(v) for v in work["stat"]],
        "pct_source": pct_source,
        "pct_reference": pct_reference,
        "base_mean": [_finite_or_none(v) for v in work["baseMean"]],
        **base_payload,
    }
    return result


def _truncate_gene_result(result: Dict[str, Any], top_n: int) -> Dict[str, Any]:
    if int(top_n) < 1:
        return result
    n = min(int(top_n), len(result.get("genes") or []))
    fields = [
        "genes",
        "log2foldchanges",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "scores",
        "pct_source",
        "pct_reference",
        "base_mean",
    ]
    for field in fields:
        if isinstance(result.get(field), list):
            result[field] = result[field][:n]
    result["table_top_n"] = int(top_n)
    return result


def _target_zscore_value(zscore: Optional[np.ndarray], source_idx: int, target_idx: int) -> Optional[float]:
    if not isinstance(zscore, np.ndarray):
        return None
    if source_idx >= zscore.shape[0] or target_idx >= zscore.shape[1]:
        return None
    value = float(zscore[source_idx, target_idx])
    return value if np.isfinite(value) else None


def _interaction_meta(
    *,
    target_neighbor_counts: np.ndarray,
    pos_mask: np.ndarray,
    neg_mask: np.ndarray,
    n_contact: int,
    n_non_contact: int,
    edge_count: float,
    target_zscore: Optional[float],
) -> Dict[str, Any]:
    return {
        "pct_contact": float((100.0 * n_contact) / max(1, n_contact + n_non_contact)),
        "mean_target_neighbors_contact": float(
            np.mean(target_neighbor_counts[pos_mask]) if n_contact > 0 else 0.0
        ),
        "mean_target_neighbors_non_contact": float(
            np.mean(target_neighbor_counts[neg_mask]) if n_non_contact > 0 else 0.0
        ),
        "target_edge_count": float(edge_count),
        "target_zscore": target_zscore,
    }


def compute_pseudobulk_interaction_markers(
    adata,
    annotation: str,
    *,
    replicate: str,
    graph,
    obs_idx: Sequence[int],
    labels: Sequence[int],
    categories: Sequence[str],
    neighbor_counts: np.ndarray,
    neighbor_zscore: Optional[np.ndarray] = None,
    neighbor_n_cells: Optional[Sequence[int]] = None,
    counts_layer: Optional[str] = "counts",
    top_targets: int = 8,
    top_genes: int = 20,
    min_cells: int = 30,
    min_neighbors: int = 1,
    min_replicates: int = 2,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.5,
    fit_type: str = "parametric",
) -> Optional[Dict[str, Dict[str, Dict[str, Any]]]]:
    """Compute replicate-aware contact-conditioned markers for one annotation.

    For each directional source -> target pair, source cells are split into
    contact+ and contact- within each replicate using the neighbor graph
    restricted to that replicate. Counts are then aggregated by
    replicate/contact status and tested with the same paired DESeq2 design used
    for category pseudobulk DE.
    """
    if replicate not in adata.obs.columns:
        print(f"  Warning: interaction pseudobulk replicate '{replicate}' not found in obs.")
        return None

    graph = graph.tocsr() if sp.issparse(graph) else sp.csr_matrix(graph)
    contact_graph = graph.copy()
    contact_graph.data = np.ones(contact_graph.nnz, dtype=np.float32)
    contact_graph.eliminate_zeros()
    obs_idx = np.asarray(obs_idx, dtype=np.int64)
    labels = np.asarray(labels, dtype=np.int32)
    categories = [str(cat) for cat in categories]
    if sp.issparse(neighbor_counts):
        counts = np.asarray(neighbor_counts.toarray(), dtype=float)
    else:
        counts = np.asarray(neighbor_counts, dtype=float)
    zscore = np.asarray(neighbor_zscore, dtype=float) if neighbor_zscore is not None else None
    if neighbor_n_cells is not None:
        n_cells = np.asarray(neighbor_n_cells, dtype=int)
    else:
        n_cells = np.bincount(labels[labels >= 0], minlength=len(categories)).astype(int)
    if graph.shape[0] != len(labels) or len(labels) != len(obs_idx):
        print(f"  Warning: interaction markers '{annotation}' graph/label size mismatch; skipping.")
        return None

    expression_matrix, counts_layer_used, warning = _as_count_matrix(adata, counts_layer)
    rep_values = adata.obs[replicate].astype(str).to_numpy()
    ctx_reps = rep_values[obs_idx]
    valid_reps = np.asarray(pd.notna(ctx_reps), dtype=bool)
    if not valid_reps.any():
        return None

    rep_to_positions: Dict[str, np.ndarray] = {
        str(rep): np.flatnonzero((ctx_reps == rep) & valid_reps)
        for rep in sorted(set(ctx_reps[valid_reps].astype(str)))
    }
    rep_subgraph_cache: Dict[str, Any] = {
        rep: contact_graph[positions][:, positions]
        for rep, positions in rep_to_positions.items()
        if positions.size > 0
    }
    top_targets = int(top_targets)
    top_genes = int(top_genes)
    min_cells = int(min_cells)
    min_neighbors = int(min_neighbors)
    min_replicates = int(min_replicates)
    print(
        f"    - {len(categories)} categories -> contact pseudobulk by {replicate} "
        f"(top {top_targets} targets/source)",
        flush=True,
    )

    results: Dict[str, Dict[str, Dict[str, Any]]] = {}
    target_neighbor_counts_cache: Dict[int, np.ndarray] = {}
    for source_idx, source_name in enumerate(categories):
        if source_idx >= len(n_cells) or int(n_cells[source_idx]) <= 0:
            continue
        source_mask = labels == source_idx
        if not source_mask.any():
            continue
        row = counts[source_idx] if source_idx < counts.shape[0] else None
        if row is None:
            continue

        candidate_targets = [
            t for t in range(len(categories))
            if t != source_idx and t < len(row) and float(row[t]) > 0
        ]
        if not candidate_targets:
            continue

        def _target_sort_key(tidx: int):
            zval = _target_zscore_value(zscore, source_idx, tidx)
            if zval is not None:
                return (0, -float(zval), -float(row[tidx]), categories[tidx])
            return (1, 0.0, -float(row[tidx]), categories[tidx])

        ranked_targets = sorted(candidate_targets, key=_target_sort_key)[:top_targets]
        source_result: Dict[str, Dict[str, Any]] = {}

        for target_idx in ranked_targets:
            target_name = categories[target_idx]
            target_neighbor_counts = target_neighbor_counts_cache.get(target_idx)
            if target_neighbor_counts is None:
                target_neighbor_counts = np.zeros(len(labels), dtype=float)
                target_mask = labels == target_idx

                # Classify contact status within each replicate using the induced
                # replicate subgraph. This prevents cross-section/mouse edges in a
                # global graph from defining contact status.
                for rep, positions in rep_to_positions.items():
                    target_vec = target_mask[positions].astype(np.float32, copy=False)
                    if not np.any(target_vec):
                        continue
                    subgraph = rep_subgraph_cache.get(rep)
                    if subgraph is None:
                        continue
                    target_neighbor_counts[positions] = np.asarray(subgraph.dot(target_vec)).ravel()
                target_neighbor_counts_cache[target_idx] = target_neighbor_counts
            target_mask = labels == target_idx

            pos_mask = source_mask & (target_neighbor_counts >= min_neighbors)
            neg_mask = source_mask & (target_neighbor_counts == 0)
            n_pos = int(np.count_nonzero(pos_mask))
            n_neg = int(np.count_nonzero(neg_mask))
            target_z = _target_zscore_value(zscore, source_idx, target_idx)
            meta = _interaction_meta(
                target_neighbor_counts=target_neighbor_counts,
                pos_mask=pos_mask,
                neg_mask=neg_mask,
                n_contact=n_pos,
                n_non_contact=n_neg,
                edge_count=float(row[target_idx]),
                target_zscore=target_z,
            )
            source_cell_mask = np.zeros(adata.n_obs, dtype=bool)
            reference_cell_mask = np.zeros(adata.n_obs, dtype=bool)
            source_cell_mask[obs_idx[pos_mask]] = True
            reference_cell_mask[obs_idx[neg_mask]] = True

            sample_keys: List[Tuple[str, str]] = []
            sample_index: Dict[Tuple[str, str], int] = {}
            selected_positions: List[int] = []
            selected_rows: List[int] = []
            for rep, positions in rep_to_positions.items():
                rep_pos = positions[pos_mask[positions]]
                rep_neg = positions[neg_mask[positions]]
                for status, status_positions in (("contact+", rep_pos), ("contact-", rep_neg)):
                    if status_positions.size < min_cells:
                        continue
                    key = (str(rep), status)
                    row_idx = sample_index.get(key)
                    if row_idx is None:
                        row_idx = len(sample_keys)
                        sample_index[key] = row_idx
                        sample_keys.append(key)
                    selected_positions.extend(status_positions.astype(int).tolist())
                    selected_rows.extend([row_idx] * int(status_positions.size))

            contact_reps = {rep for rep, status in sample_keys if status == "contact+"}
            non_contact_reps = {rep for rep, status in sample_keys if status == "contact-"}
            paired_reps = sorted(contact_reps & non_contact_reps)
            if len(rep_to_positions) == 1 and n_pos >= min_cells and n_neg >= min_cells:
                try:
                    print(
                        f"      - {source_name} -> {target_name}: fitting Welch t-test "
                        f"({n_pos} contact+ vs {n_neg} contact- cells, {adata.n_vars} genes)",
                        flush=True,
                    )
                    pair_result = _fit_welch_pair(
                        expression_matrix,
                        source_cell_mask,
                        reference_cell_mask,
                        adata.var_names,
                    )
                    formatted = _format_result(
                        pair_result,
                        source_mask=source_cell_mask,
                        reference_mask=reference_cell_mask,
                        expression_matrix=expression_matrix,
                        var_names=adata.var_names,
                        top_n=top_genes,
                        n_source=n_pos,
                        n_reference=n_neg,
                        n_replicates=1,
                        counts_layer_used=counts_layer_used,
                        warning=warning,
                        p_adjust_method=p_adjust_method,
                        min_pct_expressed=min_pct_expressed,
                        padj_cutoff=padj_cutoff,
                        log2fc_cutoff=log2fc_cutoff,
                        method="welch-t-test-contact",
                    )
                    formatted["method"] = "welch-t-test-contact"
                    formatted["n_contact"] = int(n_pos)
                    formatted["n_non_contact"] = int(n_neg)
                    formatted["min_cells_required"] = int(min_cells)
                    formatted["min_replicates_required"] = int(min_replicates)
                    formatted.update(meta)
                    source_result[target_name] = _truncate_gene_result(formatted, top_genes)
                except Exception as exc:
                    print(f"      - {source_name} -> {target_name}: failed ({exc})", flush=True)
                    source_result[target_name] = _empty_interaction_result(
                        "de_failed",
                        n_contact=n_pos,
                        n_non_contact=n_neg,
                        min_cells=min_cells,
                        min_replicates=min_replicates,
                        details=str(exc),
                        **meta,
                    )
                continue
            if len(paired_reps) < min_replicates:
                print(
                    f"      - {source_name} -> {target_name}: skipped, insufficient paired replicates "
                    f"({len(paired_reps)}; need >= {min_replicates})",
                    flush=True,
                )
                source_result[target_name] = _empty_interaction_result(
                    "insufficient_replicates",
                    n_contact=n_pos,
                    n_non_contact=n_neg,
                    min_cells=min_cells,
                    min_replicates=min_replicates,
                    details=f"{len(paired_reps)} paired replicate(s) available",
                    **meta,
                )
                continue

            kept_sample_keys = [
                key
                for key in sample_keys
                if key[0] in paired_reps and key[1] in {"contact+", "contact-"}
            ]
            keep_sample = {key: idx for idx, key in enumerate(kept_sample_keys)}
            remapped_rows = []
            kept_obs_positions = []
            for position, row_idx in zip(selected_positions, selected_rows):
                key = sample_keys[row_idx]
                new_idx = keep_sample.get(key)
                if new_idx is None:
                    continue
                kept_obs_positions.append(int(obs_idx[position]))
                remapped_rows.append(int(new_idx))

            if not kept_obs_positions:
                source_result[target_name] = _empty_interaction_result(
                    "insufficient_cells",
                    n_contact=n_pos,
                    n_non_contact=n_neg,
                    min_cells=min_cells,
                    min_replicates=min_replicates,
                    **meta,
                )
                continue

            incidence = sp.csr_matrix(
                (
                    np.ones(len(kept_obs_positions), dtype=np.float64),
                    (np.asarray(remapped_rows, dtype=np.int64), np.asarray(kept_obs_positions, dtype=np.int64)),
                ),
                shape=(len(kept_sample_keys), adata.n_obs),
            )
            pair_counts = _to_dense_counts(incidence @ expression_matrix)
            cell_counts = np.bincount(np.asarray(remapped_rows, dtype=np.int64), minlength=len(keep_sample)).astype(int)
            ordered_keys = kept_sample_keys
            pair_meta = pd.DataFrame(
                {
                    "_pb_replicate": [key[0] for key in ordered_keys],
                    "_pb_group": [key[1] for key in ordered_keys],
                    "n_cells": cell_counts,
                },
                index=[f"pb_{i}" for i in range(len(ordered_keys))],
            )
            pair_meta.attrs["gene_names"] = [str(g) for g in adata.var_names]

            try:
                print(
                    f"      - {source_name} -> {target_name}: fitting DESeq2 "
                    f"({len(paired_reps)} paired replicate"
                    f"{'s' if len(paired_reps) != 1 else ''}, "
                    f"{pair_counts.shape[0]} pseudobulk samples, "
                    f"{pair_counts.shape[1]} genes)",
                    flush=True,
                )
                sample_diagnostics = _compute_pseudobulk_sample_diagnostics(pair_counts, pair_meta)
                pair_result = _fit_deseq2_pair(
                    pair_counts,
                    pair_meta,
                    "contact+",
                    "contact-",
                    fit_type=str(fit_type or "parametric"),
                )
                formatted = _format_result(
                    pair_result,
                    source_mask=source_cell_mask,
                    reference_mask=reference_cell_mask,
                    expression_matrix=expression_matrix,
                    var_names=adata.var_names,
                    top_n=top_genes,
                    n_source=n_pos,
                    n_reference=n_neg,
                    n_replicates=len(paired_reps),
                    counts_layer_used=counts_layer_used,
                    warning=warning,
                    p_adjust_method=p_adjust_method,
                    min_pct_expressed=min_pct_expressed,
                    padj_cutoff=padj_cutoff,
                    log2fc_cutoff=log2fc_cutoff,
                    sample_diagnostics=sample_diagnostics,
                )
                formatted["method"] = "pseudobulk-deseq2-contact"
                formatted["n_contact"] = int(n_pos)
                formatted["n_non_contact"] = int(n_neg)
                formatted["min_cells_required"] = int(min_cells)
                formatted["min_replicates_required"] = int(min_replicates)
                formatted.update(meta)
                source_result[target_name] = _truncate_gene_result(formatted, top_genes)
            except Exception as exc:
                print(f"      - {source_name} -> {target_name}: failed ({exc})", flush=True)
                source_result[target_name] = _empty_interaction_result(
                    "de_failed",
                    n_contact=n_pos,
                    n_non_contact=n_neg,
                    min_cells=min_cells,
                    min_replicates=min_replicates,
                    details=str(exc),
                    **meta,
                )

        if source_result:
            results[source_name] = source_result

    return results or None


def compute_pseudobulk_group_de(
    adata,
    groupby: str,
    *,
    replicate: str,
    counts_layer: Optional[str] = "counts",
    min_cells: int = 20,
    min_replicates: int = 2,
    min_pct_expressed: float = 0.0,
    p_adjust_method: str = "fdr_bh",
    padj_cutoff: float = 0.05,
    log2fc_cutoff: float = 0.5,
    fit_type: str = "parametric",
) -> Optional[Dict[str, Dict[str, Dict[str, Any]]]]:
    """Compute pairwise and category-vs-rest pseudobulk DE for one categorical column."""
    if groupby not in adata.obs.columns:
        print(f"  Warning: pseudobulk DE groupby '{groupby}' not found in obs.")
        return None
    if replicate not in adata.obs.columns:
        print(f"  Warning: pseudobulk DE replicate '{replicate}' not found in obs.")
        return None

    col = adata.obs[groupby]
    if pd.api.types.is_numeric_dtype(col):
        print(f"  Warning: pseudobulk DE '{groupby}' is numeric; skipping.")
        return None
    if not isinstance(col.dtype, CategoricalDtype):
        col = col.astype("category")

    categories = [str(cat) for cat in col.cat.categories]
    if len(categories) < 2:
        return None

    rep_values = adata.obs[replicate].astype(str)
    group_values = col.astype(str)
    valid = rep_values.notna().to_numpy() & group_values.notna().to_numpy()
    valid &= np.asarray(col.cat.codes.to_numpy() >= 0, dtype=bool)
    if not valid.any():
        return None

    expression_matrix, counts_layer_used, warning = _as_count_matrix(adata, counts_layer)
    valid_indices = np.flatnonzero(valid)
    rep_valid = rep_values.iloc[valid_indices].astype(str).to_numpy()
    group_valid = group_values.iloc[valid_indices].astype(str).to_numpy()
    unique_reps = sorted(set(rep_valid.astype(str)))
    use_welch = len(unique_reps) == 1

    sample_keys = []
    sample_index = {}
    rows = np.empty(valid_indices.size, dtype=np.int64)
    for i, (rep, grp) in enumerate(zip(rep_valid, group_valid)):
        key = (str(rep), str(grp))
        idx = sample_index.get(key)
        if idx is None:
            idx = len(sample_keys)
            sample_index[key] = idx
            sample_keys.append(key)
        rows[i] = idx

    incidence = sp.csr_matrix(
        (np.ones(valid_indices.size, dtype=np.float64), (rows, valid_indices)),
        shape=(len(sample_keys), adata.n_obs),
    )
    aggregate = incidence @ expression_matrix
    aggregate = _to_dense_counts(aggregate)
    cell_counts = np.bincount(rows, minlength=len(sample_keys)).astype(int)

    pb_meta = pd.DataFrame(
        {
            "_pb_replicate": [key[0] for key in sample_keys],
            "_pb_group": [key[1] for key in sample_keys],
            "n_cells": cell_counts,
        },
        index=[f"pb_{i}" for i in range(len(sample_keys))],
    )
    pb_meta.attrs["gene_names"] = [str(g) for g in adata.var_names]
    aggregate_summary = _compute_category_gene_means_from_aggregate(
        aggregate,
        pb_meta,
        categories,
        adata.var_names,
    )

    results: Dict[str, Dict[str, Dict[str, Any]]] = {}
    category_index = {category: idx for idx, category in enumerate(categories)}
    fit_cache: Dict[Tuple[str, str], Dict[str, Any]] = {}
    n_directed_comparisons = len(categories) * max(len(categories) - 1, 0)
    n_unique_fits = n_directed_comparisons // 2
    if use_welch:
        print(
            f"    - {len(categories)} categories -> {n_directed_comparisons} directed comparisons "
            f"plus {len(categories)} category-vs-rest contrasts "
            "(single groupby replicate; using Welch t-test)",
            flush=True,
        )
    else:
        print(
            f"    - {len(categories)} categories -> {n_directed_comparisons} directed comparisons "
            f"plus {len(categories)} category-vs-rest contrasts "
            f"({n_unique_fits} pairwise DESeq2 fit{'s' if n_unique_fits != 1 else ''})",
            flush=True,
        )
    rest_reference = "__rest__"
    for source in categories:
        source_results: Dict[str, Dict[str, Any]] = {}
        for reference in categories:
            if reference == source:
                continue
            comparison_label = f"{source} vs {reference}"

            source_cell_mask = (group_values.to_numpy() == source) & valid
            reference_cell_mask = (group_values.to_numpy() == reference) & valid
            n_source = int(np.count_nonzero(source_cell_mask))
            n_reference = int(np.count_nonzero(reference_cell_mask))
            if n_source < int(min_cells) or n_reference < int(min_cells):
                print(
                    f"      - {comparison_label}: skipped, insufficient cells "
                    f"({n_source} vs {n_reference}; need >= {int(min_cells)})",
                    flush=True,
                )
                source_results[reference] = _empty_result(
                    "insufficient_cells",
                    n_source=n_source,
                    n_reference=n_reference,
                    min_cells=int(min_cells),
                    min_replicates=int(min_replicates),
                )
                continue

            if use_welch:
                print(
                    f"      - {comparison_label}: fitting Welch t-test "
                    f"({n_source} vs {n_reference} cells, {adata.n_vars} genes)",
                    flush=True,
                )
                pair_result = _fit_welch_pair(
                    expression_matrix,
                    source_cell_mask,
                    reference_cell_mask,
                    adata.var_names,
                )
                source_results[reference] = _format_result(
                    pair_result,
                    source_mask=source_cell_mask,
                    reference_mask=reference_cell_mask,
                    expression_matrix=expression_matrix,
                    var_names=adata.var_names,
                    top_n=0,
                    n_source=n_source,
                    n_reference=n_reference,
                    n_replicates=1,
                    counts_layer_used=counts_layer_used,
                    warning=warning,
                    p_adjust_method=p_adjust_method,
                    min_pct_expressed=min_pct_expressed,
                    padj_cutoff=padj_cutoff,
                    log2fc_cutoff=log2fc_cutoff,
                    method="welch-t-test",
                )
                continue

            source_pb = pb_meta[(pb_meta["_pb_group"] == source) & (pb_meta["n_cells"] >= int(min_cells))]
            reference_pb = pb_meta[(pb_meta["_pb_group"] == reference) & (pb_meta["n_cells"] >= int(min_cells))]
            paired_reps = sorted(set(source_pb["_pb_replicate"]) & set(reference_pb["_pb_replicate"]))
            if len(paired_reps) < int(min_replicates):
                print(
                    f"      - {comparison_label}: skipped, insufficient paired replicates "
                    f"({len(paired_reps)}; need >= {int(min_replicates)})",
                    flush=True,
                )
                source_results[reference] = _empty_result(
                    "insufficient_replicates",
                    n_source=n_source,
                    n_reference=n_reference,
                    min_cells=int(min_cells),
                    min_replicates=int(min_replicates),
                    details=f"{len(paired_reps)} paired replicate(s) available",
                )
                continue

            keep_mask = (
                pb_meta["_pb_replicate"].isin(paired_reps)
                & pb_meta["_pb_group"].isin([source, reference])
                & (pb_meta["n_cells"] >= int(min_cells))
            )
            keep_positions = np.flatnonzero(keep_mask.to_numpy())
            pair_counts = aggregate[keep_positions]
            pair_meta = pb_meta.iloc[keep_positions].copy()
            pair_gene_count = int(pair_counts.shape[1])

            try:
                pair_key = tuple(
                    sorted(
                        (source, reference),
                        key=lambda category: category_index.get(category, 0),
                    )
                )
                cached_fit = fit_cache.get(pair_key)
                if cached_fit is None:
                    print(
                        f"      - {comparison_label}: fitting DESeq2 "
                        f"({len(paired_reps)} paired replicate"
                        f"{'s' if len(paired_reps) != 1 else ''}, "
                        f"{pair_counts.shape[0]} pseudobulk samples, "
                        f"{pair_gene_count} genes)",
                        flush=True,
                    )
                    sample_diagnostics = _compute_pseudobulk_sample_diagnostics(
                        pair_counts,
                        pair_meta,
                    )
                    try:
                        pair_result = _fit_deseq2_pair(
                            pair_counts,
                            pair_meta,
                            source,
                            reference,
                            fit_type=str(fit_type or "parametric"),
                        )
                    except TypeError as exc:
                        if "fit_type" not in str(exc):
                            raise
                        pair_result = _fit_deseq2_pair(
                            pair_counts,
                            pair_meta,
                            source,
                            reference,
                        )
                    cached_fit = {
                        "source": source,
                        "reference": reference,
                        "result": pair_result,
                        "sample_diagnostics": sample_diagnostics,
                    }
                    fit_cache[pair_key] = cached_fit
                else:
                    print(
                        f"      - {comparison_label}: using cached reverse fit "
                        f"({len(paired_reps)} paired replicate"
                        f"{'s' if len(paired_reps) != 1 else ''}, "
                        f"{pair_counts.shape[0]} pseudobulk samples, "
                        f"{pair_gene_count} genes)",
                        flush=True,
                    )
                    pair_result = cached_fit["result"].copy()
                    sample_diagnostics = cached_fit.get("sample_diagnostics")
                    if cached_fit["source"] != source:
                        if "log2FoldChange" in pair_result.columns:
                            pair_result["log2FoldChange"] = -pair_result["log2FoldChange"]
                        if "stat" in pair_result.columns:
                            pair_result["stat"] = -pair_result["stat"]
                source_results[reference] = _format_result(
                    pair_result,
                    source_mask=source_cell_mask,
                    reference_mask=reference_cell_mask,
                    expression_matrix=expression_matrix,
                    var_names=adata.var_names,
                    top_n=0,
                    n_source=n_source,
                    n_reference=n_reference,
                    n_replicates=len(paired_reps),
                    counts_layer_used=counts_layer_used,
                    warning=warning,
                    p_adjust_method=p_adjust_method,
                    min_pct_expressed=min_pct_expressed,
                    padj_cutoff=padj_cutoff,
                    log2fc_cutoff=log2fc_cutoff,
                    sample_diagnostics=sample_diagnostics,
                )
            except Exception as exc:
                print(
                    f"      - {comparison_label}: failed ({exc})",
                    flush=True,
                )
                source_results[reference] = _empty_result(
                    "de_failed",
                    n_source=n_source,
                    n_reference=n_reference,
                    min_cells=int(min_cells),
                    min_replicates=int(min_replicates),
                    details=str(exc),
                )

        comparison_label = f"{source} vs rest"
        source_cell_mask = (group_values.to_numpy() == source) & valid
        reference_cell_mask = (group_values.to_numpy() != source) & valid
        n_source = int(np.count_nonzero(source_cell_mask))
        n_reference = int(np.count_nonzero(reference_cell_mask))
        if n_source < int(min_cells) or n_reference < int(min_cells):
            print(
                f"      - {comparison_label}: skipped, insufficient cells "
                f"({n_source} vs {n_reference}; need >= {int(min_cells)})",
                flush=True,
            )
            source_results[rest_reference] = _empty_result(
                "insufficient_cells",
                n_source=n_source,
                n_reference=n_reference,
                min_cells=int(min_cells),
                min_replicates=int(min_replicates),
            )
        elif use_welch:
            print(
                f"      - {comparison_label}: fitting Welch t-test "
                f"({n_source} vs {n_reference} cells, {adata.n_vars} genes)",
                flush=True,
            )
            rest_result = _fit_welch_pair(
                expression_matrix,
                source_cell_mask,
                reference_cell_mask,
                adata.var_names,
            )
            source_results[rest_reference] = _format_result(
                rest_result,
                source_mask=source_cell_mask,
                reference_mask=reference_cell_mask,
                expression_matrix=expression_matrix,
                var_names=adata.var_names,
                top_n=0,
                n_source=n_source,
                n_reference=n_reference,
                n_replicates=1,
                counts_layer_used=counts_layer_used,
                warning=warning,
                p_adjust_method=p_adjust_method,
                min_pct_expressed=min_pct_expressed,
                padj_cutoff=padj_cutoff,
                log2fc_cutoff=log2fc_cutoff,
                method="welch-t-test",
            )
        else:
            rest_rows: List[np.ndarray] = []
            rest_meta_rows: List[Dict[str, Any]] = []
            paired_reps: List[str] = []
            source_pb = pb_meta[
                (pb_meta["_pb_group"] == source)
                & (pb_meta["n_cells"] >= int(min_cells))
            ]
            for rep in sorted(set(source_pb["_pb_replicate"].astype(str))):
                source_mask_pb = (
                    (pb_meta["_pb_replicate"].astype(str) == rep)
                    & (pb_meta["_pb_group"] == source)
                    & (pb_meta["n_cells"] >= int(min_cells))
                )
                source_positions = np.flatnonzero(source_mask_pb.to_numpy())
                if source_positions.size == 0:
                    continue
                rest_mask_pb = (
                    (pb_meta["_pb_replicate"].astype(str) == rep)
                    & (pb_meta["_pb_group"] != source)
                )
                rest_positions = np.flatnonzero(rest_mask_pb.to_numpy())
                if rest_positions.size == 0:
                    continue
                rest_cells = int(pb_meta.iloc[rest_positions]["n_cells"].sum())
                if rest_cells < int(min_cells):
                    continue
                source_pos = int(source_positions[0])
                paired_reps.append(rep)
                rest_rows.append(np.asarray(aggregate[source_pos], dtype=np.int64))
                rest_meta_rows.append(
                    {
                        "_pb_replicate": rep,
                        "_pb_group": source,
                        "n_cells": int(pb_meta.iloc[source_pos]["n_cells"]),
                    }
                )
                rest_rows.append(
                    np.asarray(
                        aggregate[rest_positions].sum(axis=0),
                        dtype=np.int64,
                    ).ravel()
                )
                rest_meta_rows.append(
                    {
                        "_pb_replicate": rep,
                        "_pb_group": rest_reference,
                        "n_cells": rest_cells,
                    }
                )

            if len(paired_reps) < int(min_replicates):
                print(
                    f"      - {comparison_label}: skipped, insufficient paired replicates "
                    f"({len(paired_reps)}; need >= {int(min_replicates)})",
                    flush=True,
                )
                source_results[rest_reference] = _empty_result(
                    "insufficient_replicates",
                    n_source=n_source,
                    n_reference=n_reference,
                    min_cells=int(min_cells),
                    min_replicates=int(min_replicates),
                    details=f"{len(paired_reps)} paired replicate(s) available",
                )
            else:
                pair_counts = np.vstack(rest_rows)
                pair_meta = pd.DataFrame(
                    rest_meta_rows,
                    index=[f"pb_rest_{i}" for i in range(len(rest_meta_rows))],
                )
                pair_meta.attrs["gene_names"] = [str(g) for g in adata.var_names]
                pair_gene_count = int(pair_counts.shape[1])
                try:
                    print(
                        f"      - {comparison_label}: fitting DESeq2 "
                        f"({len(paired_reps)} paired replicate"
                        f"{'s' if len(paired_reps) != 1 else ''}, "
                        f"{pair_counts.shape[0]} pseudobulk samples, "
                        f"{pair_gene_count} genes)",
                        flush=True,
                    )
                    sample_diagnostics = _compute_pseudobulk_sample_diagnostics(
                        pair_counts,
                        pair_meta,
                    )
                    try:
                        rest_result = _fit_deseq2_pair(
                            pair_counts,
                            pair_meta,
                            source,
                            rest_reference,
                            fit_type=str(fit_type or "parametric"),
                        )
                    except TypeError as exc:
                        if "fit_type" not in str(exc):
                            raise
                        rest_result = _fit_deseq2_pair(
                            pair_counts,
                            pair_meta,
                            source,
                            rest_reference,
                        )
                    source_results[rest_reference] = _format_result(
                        rest_result,
                        source_mask=source_cell_mask,
                        reference_mask=reference_cell_mask,
                        expression_matrix=expression_matrix,
                        var_names=adata.var_names,
                        top_n=0,
                        n_source=n_source,
                        n_reference=n_reference,
                        n_replicates=len(paired_reps),
                        counts_layer_used=counts_layer_used,
                        warning=warning,
                        p_adjust_method=p_adjust_method,
                        min_pct_expressed=min_pct_expressed,
                        padj_cutoff=padj_cutoff,
                        log2fc_cutoff=log2fc_cutoff,
                        sample_diagnostics=sample_diagnostics,
                    )
                except Exception as exc:
                    print(
                        f"      - {comparison_label}: failed ({exc})",
                        flush=True,
                    )
                    source_results[rest_reference] = _empty_result(
                        "de_failed",
                        n_source=n_source,
                        n_reference=n_reference,
                        min_cells=int(min_cells),
                        min_replicates=int(min_replicates),
                        details=str(exc),
                    )

        if source_results:
            results[source] = source_results

    if results:
        results["_summary"] = {
            "category_gene_means": aggregate_summary,
            "replicate": str(replicate),
            "groupby": str(groupby),
            "counts_layer": counts_layer_used,
        }
    return results or None

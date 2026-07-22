import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from karospace.pseudobulk import compute_pseudobulk_interaction_markers


def test_pseudobulk_interaction_markers_aggregate_contact_status_by_replicate(monkeypatch):
    obs = pd.DataFrame(
        {
            "mouse_id": ["M1", "M1", "M1", "M1", "M2", "M2", "M2", "M2"],
            "celltype": pd.Categorical(["A", "B", "A", "C", "A", "B", "A", "C"]),
        },
        index=[f"cell_{i}" for i in range(8)],
    )
    var = pd.DataFrame(index=["G1", "G2"])
    counts = np.array(
        [
            [10, 1],
            [0, 0],
            [1, 8],
            [0, 0],
            [12, 2],
            [0, 0],
            [2, 9],
            [0, 0],
        ],
        dtype=int,
    )
    adata = AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()

    graph = sp.lil_matrix((8, 8), dtype=float)
    graph[0, 1] = 1
    graph[1, 0] = 1
    graph[4, 5] = 1
    graph[5, 4] = 1
    graph = graph.tocsr()

    categories = ["A", "B", "C"]
    labels = np.array([0, 1, 0, 2, 0, 1, 0, 2], dtype=np.int32)
    onehot = sp.csr_matrix(
        (np.ones(labels.size), (np.arange(labels.size), labels)),
        shape=(labels.size, len(categories)),
    )
    neighbor_counts = np.asarray(onehot.T.dot(graph).dot(onehot))

    captured = {}

    def fake_fit(pair_counts, metadata, source, reference, fit_type="parametric"):
        captured["pair_counts"] = pair_counts.copy()
        captured["metadata"] = metadata.copy()
        captured["source"] = source
        captured["reference"] = reference
        return pd.DataFrame(
            {
                "baseMean": [6.0, 5.0],
                "log2FoldChange": [2.0, -1.5],
                "stat": [4.0, -3.0],
                "pvalue": [0.001, 0.2],
                "padj": [0.002, 0.2],
            },
            index=["G1", "G2"],
        )

    monkeypatch.setattr("karospace.pseudobulk._fit_deseq2_pair", fake_fit)

    result = compute_pseudobulk_interaction_markers(
        adata,
        "celltype",
        replicate="mouse_id",
        graph=graph,
        obs_idx=np.arange(adata.n_obs),
        labels=labels,
        categories=categories,
        neighbor_counts=neighbor_counts,
        counts_layer="counts",
        top_targets=1,
        top_genes=2,
        min_cells=1,
        min_neighbors=1,
        min_replicates=2,
        fit_type="mean",
    )

    pair = result["A"]["B"]
    assert pair["available"] is True
    assert pair["method"] == "pseudobulk-deseq2-contact"
    assert pair["n_contact"] == 2
    assert pair["n_non_contact"] == 2
    assert pair["n_replicates"] == 2
    assert pair["genes"] == ["G1", "G2"]
    assert captured["source"] == "contact+"
    assert captured["reference"] == "contact-"
    assert captured["metadata"]["_pb_replicate"].tolist() == ["M1", "M1", "M2", "M2"]
    assert captured["metadata"]["_pb_group"].tolist() == ["contact+", "contact-", "contact+", "contact-"]
    assert captured["pair_counts"].tolist() == [[10, 1], [1, 8], [12, 2], [2, 9]]


def test_pseudobulk_interaction_markers_count_weighted_edges_as_neighbors(monkeypatch):
    obs = pd.DataFrame(
        {
            "mouse_id": ["M1", "M1", "M1", "M1", "M2", "M2", "M2", "M2"],
            "celltype": pd.Categorical(["A", "B", "A", "C", "A", "B", "A", "C"]),
        },
        index=[f"cell_{i}" for i in range(8)],
    )
    var = pd.DataFrame(index=["G1"])
    counts = np.array([[10], [0], [1], [0], [12], [0], [2], [0]], dtype=int)
    adata = AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()

    graph = sp.lil_matrix((8, 8), dtype=float)
    graph[0, 1] = 0.25
    graph[1, 0] = 0.25
    graph[4, 5] = 0.5
    graph[5, 4] = 0.5
    graph = graph.tocsr()

    categories = ["A", "B", "C"]
    labels = np.array([0, 1, 0, 2, 0, 1, 0, 2], dtype=np.int32)
    onehot = sp.csr_matrix(
        (np.ones(labels.size), (np.arange(labels.size), labels)),
        shape=(labels.size, len(categories)),
    )
    neighbor_counts = np.asarray(onehot.T.dot(graph).dot(onehot))

    captured = {}

    def fake_fit(pair_counts, metadata, source, reference, fit_type="parametric"):
        captured["metadata"] = metadata.copy()
        return pd.DataFrame(
            {
                "baseMean": [6.0],
                "log2FoldChange": [2.0],
                "stat": [4.0],
                "pvalue": [0.001],
                "padj": [0.002],
            },
            index=["G1"],
        )

    monkeypatch.setattr("karospace.pseudobulk._fit_deseq2_pair", fake_fit)

    result = compute_pseudobulk_interaction_markers(
        adata,
        "celltype",
        replicate="mouse_id",
        graph=graph,
        obs_idx=np.arange(adata.n_obs),
        labels=labels,
        categories=categories,
        neighbor_counts=neighbor_counts,
        counts_layer="counts",
        top_targets=1,
        top_genes=1,
        min_cells=1,
        min_neighbors=1,
        min_replicates=2,
        fit_type="mean",
    )

    assert result["A"]["B"]["available"] is True
    assert result["A"]["B"]["n_contact"] == 2
    assert result["A"]["B"]["n_replicates"] == 2
    assert captured["metadata"]["_pb_group"].tolist() == ["contact+", "contact-", "contact+", "contact-"]


def test_pseudobulk_interaction_markers_compact_kept_paired_samples(monkeypatch):
    obs = pd.DataFrame(
        {
            "mouse_id": [
                "M0", "M0",
                "M1", "M1", "M1",
                "M2", "M2", "M2",
            ],
            "celltype": pd.Categorical(["A", "B", "A", "B", "A", "A", "B", "A"]),
        },
        index=[f"cell_{i}" for i in range(8)],
    )
    var = pd.DataFrame(index=["G1"])
    counts = np.array([[5], [0], [10], [0], [1], [12], [0], [2]], dtype=int)
    adata = AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()

    graph = sp.lil_matrix((8, 8), dtype=float)
    for a_idx, b_idx in [(0, 1), (2, 3), (5, 6)]:
        graph[a_idx, b_idx] = 1
        graph[b_idx, a_idx] = 1
    graph = graph.tocsr()

    categories = ["A", "B"]
    labels = np.array([0, 1, 0, 1, 0, 0, 1, 0], dtype=np.int32)
    onehot = sp.csr_matrix(
        (np.ones(labels.size), (np.arange(labels.size), labels)),
        shape=(labels.size, len(categories)),
    )
    neighbor_counts = np.asarray(onehot.T.dot(graph).dot(onehot))

    captured = {}

    def fake_fit(pair_counts, metadata, source, reference, fit_type="parametric"):
        captured["pair_counts"] = pair_counts.copy()
        captured["metadata"] = metadata.copy()
        return pd.DataFrame(
            {
                "baseMean": [6.0],
                "log2FoldChange": [2.0],
                "stat": [4.0],
                "pvalue": [0.001],
                "padj": [0.002],
            },
            index=["G1"],
        )

    monkeypatch.setattr("karospace.pseudobulk._fit_deseq2_pair", fake_fit)

    result = compute_pseudobulk_interaction_markers(
        adata,
        "celltype",
        replicate="mouse_id",
        graph=graph,
        obs_idx=np.arange(adata.n_obs),
        labels=labels,
        categories=categories,
        neighbor_counts=neighbor_counts,
        counts_layer="counts",
        top_targets=1,
        top_genes=1,
        min_cells=1,
        min_neighbors=1,
        min_replicates=2,
        fit_type="mean",
    )

    assert result["A"]["B"]["available"] is True
    assert result["A"]["B"]["n_replicates"] == 2
    assert captured["metadata"]["_pb_replicate"].tolist() == ["M1", "M1", "M2", "M2"]
    assert captured["metadata"]["_pb_group"].tolist() == ["contact+", "contact-", "contact+", "contact-"]
    assert captured["pair_counts"].shape == (4, 1)

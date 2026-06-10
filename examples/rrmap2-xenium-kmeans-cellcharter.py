"""
KaroSpace export for the RRMAP2 Xenium "all samples" kmeans-separated object,
annotated + filtered + processed with CellCharter (companion-ready).

Colours by ALL leiden resolutions and ALL CellCharter resolutions; because the
file is companion-ready, the per-cluster analytics (marker genes / DE / neighbor
enrichment) for every clustering are cheap precomputed lookups.

~1.42M cells across 54 samples — writes a binary-sidecar viewer bundle and a
matching packaged .karospace.

Usage:
    python examples/rrmap2-xenium-kmeans-cellcharter.py
    RRMAP2_KMEANS_CELLCHARTER_H5AD=/path/to/file.h5ad python examples/rrmap2-xenium-kmeans-cellcharter.py
"""

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("KMP_WARNINGS", "0")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/karospace-mpl-cache")

from karospace import export_to_html, load_spatial_data

H5AD_PATH = os.environ.get(
    "RRMAP2_KMEANS_CELLCHARTER_H5AD",
    "/Volumes/moldiassd/RRMAP2_xenium_adata/kmeans_separated/"
    "RRMAP2_xenium_all_samples.cellcharter.companion.ready.with_metadata.rerun.h5ad",
)

# Sample-level clinical/experimental metadata — exposed BOTH as section metadata
# and as colour options.
METADATA_COLS = ["stage", "condition", "region", "sex", "model"]

# All clusterings = every CellCharter resolution + every leiden resolution.
CELLCHARTER = [
    "CellCharter_5", "CellCharter_10", "CellCharter_15", "CellCharter_20",
    "CellCharter_25", "CellCharter_30", "CellCharter_35", "CellCharter_40",
    "CellCharter_45", "CellCharter_50",
]
LEIDEN = [
    "leiden_0.5", "leiden_1", "leiden_1.5", "leiden_2",
    "leiden_2.5", "leiden_3", "leiden_3.5", "leiden_4.0",
]
ALL_CLUSTERINGS = LEIDEN + CELLCHARTER

PRIMARY_COLOR = "leiden_2.5"
ADDITIONAL_COLORS = ALL_CLUSTERINGS + METADATA_COLS

SIDECAR_OUTPUT = "rrmap2-xenium-kmeans-cellcharter.html"
PACKAGE_OUTPUT = "rrmap2-xenium-kmeans-cellcharter.karospace"
GENE_AUX_PATH = "rrmap2-xenium-kmeans-cellcharter.genes.json"


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet RRMAP2_KMEANS_CELLCHARTER_H5AD or edit the path above."
        )

    print("Loading spatial data (15 GB / ~1.4M cells — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        groupby="kmeans_split_id",
        spatial_key="spatial",
        metadata_columns=METADATA_COLS,
    )
    print(f"  Loaded {dataset.n_sections} sections with {dataset.n_cells:,} total cells")

    common_kwargs = dict(
        color=PRIMARY_COLOR,
        title="RRMAP2 Xenium All Samples — kmeans / CellCharter",
        min_panel_size=120,
        spot_size="auto",
        downsample=10_000_000,
        theme="light",
        outline_by=None,
        additional_colors=ADDITIONAL_COLORS,
        genes=[],
        use_hvgs=False,
        gene_storage="sidecar",
        gene_sidecar_format="binary-v1",
        gene_encoding="auto",
        gene_value_encoding="uint8",
        gene_aux_path=GENE_AUX_PATH,
        gene_sidecar_shard_size=128,
        # All clusterings get analytics (cheap precomputed companion lookups).
        marker_genes_groupby=ALL_CLUSTERINGS,
        marker_genes_top_n=30,
        neighbor_stats_groupby=ALL_CLUSTERINGS,
        neighbor_stats_permutations=0,
        neighbor_stats_seed=42,
        cluster_de_groupby=ALL_CLUSTERINGS,
        cluster_de_top_n=20,
        cluster_de_method="t-test",
        cluster_de_layer=None,
        cluster_de_min_cells=20,
        interaction_markers_groupby=None,
    )

    print("Exporting binary-sidecar viewer...")
    export_to_html(dataset, output_path=SIDECAR_OUTPUT, **common_kwargs)

    print("Packaging .karospace archive...")
    export_to_html(dataset, output_path=PACKAGE_OUTPUT, **common_kwargs)

    print(f"\nDone!\n  sidecar: {SIDECAR_OUTPUT} (+ {GENE_AUX_PATH} + shard dir)")
    print(f"  package: {PACKAGE_OUTPUT} (+ {Path(PACKAGE_OUTPUT).with_suffix('.loader.html')})")
    print("Coloured by leiden_2.5; switch between all leiden / CellCharter resolutions in the colour selector.")


if __name__ == "__main__":
    main()

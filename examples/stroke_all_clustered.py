"""
Export a KaroSpace viewer for the stroke (dMCAO) dataset stroke_all_clustered —
3 samples, ~261k cells, split into one section per sample.

No embedded morphology images (no uns/spatial); it's a straight spatial cell map
coloured by the leiden clustering, with cell morphology / QC metrics available as
alternative continuous colourings (area, density, elongation, counts, ...).

Usage:
    python examples/stroke_all_clustered.py
    STROKE_H5AD=/path/to/file.h5ad python examples/stroke_all_clustered.py
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
    "STROKE_H5AD",
    "/Users/chrislangseth/Downloads/stroke_all_clustered.companion.ready.h5ad",
)

OUTPUT_DIR = REPO_ROOT
SIDECAR_OUTPUT = str(OUTPUT_DIR / "stroke_all_clustered.html")
PACKAGE_OUTPUT = str(OUTPUT_DIR / "stroke_all_clustered.karospace")
GENE_AUX_PATH = str(OUTPUT_DIR / "stroke_all_clustered.genes.json")

PRIMARY_COLOR = "leiden"
ADDITIONAL_COLORS = [
    "leiden",
    "total_counts",
    "n_genes_by_counts",
    "density",
    "area",
    "elongation",
    "avg_confidence",
]
CLUSTER_COLUMNS = ["leiden"]


def main() -> None:
    if not Path(H5AD_PATH).exists():
        raise SystemExit(
            f"H5AD not found: {H5AD_PATH}\nSet STROKE_H5AD or edit the path above."
        )

    print("Loading spatial data (large file — this can take a while)...")
    dataset = load_spatial_data(
        H5AD_PATH,
        groupby="sample",
        spatial_key="spatial",
        metadata_section=["sample"],
    )
    print(f"  {dataset.n_sections} section(s), {dataset.n_cells:,} cells")
    print(f"  section IDs: {[s.section_id for s in dataset.sections]}")

    common_kwargs = dict(
        annotation=PRIMARY_COLOR,
        title="Stroke (dMCAO) — leiden clusters",
        min_panel_size=140,
        spot_size="auto",
        outline_by=None,
        cells_annotations=ADDITIONAL_COLORS,
        genes=[],
        use_hvgs=False,
        gene_storage="sidecar",
        gene_encoding="auto",
        gene_value_encoding="uint8",
        gene_sidecar_shard_size=128,
        marker_genes_groupby=CLUSTER_COLUMNS,
        marker_genes_top_n=30,
        neighbor_stats_groupby=CLUSTER_COLUMNS,
        neighbor_stats_permutations=0,
        cluster_de_groupby=CLUSTER_COLUMNS,
        cluster_de_top_n=20,
        cluster_de_method="t-test",
        cluster_de_layer=None,
        interaction_markers_groupby=None,
    )

    # The sidecar writes its gene manifest next to the HTML (full path is fine);
    # the .karospace packager requires a bare filename (it lives inside the archive).
    print("Exporting sidecar HTML...")
    export_to_html(
        dataset, output_path=SIDECAR_OUTPUT, gene_aux_path=GENE_AUX_PATH, **common_kwargs
    )

    print("Packaging .karospace archive...")
    export_to_html(
        dataset,
        output_path=PACKAGE_OUTPUT,
        gene_aux_path=Path(GENE_AUX_PATH).name,
        **common_kwargs,
    )

    print(f"Done!\n  {SIDECAR_OUTPUT}\n  {PACKAGE_OUTPUT}")
    print("Coloured by leiden; switch to morphology/QC metrics via the colour selector.")


if __name__ == "__main__":
    main()

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon


def read_results(results_csv: Path):
    df = pd.read_csv(results_csv)
    required = {"node_id", "x", "y", "ux", "uy", "umag"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Results CSV missing columns: {sorted(missing)}")
    return df.sort_values("node_id").reset_index(drop=True)


def read_connectivity(connectivity_csv: Path):
    df = pd.read_csv(connectivity_csv)
    required = {"type", "n1", "n2", "n3", "n4"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Connectivity CSV missing columns: {sorted(missing)}")
    return df


def get_polygon_node_ids(row):
    etype = str(row["type"]).strip()
    if etype == "Tri3":
        return [int(row["n1"]), int(row["n2"]), int(row["n3"])]
    elif etype == "Quad4":
        return [int(row["n1"]), int(row["n2"]), int(row["n3"]), int(row["n4"])]
    else:
        raise ValueError(f"Unsupported element type: {etype}")


def polygon_segments(points):
    segs = []
    n = len(points)
    for i in range(n):
        segs.append([points[i], points[(i + 1) % n]])
    return segs


def main():
    parser = argparse.ArgumentParser(description="Plot static FEM results for Tri3/Quad4 meshes.")
    parser.add_argument("results_csv", type=str)
    parser.add_argument("connectivity_csv", type=str)
    parser.add_argument("--scale", type=float, default=1.0, help="Deformation scale factor")
    parser.add_argument("--output", type=str, default="")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--cmap", type=str, default="viridis")
    parser.add_argument("--alpha", type=float, default=0.85)
    parser.add_argument("--coord-unit", type=str, default="m", choices=["m", "mm"], help="Displayed coordinate unit")
    parser.add_argument("--disp-unit", type=str, default="mm", choices=["m", "mm"], help="Displayed displacement unit")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    results_df = read_results(Path(args.results_csv))
    connectivity_df = read_connectivity(Path(args.connectivity_csv))

    coord_scale = 1.0 if args.coord_unit == "m" else 1000.0
    disp_scale = 1.0 if args.disp_unit == "m" else 1000.0

    node_ids = results_df["node_id"].to_numpy(dtype=int)

    coords = results_df[["x", "y"]].to_numpy(dtype=float) * coord_scale
    disp = results_df[["ux", "uy"]].to_numpy(dtype=float) * coord_scale
    umag = results_df["umag"].to_numpy(dtype=float) * disp_scale

    coords_map = {int(node_ids[i]): coords[i] for i in range(len(node_ids))}
    disp_map = {int(node_ids[i]): disp[i] for i in range(len(node_ids))}
    umag_map = {int(node_ids[i]): float(umag[i]) for i in range(len(node_ids))}

    undeformed_lines = []
    deformed_lines = []
    patches = []
    patch_values = []

    for _, row in connectivity_df.iterrows():
        node_list = get_polygon_node_ids(row)

        pts_undef = [coords_map[nid] for nid in node_list]
        pts_def = [coords_map[nid] + args.scale * disp_map[nid] for nid in node_list]

        undeformed_lines.extend(polygon_segments(pts_undef))
        deformed_lines.extend(polygon_segments(pts_def))

        patches.append(Polygon(np.array(pts_def), closed=True))
        patch_values.append(np.mean([umag_map[nid] for nid in node_list]))

    undeformed_lines = np.array(undeformed_lines, dtype=float)
    deformed_lines = np.array(deformed_lines, dtype=float)
    deformed_coords = coords + args.scale * disp

    fig, ax = plt.subplots(figsize=(11, 7))

    pc = PatchCollection(
        patches,
        cmap=args.cmap,
        alpha=args.alpha,
        edgecolor="none",
        linewidth=0.0,
        zorder=1
    )
    pc.set_array(np.array(patch_values))
    ax.add_collection(pc)

    lc_undeformed = LineCollection(
        undeformed_lines,
        colors="0.75",
        linewidths=0.8,
        zorder=2
    )
    ax.add_collection(lc_undeformed)

    lc_deformed = LineCollection(
        deformed_lines,
        colors="k",
        linewidths=1.0,
        zorder=3
    )
    ax.add_collection(lc_deformed)

    cbar = fig.colorbar(pc, ax=ax, pad=0.02)
    cbar.set_label(f"Element-averaged displacement magnitude [{args.disp_unit}]")

    legend_handles = [
        Line2D([0], [0], color="0.75", lw=1.5, label="Undeformed mesh"),
        Line2D([0], [0], color="k", lw=1.5, label=f"Deformed mesh (scale = {args.scale:g})"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=True)

    all_pts = np.vstack([coords, deformed_coords])
    xmin, ymin = np.min(all_pts, axis=0)
    xmax, ymax = np.max(all_pts, axis=0)

    dx = xmax - xmin
    dy = ymax - ymin
    pad_x = 0.05 * dx if dx > 0 else 1.0
    pad_y = 0.05 * dy if dy > 0 else 1.0

    ax.set_xlim(xmin - pad_x, xmax + pad_x)
    ax.set_ylim(ymin - pad_y, ymax + pad_y)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"x [{args.coord_unit}]")
    ax.set_ylabel(f"y [{args.coord_unit}]")
    ax.set_title("Static FEM Result")
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.35)

    fig.tight_layout()

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight")

    if not args.no_show:
        plt.show()

    plt.close(fig)


if __name__ == "__main__":
    main()
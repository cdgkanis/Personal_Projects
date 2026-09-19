import argparse
from pathlib import Path

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colors


def build_segments(results_df, conn_df, scale=1.0):
    undeformed_segments = []
    deformed_segments = []
    deformed_values = []

    node_map = {}
    for _, row in results_df.iterrows():
        node_map[int(row["node_id"])] = {
            "x": float(row["x"]),
            "y": float(row["y"]),
            "ux": float(row["ux"]),
            "uy": float(row["uy"]),
            "umag": float(row["umag"]),
        }

    for _, elem in conn_df.iterrows():
        elem_type = str(elem["type"]).strip()

        node_ids = []
        for key in ["n1", "n2", "n3", "n4"]:
            if key in elem and pd.notna(elem[key]):
                nid = int(elem[key])
                if nid >= 0:
                    node_ids.append(nid)

        if elem_type == "Tri3" and len(node_ids) == 3:
            node_ids.append(node_ids[0])
        elif elem_type == "Quad4" and len(node_ids) == 4:
            node_ids.append(node_ids[0])
        else:
            continue

        undeformed_poly = []
        deformed_poly = []
        poly_umag = []

        for nid in node_ids:
            if nid not in node_map:
                raise RuntimeError(f"Node ID {nid} from connectivity not found in results CSV.")

            nd = node_map[nid]
            x = nd["x"]
            y = nd["y"]
            ux = nd["ux"]
            uy = nd["uy"]
            umag = nd["umag"]

            undeformed_poly.append((x, y))
            deformed_poly.append((x + scale * ux, y + scale * uy))
            poly_umag.append(umag)

        for i in range(len(undeformed_poly) - 1):
            undeformed_segments.append([undeformed_poly[i], undeformed_poly[i + 1]])
            deformed_segments.append([deformed_poly[i], deformed_poly[i + 1]])
            deformed_values.append(0.5 * (poly_umag[i] + poly_umag[i + 1]))

    return undeformed_segments, deformed_segments, deformed_values


def plot_mesh(ax, results_df, conn_df, scale):
    undeformed_segments, deformed_segments, deformed_values = build_segments(results_df, conn_df, scale)

    if not deformed_segments:
        raise RuntimeError("No valid element segments were generated for plotting.")

    undeformed_lc = LineCollection(
        undeformed_segments,
        colors="0.7",
        linewidths=1.0,
        linestyles="dashed",
        label="Undeformed"
    )
    ax.add_collection(undeformed_lc)

    vmin = min(deformed_values)
    vmax = max(deformed_values)
    if abs(vmax - vmin) < 1e-16:
        vmax = vmin + 1e-16

    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.colormaps["viridis"]

    deformed_lc = LineCollection(
        deformed_segments,
        cmap=cmap,
        norm=norm,
        linewidths=2.0
    )
    deformed_lc.set_array(pd.Series(deformed_values).to_numpy())
    ax.add_collection(deformed_lc)

    cbar = plt.colorbar(deformed_lc, ax=ax)
    cbar.set_label("Displacement magnitude")

    xs = []
    ys = []
    for _, row in results_df.iterrows():
        xs.append(float(row["x"]))
        ys.append(float(row["y"]))
        xs.append(float(row["x"]) + scale * float(row["ux"]))
        ys.append(float(row["y"]) + scale * float(row["uy"]))

    ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(min(ys), max(ys))
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"Deformed shape (scale = {scale:g})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, linestyle="--", alpha=0.35)


def plot_history(ax, history_df):
    if "displacement" not in history_df.columns or "load" not in history_df.columns:
        raise RuntimeError("History CSV must contain 'displacement' and 'load' columns.")

    ax.plot(
        history_df["displacement"],
        history_df["load"],
        marker="o",
        linewidth=1.8,
        markersize=4.5,
        color="tab:red"
    )
    ax.set_title("Load-displacement diagram")
    ax.set_xlabel("Displacement")
    ax.set_ylabel("Load")
    ax.grid(True, linestyle="--", alpha=0.35)

    if "step" in history_df.columns:
        for _, row in history_df.iterrows():
            ax.annotate(
                str(int(row["step"])),
                (row["displacement"], row["load"]),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=8,
                alpha=0.75
            )


def main():
    parser = argparse.ArgumentParser(description="Plot static FE results and optional load-displacement history.")
    parser.add_argument("results_csv", help="CSV file with nodal results.")
    parser.add_argument("connectivity_csv", help="CSV file with element connectivity.")
    parser.add_argument("--scale", type=float, default=10.0, help="Deformation scale factor.")
    parser.add_argument("--history", type=str, default=None, help="Optional load-history CSV.")
    parser.add_argument("--output", type=str, default="static_plot.png", help="Output image path.")
    parser.add_argument("--no-show", action="store_true", help="Do not display the plot window.")
    args = parser.parse_args()

    results_path = Path(args.results_csv)
    conn_path = Path(args.connectivity_csv)

    if not results_path.exists():
        raise FileNotFoundError(f"Results CSV not found: {results_path}")
    if not conn_path.exists():
        raise FileNotFoundError(f"Connectivity CSV not found: {conn_path}")
    if args.scale <= 0.0:
        raise ValueError("Scale factor must be positive.")

    results_df = pd.read_csv(results_path)
    conn_df = pd.read_csv(conn_path)

    required_results_cols = {"node_id", "x", "y", "ux", "uy", "umag"}
    required_conn_cols = {"element_id", "type", "n1", "n2", "n3", "n4"}

    if not required_results_cols.issubset(results_df.columns):
        missing = required_results_cols - set(results_df.columns)
        raise RuntimeError(f"Results CSV missing columns: {sorted(missing)}")

    if not required_conn_cols.issubset(conn_df.columns):
        missing = required_conn_cols - set(conn_df.columns)
        raise RuntimeError(f"Connectivity CSV missing columns: {sorted(missing)}")

    if args.history:
        history_path = Path(args.history)
        if not history_path.exists():
            raise FileNotFoundError(f"History CSV not found: {history_path}")

        history_df = pd.read_csv(history_path)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        plot_mesh(axes[0], results_df, conn_df, args.scale)
        plot_history(axes[1], history_df)
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_mesh(ax, results_df, conn_df, args.scale)

    fig.tight_layout()
    fig.savefig(args.output, dpi=200, bbox_inches="tight")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter
from matplotlib.collections import LineCollection


def parse_args():
    parser = argparse.ArgumentParser(description="Animate transient FEM solution.")
    parser.add_argument("folder", help="Folder containing results and geometry CSV files")
    parser.add_argument("--results", default="transient_results.csv", help="Results CSV filename")
    parser.add_argument("--nodes", default="nodes.csv", help="Nodes CSV filename")
    parser.add_argument("--connectivity", default="connectivity.csv", help="Connectivity CSV filename")
    parser.add_argument("--scale", type=float, default=1.0, help="Deformation scale factor")
    parser.add_argument("--fps", type=int, default=30, help="Animation frames per second")
    parser.add_argument("--max-frames", type=int, default=300, help="Maximum number of frames to export/display")
    parser.add_argument("--dpi", type=int, default=120, help="Output DPI for saved animation")
    parser.add_argument("--figsize", type=float, nargs=2, default=[8.0, 8.0], metavar=("W", "H"),
                        help="Figure size in inches")
    parser.add_argument("--save-gif", action="store_true", help="Save animation as GIF")
    parser.add_argument("--save-mp4", action="store_true", help="Save animation as MP4")
    parser.add_argument("--no-show", action="store_true", help="Do not open interactive window")
    parser.add_argument("--color-by-amplitude", action="store_true",
                        help="Color deformed edges by average nodal displacement magnitude")
    return parser.parse_args()


def validate_inputs(folder: Path, results_name: str, nodes_name: str, conn_name: str):
    results_path = folder / results_name
    nodes_path = folder / nodes_name
    conn_path = folder / conn_name

    if not results_path.exists():
        raise FileNotFoundError(f"Results file not found: {results_path}")
    if not nodes_path.exists():
        raise FileNotFoundError(f"Nodes file not found: {nodes_path}")
    if not conn_path.exists():
        raise FileNotFoundError(f"Connectivity file not found: {conn_path}")

    return results_path, nodes_path, conn_path


def load_data(results_path: Path, nodes_path: Path, conn_path: Path):
    res = pd.read_csv(results_path)
    nodes = pd.read_csv(nodes_path)
    conn = pd.read_csv(conn_path)

    if "time" not in res.columns:
        raise ValueError("Results CSV must contain 'time' column.")
    if not {"node_id", "x", "y"}.issubset(nodes.columns):
        raise ValueError("Nodes CSV must contain node_id, x, y columns.")
    if not {"n1", "n2"}.issubset(conn.columns):
        raise ValueError("Connectivity CSV must contain n1, n2 columns.")

    nodes = nodes.sort_values("node_id").reset_index(drop=True)

    time = res["time"].to_numpy(dtype=float)
    xy0 = nodes[["x", "y"]].to_numpy(dtype=float)
    node_ids = nodes["node_id"].to_numpy(dtype=int)
    n_nodes = len(node_ids)

    expected_u_cols = [f"u_{i}" for i in range(2 * n_nodes)]
    missing = [c for c in expected_u_cols if c not in res.columns]
    if missing:
        raise ValueError(
            f"Missing displacement columns. Example missing columns: {missing[:10]}"
        )

    U = res[expected_u_cols].to_numpy(dtype=float)

    edge_pairs = conn[["n1", "n2"]].to_numpy(dtype=int)

    if np.any(edge_pairs < 0) or np.any(edge_pairs >= n_nodes):
        raise ValueError("Connectivity contains node indices outside node table range.")

    return time, xy0, U, edge_pairs


def build_frame_indices(n_steps: int, max_frames: int):
    if n_steps <= 0:
        raise ValueError("No time steps found in results.")
    if max_frames <= 0:
        raise ValueError("max_frames must be positive.")

    if n_steps <= max_frames:
        return list(range(n_steps))

    stride = int(np.ceil(n_steps / max_frames))
    frames = list(range(0, n_steps, stride))
    if frames[-1] != n_steps - 1:
        frames.append(n_steps - 1)
    return frames


def compute_bounds(xy0: np.ndarray, U: np.ndarray, scale: float):
    xmin0, ymin0 = np.min(xy0, axis=0)
    xmax0, ymax0 = np.max(xy0, axis=0)

    n_nodes = xy0.shape[0]
    ux = U[:, 0:2 * n_nodes:2]
    uy = U[:, 1:2 * n_nodes:2]

    x_all_min = min(xmin0, np.min(xy0[:, 0][None, :] + scale * ux))
    x_all_max = max(xmax0, np.max(xy0[:, 0][None, :] + scale * ux))
    y_all_min = min(ymin0, np.min(xy0[:, 1][None, :] + scale * uy))
    y_all_max = max(ymax0, np.max(xy0[:, 1][None, :] + scale * uy))

    span = max(x_all_max - x_all_min, y_all_max - y_all_min, 1.0)
    pad = 0.10 * span

    return (
        x_all_min - pad,
        x_all_max + pad,
        y_all_min - pad,
        y_all_max + pad,
    )


def build_segments_for_frame(xy0: np.ndarray, U: np.ndarray, edge_pairs: np.ndarray, frame: int, scale: float):
    n_nodes = xy0.shape[0]
    u = U[frame]

    xy = xy0.copy()
    xy[:, 0] += scale * u[0:2 * n_nodes:2]
    xy[:, 1] += scale * u[1:2 * n_nodes:2]

    segments = np.stack([xy[edge_pairs[:, 0]], xy[edge_pairs[:, 1]]], axis=1)
    return xy, segments


def nodal_magnitude_for_frame(U: np.ndarray, frame: int, n_nodes: int):
    u = U[frame]
    ux = u[0:2 * n_nodes:2]
    uy = u[1:2 * n_nodes:2]
    return np.sqrt(ux * ux + uy * uy)


def main():
    args = parse_args()

    if args.scale <= 0.0:
        raise ValueError("scale must be positive.")
    if args.fps <= 0:
        raise ValueError("fps must be positive.")
    if args.dpi <= 0:
        raise ValueError("dpi must be positive.")

    folder = Path(args.folder)
    results_path, nodes_path, conn_path = validate_inputs(folder, args.results, args.nodes, args.connectivity)
    time, xy0, U, edge_pairs = load_data(results_path, nodes_path, conn_path)

    n_nodes = xy0.shape[0]
    frame_indices = build_frame_indices(len(time), args.max_frames)
    undeformed_segments = np.stack([xy0[edge_pairs[:, 0]], xy0[edge_pairs[:, 1]]], axis=1)

    xlim_l, xlim_r, ylim_b, ylim_t = compute_bounds(xy0, U, args.scale)

    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    ax.set_aspect("equal")
    ax.set_xlim(xlim_l, xlim_r)
    ax.set_ylim(ylim_b, ylim_t)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.25)

    undeformed_lc = LineCollection(
        undeformed_segments,
        colors="0.75",
        linewidths=1.0,
        linestyles="dashed"
    )
    ax.add_collection(undeformed_lc)

    if args.color_by_amplitude:
        mags0 = nodal_magnitude_for_frame(U, frame_indices[0], n_nodes)
        edge_mag0 = 0.5 * (mags0[edge_pairs[:, 0]] + mags0[edge_pairs[:, 1]])
        deformed_lc = LineCollection(
            undeformed_segments,
            cmap="viridis",
            linewidths=2.0
        )
        deformed_lc.set_array(edge_mag0)
        deformed_lc.set_clim(vmin=0.0, vmax=max(1e-16, np.max(np.sqrt(U[:, 0:2*n_nodes:2]**2 + U[:, 1:2*n_nodes:2]**2))))
        ax.add_collection(deformed_lc)
        cbar = fig.colorbar(deformed_lc, ax=ax, shrink=0.85)
        cbar.set_label("|u|")
    else:
        deformed_lc = LineCollection(
            undeformed_segments,
            colors="tab:blue",
            linewidths=2.0
        )
        ax.add_collection(deformed_lc)

    title = ax.set_title("")
    info_text = ax.text(
        0.02, 0.98, "",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="0.8")
    )

    def set_frame(step_index):
        frame = frame_indices[step_index]
        xy, segments = build_segments_for_frame(xy0, U, edge_pairs, frame, args.scale)
        deformed_lc.set_segments(segments)

        mags = nodal_magnitude_for_frame(U, frame, n_nodes)
        max_mag = float(np.max(mags))
        rms_mag = float(np.sqrt(np.mean(mags ** 2)))

        if args.color_by_amplitude:
            edge_mag = 0.5 * (mags[edge_pairs[:, 0]] + mags[edge_pairs[:, 1]])
            deformed_lc.set_array(edge_mag)

        title.set_text(f"Transient deformation   t = {time[frame]:.6f} s")
        info_text.set_text(
            f"frame {step_index + 1}/{len(frame_indices)}\n"
            f"scale = {args.scale:.3g}\n"
            f"max|u| = {max_mag:.6e}\n"
            f"rms|u| = {rms_mag:.6e}"
        )

        return deformed_lc, title, info_text

    def init():
        return set_frame(0)

    def update(step_index):
        return set_frame(step_index)

    ani = FuncAnimation(
        fig,
        update,
        frames=len(frame_indices),
        init_func=init,
        interval=1000.0 / args.fps,
        blit=True,
        repeat=False
    )

    saved_any = False

    if args.save_mp4:
        ffmpeg_ok = shutil.which("ffmpeg") is not None
        if not ffmpeg_ok:
            raise RuntimeError("ffmpeg was not found in PATH. Install ffmpeg or use --save-gif.")

        mp4_path = folder / "transient_animation.mp4"
        writer = FFMpegWriter(fps=args.fps, bitrate=1800)
        ani.save(mp4_path, writer=writer, dpi=args.dpi)
        print(f"Saved MP4: {mp4_path}")
        saved_any = True

    if args.save_gif:
        gif_path = folder / "transient_animation.gif"
        gif_fps = min(args.fps, 15)
        ani.save(gif_path, writer=PillowWriter(fps=gif_fps), dpi=min(args.dpi, 100))
        print(f"Saved GIF: {gif_path}")
        saved_any = True

    if not args.no_show and not saved_any:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
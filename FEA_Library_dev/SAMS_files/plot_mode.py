import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def read_mode_csv(csv_path):
    meta = {}
    rows = []

    with open(csv_path, "r", newline="") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("#"):
                s = s[1:].strip()
                parts = [p.strip() for p in s.split(",", 1)]
                if len(parts) == 2:
                    meta[parts[0]] = parts[1]
                continue

            if s.startswith("node_id,"):
                header = [h.strip() for h in s.split(",")]
                break
        else:
            raise RuntimeError("Mode CSV header not found.")

        reader = csv.DictReader(f, fieldnames=header)
        for row in reader:
            rows.append(row)

    data = {
        "node_id": np.array([int(r["node_id"]) for r in rows], dtype=int),
        "x": np.array([float(r["x"]) for r in rows], dtype=float),
        "y": np.array([float(r["y"]) for r in rows], dtype=float),
        "ux": np.array([float(r["ux"]) for r in rows], dtype=float),
        "uy": np.array([float(r["uy"]) for r in rows], dtype=float),
        "umag": np.array([float(r["umag"]) for r in rows], dtype=float),
    }
    return meta, data


def read_connectivity_csv(csv_path):
    elems = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            typ = row["type"].strip()
            n1 = int(row["n1"])
            n2 = int(row["n2"])
            n3 = int(row["n3"])
            n4 = int(row["n4"])
            elems.append((typ, [n1, n2, n3, n4]))
    return elems


def build_triangles(elems, id_to_idx):
    tris = []

    for typ, conn in elems:
        if typ == "Tri3":
            tris.append([
                id_to_idx[conn[0]],
                id_to_idx[conn[1]],
                id_to_idx[conn[2]]
            ])
        elif typ == "Quad4":
            tris.append([
                id_to_idx[conn[0]],
                id_to_idx[conn[1]],
                id_to_idx[conn[2]]
            ])
            tris.append([
                id_to_idx[conn[0]],
                id_to_idx[conn[2]],
                id_to_idx[conn[3]]
            ])
        else:
            raise RuntimeError(f"Unsupported element type: {typ}")

    return np.array(tris, dtype=int)


def compute_scale(x, y, ux, uy, scale_factor=0.10):
    span = max(np.max(x) - np.min(x), np.max(y) - np.min(y), 1e-12)
    peak = max(np.max(np.abs(ux)), np.max(np.abs(uy)), 1e-12)
    return scale_factor * span / peak


def plot_mode(mode_csv, conn_csv, scale_factor=0.10):
    meta, data = read_mode_csv(mode_csv)
    elems = read_connectivity_csv(conn_csv)

    node_id = data["node_id"]
    x = data["x"]
    y = data["y"]
    ux = data["ux"]
    uy = data["uy"]
    umag = data["umag"]

    id_to_idx = {nid: i for i, nid in enumerate(node_id)}
    triangles = build_triangles(elems, id_to_idx)

    scale = compute_scale(x, y, ux, uy, scale_factor)
    xd = x + scale * ux
    yd = y + scale * uy

    triang = mtri.Triangulation(xd, yd, triangles)

    fig, ax = plt.subplots(figsize=(10, 8))

    tpc = ax.tripcolor(triang, umag, shading="gouraud", cmap="viridis")
    ax.triplot(triang, color="k", linewidth=0.5, alpha=0.5)

    cbar = fig.colorbar(tpc, ax=ax)
    cbar.set_label("Modal displacement magnitude")

    ax.set_title(f"Mode {meta.get('mode_index', '?')}   f = {meta.get('frequency_hz', '?')} Hz")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.2)

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode_csv")
    parser.add_argument("conn_csv")
    parser.add_argument("--scale", type=float, default=0.10)
    args = parser.parse_args()

    plot_mode(args.mode_csv, args.conn_csv, args.scale)
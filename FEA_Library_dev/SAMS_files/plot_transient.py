import sys
print(sys.executable)
print(sys.path)

import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def plot_transient_results(csv_file, dof):
    df = pd.read_csv(csv_file)

    time_col = "time"
    u_col = f"u_{dof}"
    v_col = f"v_{dof}"
    a_col = f"a_{dof}"

    required = [time_col, u_col, v_col, a_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in CSV: {missing}")

    t = df[time_col]

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(t, df[u_col], color="tab:blue", linewidth=1.5)
    axes[0].set_ylabel("Displacement")
    axes[0].set_title(f"Transient response for DOF {dof}")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, df[v_col], color="tab:orange", linewidth=1.5)
    axes[1].set_ylabel("Velocity")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, df[a_col], color="tab:green", linewidth=1.5)
    axes[2].set_ylabel("Acceleration")
    axes[2].set_xlabel("Time [s]")
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot transient FEM results from CSV.")
    parser.add_argument("csv_file", type=str, help="Path to CSV file")
    parser.add_argument("--dof", type=int, default=0, help="DOF index to plot")
    args = parser.parse_args()

    plot_transient_results(args.csv_file, args.dof)
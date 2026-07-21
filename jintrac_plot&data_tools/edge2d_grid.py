"""Plot the EDGE2D grid and separatrix stored in a run's TRAN file.

Each EDGE2D cell is stored in RVERTP/ZVERTP as a group of five points:
four vertices followed by the cell centroid. The separatrix is read from
RSEPX/ZSEPX.
"""

import os
import sys

sys.path.append('/home/lgarzot/Work/Python/JET') #jintrac python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #jintrac python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import eproc as ep
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from available_data_run import find_run_directory


def ask_run_name(prompt="Run name: "):
    """Ask for a non-empty run name."""
    while True:
        run_name = input(prompt).strip()
        if run_name:
            return run_name
        print("The run name cannot be empty.")


def read_ep_array(tran_file, signal):
    """Read the valid points of one EPROC signal as a flat NumPy array."""
    signal_data = ep.data(tran_file, signal)
    return np.asarray(signal_data.data[:signal_data.nPts]).reshape(-1)


def load_edge2d_geometry(tran_file):
    """Load cell vertices, centroids and separatrix coordinates."""
    rvertp = read_ep_array(tran_file, "RVERTP")
    zvertp = -read_ep_array(tran_file, "ZVERTP")

    number_of_cells = min(len(rvertp), len(zvertp)) // 5
    if number_of_cells == 0:
        raise ValueError("No EDGE2D cells were found in RVERTP/ZVERTP.")

    # Ignore incomplete trailing data, if present, and reshape each quintuple.
    cell_points = np.column_stack(
        [rvertp[: 5 * number_of_cells], zvertp[: 5 * number_of_cells]]
    ).reshape(number_of_cells, 5, 2)
    vertices = cell_points[:, :4, :]
    centroids = cell_points[:, 4, :]

    r_sep = read_ep_array(tran_file, "RSEPX")
    z_sep = -read_ep_array(tran_file, "ZSEPX")
    number_of_separatrix_points = min(len(r_sep), len(z_sep))
    separatrix = np.column_stack(
        [
            r_sep[:number_of_separatrix_points],
            z_sep[:number_of_separatrix_points],
        ]
    )

    valid_cells = np.isfinite(vertices).all(axis=(1, 2)) & np.isfinite(
        centroids
    ).all(axis=1)
    vertices = vertices[valid_cells]
    centroids = centroids[valid_cells]
    separatrix = separatrix[np.isfinite(separatrix).all(axis=1)]

    if len(vertices) == 0:
        raise ValueError("The EDGE2D grid does not contain finite cell coordinates.")
    if len(separatrix) == 0:
        raise ValueError("RSEPX/ZSEPX do not contain finite separatrix points.")

    return vertices, centroids, separatrix


def plot_edge2d_grid(run_name):
    """Plot and return the EDGE2D grid figure for one run."""
    run_path = find_run_directory(run_name)
    tran_file = os.path.join(run_path, "tran")
    if not os.path.exists(tran_file):
        raise FileNotFoundError(
            f"TRAN file not found for run '{run_name}'; the EDGE2D grid is unavailable."
        )

    vertices, centroids, separatrix = load_edge2d_geometry(tran_file)

    # Close every quadrilateral by repeating its first vertex.
    closed_cells = np.concatenate([vertices, vertices[:, :1, :]], axis=1)
    all_vertices = vertices.reshape(-1, 2)

    figure, axis = plt.subplots(figsize=(10, 10))
    cell_collection = LineCollection(
        closed_cells,
        colors="0.35",
        linewidths=0.45,
        alpha=0.9,
        label="EDGE2D cells",
    )
    axis.add_collection(cell_collection)

    axis.scatter(
        all_vertices[:, 0],
        all_vertices[:, 1],
        s=3,
        color="tab:blue",
        alpha=0.65,
        label="Cell vertices",
        zorder=2,
    )
    axis.scatter(
        centroids[:, 0],
        centroids[:, 1],
        s=8,
        color="tab:orange",
        label="Cell centroids",
        zorder=3,
    )

    axis.plot(
        separatrix[:, 0],
        separatrix[:, 1],
        color="red",
        linewidth=1.5,
        marker="o",
        markersize=3,
        label="Separatrix (RSEPX/ZSEPX)",
        zorder=4,
    )

    axis.autoscale()
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("R (m)")
    axis.set_ylabel("Z (m)")
    axis.set_title(f"EDGE2D grid - {run_name}")
    axis.grid(True, alpha=0.2)
    axis.legend(loc="best")
    figure.tight_layout()

    print(f"Cells plotted: {len(vertices)}")
    print(f"Separatrix points plotted: {len(separatrix)}")
    return figure, axis


def main():
    run_name = sys.argv[1].strip() if len(sys.argv) > 1 else ask_run_name()
    try:
        plot_edge2d_grid(run_name)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}")
        return 1

    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

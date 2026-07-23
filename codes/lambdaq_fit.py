"""
Author: Rachele Cicioni
last update: 21-07-2026

Calculate lambda_q and plot the EDGE2D QPARTOT map.

The QPARTOT profile is extracted along the outer mid-plane Z=Z0. Only the
points outside the outer separatrix are fitted, using the radial distance

    r = R - R_sep.

Lambda_q is obtained from the exponential decay law defined in the OMP (Outer Mid Plane):
!! The OMP is defined as the horizontal line passing through the magnetic axis, Z=Z0.
    q(r) = q0 * exp(-r / lambda_q).

This is in agreement with T.Eich NF paper (2013) 53 093031 in which the definition of lambda_q
is considered in the OMP.
Taking the logarithm gives

    ln|q(r)| = ln|q0| - r / lambda_q.

A linear fit of ln|q| versus r therefore has slope m = -1/lambda_q, so the
final value is calculated as lambda_q = -1/m.
"""

import os
import sys

sys.path.append('/home/lgarzot/Work/Python/JET') #jintrac python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #jintrac python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import eproc as ep
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree

from available_data_run import ask_run_name, find_run_directory

PROFILE = "QPARTOT"

def read_signal(tran_file, signal):
    # Read the valid points of an EPROC signal.
    data = ep.data(tran_file, signal)
    return np.asarray(data.data[:data.nPts]).reshape(-1)

def load_edge2d_data(tran_file):
    # Read mesh coordinates using the plotting Z convention.
    rmesh = read_signal(tran_file, "RMESH")
    zmesh = -read_signal(tran_file, "ZMESH") #negative Z
    rvertp = read_signal(tran_file, "RVERTP")
    zvertp = -read_signal(tran_file, "ZVERTP") #negative Z
    qpartot = read_signal(tran_file, PROFILE)

    # Each cell contains four vertices and one centroid.
    number_of_cells = min(len(rvertp), len(zvertp)) // 5
    if number_of_cells == 0:
        raise ValueError("No EDGE2D cells found in RVERTP/ZVERTP.")

    points = np.column_stack(
        [rvertp[: 5 * number_of_cells], zvertp[: 5 * number_of_cells]]
    ).reshape(number_of_cells, 5, 2)
    vertices = points[:, :4, :]
    centroids = points[:, 4, :]

    # Remove cells with invalid geometry.
    valid_cells = np.isfinite(vertices).all(axis=(1, 2)) & np.isfinite(
        centroids
    ).all(axis=1)
    vertices = vertices[valid_cells]
    centroids = centroids[valid_cells]
    if len(vertices) == 0:
        raise ValueError("The EDGE2D mesh has no valid cells.")

    # Match each cell centroid to the nearest QPARTOT mesh point.
    mesh_points = np.column_stack([rmesh, zmesh])
    _, profile_indices = cKDTree(mesh_points).query(centroids)
    if np.max(profile_indices) >= len(qpartot):
        raise ValueError("QPARTOT does not contain one value per mesh point.")
    cell_values = qpartot[profile_indices]

    # Read the magnetic axis and separatrix.
    r0 = float(read_signal(tran_file, "R0")[0])
    z0 = -float(read_signal(tran_file, "Z0")[0])
    r_sep = read_signal(tran_file, "RSEPX")
    z_sep = -read_signal(tran_file, "ZSEPX")
    number_of_sep_points = min(len(r_sep), len(z_sep))
    separatrix = np.column_stack(
        [r_sep[:number_of_sep_points], z_sep[:number_of_sep_points]]
    )
    separatrix = separatrix[np.isfinite(separatrix).all(axis=1)]
    if len(separatrix) < 2:
        raise ValueError("RSEPX/ZSEPX do not define a valid separatrix.")

    return {
        "vertices": vertices,
        "centroids": centroids,
        "qpartot": cell_values,
        "separatrix": separatrix,
        "r0": r0,
        "z0": z0,
    }

def outer_separatrix_radius(separatrix, z0):
    # Intersect every separatrix segment with the line Z=Z0.
    intersections = []
    for index, point_1 in enumerate(separatrix):
        point_2 = separatrix[(index + 1) % len(separatrix)]
        r1, z1 = point_1
        r2, z2 = point_2

        if np.isclose(z1, z0) and np.isclose(z2, z0):
            intersections.extend([r1, r2])
            continue
        if np.isclose(z1, z2):
            continue

        fraction = (z0 - z1) / (z2 - z1)
        if -1e-9 <= fraction <= 1.0 + 1e-9:
            intersections.append(r1 + fraction * (r2 - r1))

    if not intersections:
        raise ValueError("The separatrix does not intersect the line Z=Z0.")

    # The largest intersection radius is the outer separatrix.
    return float(np.max(intersections))

def cross2d(vector_a, vector_b):
    # Scalar cross product of two 2D vectors.
    return vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0]

def line_segment_intersection(origin, direction, point_1, point_2):
    # Find the line parameter at an intersection with a finite segment.
    segment = point_2 - point_1
    denominator = cross2d(direction, segment)
    if abs(denominator) < 1e-12:
        return None

    offset = point_1 - origin
    t_value = cross2d(offset, segment) / denominator
    segment_fraction = cross2d(offset, direction) / denominator
    if -1e-9 <= segment_fraction <= 1.0 + 1e-9:
        return t_value
    return None

def line_cell_interval(origin, direction, vertices):
    # Find where the outer-midplane line enters and leaves one cell.
    intersections = []
    for index, point_1 in enumerate(vertices):
        point_2 = vertices[(index + 1) % len(vertices)]
        t_value = line_segment_intersection(origin, direction, point_1, point_2)
        if t_value is not None:
            intersections.append(t_value)

    if len(intersections) < 2:
        return None

    t_enter = min(intersections)
    t_exit = max(intersections)
    if t_exit - t_enter < 1e-10:
        return None
    return t_enter, t_exit

def outer_midplane_profile(data):
    # Extract QPARTOT along the horizontal line Z=Z0.
    r0 = data["r0"]
    z0 = data["z0"]
    r_max = float(np.max(data["vertices"][:, :, 0]))
    origin = np.array([r0, z0])
    direction = np.array([r_max - r0, 0.0])
    if direction[0] <= 0.0:
        raise ValueError("The outer-midplane line has non-positive length.")

    radial_positions = []
    heat_flux = []
    for vertices, value in zip(data["vertices"], data["qpartot"]):
        interval = line_cell_interval(origin, direction, vertices)
        if interval is None:
            continue

        t_enter, t_exit = interval
        radius = r0 + 0.5 * (t_enter + t_exit) * direction[0]
        if r0 <= radius <= r_max:
            radial_positions.append(radius)
            heat_flux.append(abs(value))

    if not radial_positions:
        raise ValueError("No mesh cells intersect the outer mid-plane Z=Z0.")

    order = np.argsort(radial_positions)
    return np.asarray(radial_positions)[order], np.asarray(heat_flux)[order]

def fit_lambda_q(radial_positions, heat_flux, r_sep):
    # Keep only finite SOL data outside the outer separatrix.
    sol = (
        (radial_positions > r_sep)
        & (heat_flux > 0.0)
        & np.isfinite(radial_positions)
        & np.isfinite(heat_flux)
    )
    distance = radial_positions[sol] - r_sep
    sol_heat_flux = heat_flux[sol]
    if len(distance) < 3:
        raise ValueError("At least three points beyond the separatrix are required.")

    # ln(q) = ln(q0) - x/lambda_q.
    slope, intercept = np.polyfit(distance, np.log(sol_heat_flux), 1)
    if slope >= 0.0:
        raise ValueError("QPARTOT does not decay beyond the separatrix.")

    lambda_q = -1.0 / slope
    q0 = np.exp(intercept)
    return lambda_q, q0, distance, sol_heat_flux

def plot_qpartot(data, run_name, r_sep):
    # Plot the 2D QPARTOT map, separatrix and outer mid-plane.
    figure, axis = plt.subplots(figsize=(10, 9))
    cells = PolyCollection(
        data["vertices"],
        array=np.asarray(data["qpartot"]),
        cmap="viridis",
        edgecolors="black",
        linewidths=0.25,
    )
    axis.add_collection(cells)
    figure.colorbar(cells, ax=axis, label="QPARTOT (W/m²)")

    separatrix = data["separatrix"]
    axis.plot(
        separatrix[:, 0],
        separatrix[:, 1],
        "r.-",
        linewidth=1.3,
        markersize=3,
        label="Separatrix",
    )
    axis.axhline(data["z0"], color="orange", linestyle="--", label="OMP: Z=Z0")
    axis.plot(r_sep, data["z0"], "ro", label="Outer separatrix")

    axis.autoscale()
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("R (m)")
    axis.set_ylabel("Z (m)")
    axis.set_title(f"QPARTOT - {run_name}")
    axis.legend(loc="best")
    figure.tight_layout()
    return figure, axis


def plot_lambda_q_fit(
    radial_positions,
    heat_flux,
    r_sep,
    distance,
    sol_heat_flux,
    q0,
    lambda_q,
    run_name,
):
    # Plot the OMP profile and the exponential fit in the SOL.
    relative_radius = radial_positions - r_sep
    valid_profile = (
        (heat_flux > 0.0)
        & np.isfinite(relative_radius)
        & np.isfinite(heat_flux)
    )

    fit_radius = np.linspace(0.0, np.max(distance), 300)
    fitted_heat_flux = q0 * np.exp(-fit_radius / lambda_q)

    figure, axis = plt.subplots(figsize=(9, 6))
    axis.semilogy(
        relative_radius[valid_profile],
        heat_flux[valid_profile],
        "o-",
        color="0.65",
        markersize=4,
        linewidth=1,
        label="OMP profile",
    )
    axis.semilogy(
        distance,
        sol_heat_flux,
        "o",
        color="tab:blue",
        markersize=6,
        label="SOL data used in fit",
    )
    axis.semilogy(
        fit_radius,
        fitted_heat_flux,
        "-",
        color="tab:red",
        linewidth=2,
        label=rf"Fit: $\lambda_q={lambda_q * 1000.0:.3f}$ mm",
    )
    axis.axvline(0.0, color="black", linestyle="--", label="Separatrix")
    axis.set_xlabel(r"$r=R-R_{\mathrm{sep}}$ (m)")
    axis.set_ylabel(r"$|QPARTOT|$ (W/m$^2$)")
    axis.set_title(f"Outer-midplane QPARTOT fit - {run_name}")
    axis.grid(True, which="both", alpha=0.3)
    axis.legend(loc="best")
    figure.tight_layout()
    return figure, axis


def analyze_run(run_name):
    # Load the run and calculate lambda_q.
    run_path = find_run_directory(run_name)
    tran_file = os.path.join(run_path, "tran")
    if not os.path.exists(tran_file):
        raise FileNotFoundError(f"TRAN file not found for run '{run_name}'.")

    data = load_edge2d_data(tran_file)
    r_sep = outer_separatrix_radius(data["separatrix"], data["z0"])
    radial_positions, heat_flux = outer_midplane_profile(data)
    lambda_q, q0, distance, sol_heat_flux = fit_lambda_q(
        radial_positions, heat_flux, r_sep
    )
    plot_qpartot(data, run_name, r_sep)
    plot_lambda_q_fit(
        radial_positions,
        heat_flux,
        r_sep,
        distance,
        sol_heat_flux,
        q0,
        lambda_q,
        run_name,
    )
    return lambda_q, r_sep, data["z0"]

def main():
    run_name = sys.argv[1].strip() if len(sys.argv) > 1 else ask_run_name()
    try:
        lambda_q, r_sep, z0 = analyze_run(run_name)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}")
        return 1

    print(f"Z0 = {z0:.6f} m")
    print(f"Outer separatrix at Z0: R = {r_sep:.6f} m")
    print(f"lambda_q = {lambda_q * 1000.0:.3f} mm")
    plt.show()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

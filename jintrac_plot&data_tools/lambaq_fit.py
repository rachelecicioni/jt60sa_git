"""
author: R. Cicioni
last updated: 21-07-2026

Calculate lambda_q from the outer-midplane QPARTOT (TOT.PAR.ENERGY FLUX [W/m^2]) profile of an EDGE2D run.

The script reconstructs the QPARTOT values crossed by the Z=0 outer-midplane line, keeps the scrape-off-layer points (R > R_sep), and fits
|q_parallel(x)| = q0 * exp(-x / lambda_q),  x = R - R_sep.

Only lambda_q is reported; no plots are produced.
"""
import os
import sys

sys.path.append('/home/lgarzot/Work/Python/JET') #jintrac python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #jintrac python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import eproc as ep
import numpy as np
from scipy.spatial import cKDTree

PROFILE = "QPARTOT"

def ask_run_name(prompt="Run name: "):
    #Ask for a non-empty run name
    while True:
        run_name = input(prompt).strip()
        if run_name:
            return run_name
        print("The run name cannot be empty.")


def find_run_directory(run_dir):
    #Find a run directory in the standard JDC and EFGW locations
    potential_paths = [
        f"/common/cmg/jnv7243/jetto/runs/{run_dir}",
        f"/pfs/work/g2rcicio/jetto/runs/{run_dir}",
        f"/common/cmg/jnv7243/edge2d/runs/{run_dir}",
        f"/pfs/work/g2rcicio/edge2d/runs/{run_dir}",
    ]

    for path in potential_paths:
        if os.path.isdir(path):
            print(f"Run found at: {path}")
            return path

    raise FileNotFoundError(
        f"Run '{run_dir}' not found in any of the expected locations:\n"
        + "\n".join(potential_paths)
    )


def read_ep_array(tran_file, signal):
    #Read the valid data points of one EPROC signal as a NumPy array
    signal_data = ep.data(tran_file, signal)
    return np.asarray(signal_data.data[:signal_data.nPts])


def cross2d(vector_a, vector_b):
    #Return the scalar cross product of two 2D vectors
    return vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0]


def line_segment_intersection(origin, direction, point_1, point_2):
    #Return the line parameter t at an intersection with a finite segment
    segment = point_2 - point_1
    denominator = cross2d(direction, segment)
    if abs(denominator) < 1e-12:
        return None

    offset = point_1 - origin
    t_value = cross2d(offset, segment) / denominator
    segment_parameter = cross2d(offset, direction) / denominator
    if -1e-9 <= segment_parameter <= 1.0 + 1e-9:
        return t_value
    return None


def line_quad_interval(origin, direction, quad_vertices):
    #Return the interval in which an infinite line crosses a quadrilateral
    intersections = []
    for index, point_1 in enumerate(quad_vertices):
        point_2 = quad_vertices[(index + 1) % len(quad_vertices)]
        t_value = line_segment_intersection(
            origin, direction, point_1, point_2
        )
        if t_value is not None:
            intersections.append(t_value)

    if len(intersections) < 2:
        return None

    t_enter = min(intersections)
    t_exit = max(intersections)
    if t_exit - t_enter < 1e-10:
        return None
    return t_enter, t_exit


def outer_midplane_profile(tran_file):
    #Return R, |QPARTOT| and R_sep along the Z=0 outer midplane.
    rmesh = read_ep_array(tran_file, "RMESH")
    zmesh = -read_ep_array(tran_file, "ZMESH")
    rvertp = read_ep_array(tran_file, "RVERTP")
    zvertp = -read_ep_array(tran_file, "ZVERTP")
    profile = read_ep_array(tran_file, PROFILE)

    r0 = float(read_ep_array(tran_file, "R0")[0])
    r_sep = float(read_ep_array(tran_file, "RSEPX")[0])
    r_max = float(np.max(rmesh))

    number_of_quads = min(len(rvertp), len(zvertp)) // 5
    if number_of_quads == 0:
        raise ValueError("No EDGE2D quadrilateral cells were found.")

    # RVERTP/ZVERTP use groups of five points: four vertices and one centroid.
    all_points = np.column_stack(
        [rvertp[: 5 * number_of_quads], zvertp[: 5 * number_of_quads]]
    ).reshape(number_of_quads, 5, 2)
    quad_vertices = all_points[:, :4, :]
    centroids = all_points[:, 4, :]

    # Associate every cell centroid with the nearest RMESH/ZMESH profile point.
    mesh_points = np.column_stack([rmesh, zmesh])
    _, profile_indices = cKDTree(mesh_points).query(centroids)
    if np.max(profile_indices) >= len(profile):
        raise ValueError(f"{PROFILE} does not contain one value per mesh point.")
    quad_values = profile[profile_indices]

    origin = np.array([r0, 0.0])
    direction = np.array([r_max - r0, 0.0])
    if direction[0] <= 0.0:
        raise ValueError("The outer-midplane direction has non-positive length.")

    radial_positions = []
    heat_flux_values = []
    for vertices, value in zip(quad_vertices, quad_values):
        interval = line_quad_interval(origin, direction, vertices)
        if interval is None:
            continue

        t_enter, t_exit = interval
        radial_midpoint = origin[0] + 0.5 * (t_enter + t_exit) * direction[0]
        if r0 - 0.01 <= radial_midpoint <= r_max + 0.01:
            radial_positions.append(radial_midpoint)
            heat_flux_values.append(abs(value))

    if not radial_positions:
        raise ValueError("No mesh cells intersect the Z=0 outer midplane.")

    radial_positions = np.asarray(radial_positions)
    heat_flux_values = np.asarray(heat_flux_values)
    sort_indices = np.argsort(radial_positions)
    return (
        radial_positions[sort_indices],
        heat_flux_values[sort_indices],
        r_sep,
    )


def fit_lambda_q(distance, heat_flux):
    """Fit ln(|q_parallel|) versus distance and return lambda_q in metres."""
    valid = (
        (distance > 0.0)
        & (heat_flux > 0.0)
        & np.isfinite(distance)
        & np.isfinite(heat_flux)
    )
    distance = np.asarray(distance)[valid]
    heat_flux = np.asarray(heat_flux)[valid]

    if len(distance) < 3:
        raise ValueError("At least three valid SOL points are required for the fit.")

    slope, _ = np.polyfit(distance, np.log(heat_flux), 1)
    if slope >= 0.0:
        raise ValueError(
            "The fitted heat-flux profile does not decay: the slope is non-negative."
        )

    return -1.0 / slope


def calculate_lambda_q(run_dir, skip_first_sol_point=True):
    """Calculate lambda_q for one EDGE2D run and return its value in metres."""
    run_path = find_run_directory(run_dir)
    tran_file = os.path.join(run_path, "tran")
    if not os.path.exists(tran_file):
        raise FileNotFoundError(
            f"TRAN file not found for run '{run_dir}'; lambda_q requires EDGE2D data."
        )

    radial_positions, heat_flux, r_sep = outer_midplane_profile(tran_file)
    sol_mask = radial_positions > r_sep
    distance = radial_positions[sol_mask] - r_sep
    sol_heat_flux = heat_flux[sol_mask]

    # Preserve the original script's choice to exclude the first SOL cell.
    if skip_first_sol_point and len(distance) > 0:
        distance = distance[1:]
        sol_heat_flux = sol_heat_flux[1:]

    return fit_lambda_q(distance, sol_heat_flux)


def main():
    run_name = sys.argv[1].strip() if len(sys.argv) > 1 else ask_run_name()
    try:
        lambda_q_m = calculate_lambda_q(run_name)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}")
        return 1

    print(f"lambda_q = {lambda_q_m:.6f} m")
    print(f"lambda_q = {lambda_q_m * 1000.0:.3f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

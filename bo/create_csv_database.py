"""
Author: Rachele Cicioni
last update: 23-07-2026

Build a CSV database from a manually defined list of JINTRAC runs.

Edit SIMULATIONS below, then run:

    python create_csv_database.py

An alternative output path can be passed as the first command-line argument:

    python create_csv_database.py my_database.csv

Units:
    D_puff, Ar_seeding, CHI_SOL: units used in the input files
    lambda_q: m
    qpar_max: W/m^2
    Pcore, Psol, Pin: values are combined in their native, consistent units
"""
import csv
import os
from pathlib import Path
import re
import sys

# Add JINTRAC Python tools.
sys.path.append("/home/lgarzot/Work/Python/JET")
sys.path.append("/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python")
sys.path.append("/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs")

# Import the lambda_q calculation without running its main function.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "codes"))

import eproc as ep
import jetto_binary_tools
import numpy as np

from lambdaq_fit import (
    fit_lambda_q,
    load_edge2d_data,
    outer_midplane_profile,
    outer_separatrix_radius,
)

# Add one dictionary for every simulation.
SIMULATIONS = [
    # {"run_id": "run_sa_fk_c_p20_a3_ArW_fc_farsol"},
    # {"run_id": "another_run_name"},
]


CSV_COLUMNS = [
    "run_id",
    "D_puff",
    "Ar_seeding",
    "CHI_SOL",
    "chi_ETB_over_D_ETB",
    "lambda_q",
    "qpar_max",
    "f_rad",
    "status",
]

FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
VALID_STATUSES = {
    "running": "Running",
    "completed successfully": "Completed successfully",
    "failed": "Failed",
}


def find_run_directory(run_id):
    potential_paths = [
        f"/common/cmg/jnv7243/edge2d/runs/{run_id}",
        f"/pfs/work/g2rcicio/edge2d/runs/{run_id}",
        f"/common/cmg/jnv7243/jetto/runs/{run_id}",
        f"/pfs/work/g2rcicio/jetto/runs/{run_id}",
    ]

    for path in potential_paths:
        if os.path.isdir(path):
            return path

    raise FileNotFoundError(
        f"Run '{run_id}' not found in:\n" + "\n".join(potential_paths)
    )


def parse_text_value(file_path, variable):
    # Read the first numeric value written next to a variable.
    pattern = re.compile(
        rf"^\s*{re.escape(variable)}\s*(?:=|:)?\s*({FLOAT_PATTERN})",
        re.IGNORECASE,
    )

    with open(file_path, "r", encoding="utf-8", errors="replace") as file:
        for line in file:
            line = line.split("!", 1)[0]
            match = pattern.search(line)
            if match:
                return float(match.group(1).replace("D", "E").replace("d", "e"))

    raise KeyError(f"Variable '{variable}' not found in {file_path}")


def find_f2d_file(run_path):
    # JINTRAC runs may use either f2d or jetto.f2d.
    for filename in ("f2d", "jetto.f2d"):
        file_path = os.path.join(run_path, filename)
        if os.path.isfile(file_path):
            return file_path
    raise FileNotFoundError(f"No f2d or jetto.f2d file found in {run_path}")


def read_status(run_path):
    # The first line must be: Status : <run status>.
    status_path = os.path.join(run_path, "hfps.status")
    with open(status_path, "r", encoding="utf-8", errors="replace") as file:
        first_line = file.readline().strip()

    match = re.fullmatch(r"Status\s*:\s*(.+)", first_line, re.IGNORECASE)
    if not match:
        raise ValueError(f"Invalid first line in {status_path}: {first_line!r}")

    status = match.group(1).strip()
    canonical_status = VALID_STATUSES.get(status.lower())
    if canonical_status is None:
        raise ValueError(f"Unknown run status in {status_path}: {status!r}")
    return canonical_status


def last_jst_value(jst, variable):
    # Return the value at the last JST time step.
    values = np.asarray(jst[variable], dtype=float).reshape(-1)
    if values.size == 0:
        raise ValueError(f"JST/{variable} is empty")
    return float(values[-1])


def last_tran_slice(tran_file, signal):
    # Return all values belonging to the last TRAN time step.
    result = ep.time(tran_file, signal)
    values = np.asarray(result.yData, dtype=float)
    times = np.asarray(result.xData).reshape(-1)

    if values.size == 0:
        raise ValueError(f"TRAN/{signal} is empty")
    if values.ndim == 0:
        return values.reshape(1)
    if values.ndim == 1:
        return values[-1:]

    # Identify the time axis from the number of xData points.
    time_axes = [axis for axis, size in enumerate(values.shape) if size == len(times)]
    time_axis = time_axes[0] if time_axes else 0
    return np.asarray(np.take(values, -1, axis=time_axis)).reshape(-1)


def calculate_chi_ratio(jst):
    # Calculate (JST/XIBA)/(JST/DBAR) at the last time step.
    xiba = last_jst_value(jst, "XIBA")
    dbar = last_jst_value(jst, "DBAR")
    if dbar == 0.0:
        raise ZeroDivisionError("JST/DBAR is zero at the last time step")
    return xiba / dbar


def calculate_lambda_q(tran_file):
    # Reuse the numerical calculation from codes/lambdaq_fit.py.
    data = load_edge2d_data(tran_file)
    r_sep = outer_separatrix_radius(data["separatrix"], data["z0"])
    radial_positions, heat_flux = outer_midplane_profile(data)
    lambda_q, _, _, _ = fit_lambda_q(radial_positions, heat_flux, r_sep)
    return lambda_q


def calculate_qpar_max(tran_file):
    # Use the maximum magnitude of QPARTOT at the last time step.
    final_qpartot = last_tran_slice(tran_file, "QPARTOT")
    return float(np.nanmax(np.abs(final_qpartot)))


def calculate_f_rad(jst, tran_file):
    # f_rad = (Pcore + Psol) / Pin at the last time step.
    pcore = last_jst_value(jst, "PRAD")
    pin = last_jst_value(jst, "PAUX")
    if pin == 0.0:
        raise ZeroDivisionError("JST/PAUX is zero at the last time step")

    psol = 0.0
    for signal in ("HRAD", "ZRAD_1", "ZRAD_2"):
        psol += float(np.nansum(last_tran_slice(tran_file, signal)))

    return (pcore + psol) / pin


def read_optional(run_id, field, calculation):
    # Keep building the row when an output is unavailable.
    try:
        value = calculation()
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            raise ValueError("result is not finite")
        return value
    except Exception as error:
        print(f"Warning [{run_id}] {field}: {error}")
        return None


def build_run_row(simulation):
    # Collect all requested quantities for one run.
    run_id = simulation["run_id"].strip()
    if not run_id:
        raise ValueError("run_id cannot be empty")

    run_path = find_run_directory(run_id)
    coset_file = os.path.join(run_path, "edge2d.coset")
    tran_file = os.path.join(run_path, "tran")

    row = {column: None for column in CSV_COLUMNS}
    row["run_id"] = run_id
    row["D_puff"] = read_optional(
        run_id,
        "D_puff",
        lambda: parse_text_value(coset_file, "EPuffHPanel.HPUFF"),
    )
    row["Ar_seeding"] = read_optional(
        run_id,
        "Ar_seeding",
        lambda: parse_text_value(coset_file, "EPuffHPanel.ZPUFFT[0]"),
    )
    row["CHI_SOL"] = read_optional(
        run_id,
        "CHI_SOL",
        lambda: parse_text_value(find_f2d_file(run_path), "XI_SOLETB"),
    )
    row["status"] = read_optional(
        run_id,
        "status",
        lambda: read_status(run_path),
    )

    # Read JST once and reuse it for all JST-derived quantities.
    jst = read_optional(
        run_id,
        "jetto.jst",
        lambda: jetto_binary_tools.read_binary_file(
            os.path.join(run_path, "jetto.jst")
        ),
    )
    if jst is not None:
        row["chi_ETB_over_D_ETB"] = read_optional(
            run_id,
            "chi_ETB_over_D_ETB",
            lambda: calculate_chi_ratio(jst),
        )

    if os.path.isfile(tran_file):
        row["lambda_q"] = read_optional(
            run_id,
            "lambda_q",
            lambda: calculate_lambda_q(tran_file),
        )
        row["qpar_max"] = read_optional(
            run_id,
            "qpar_max",
            lambda: calculate_qpar_max(tran_file),
        )
        if jst is not None:
            row["f_rad"] = read_optional(
                run_id,
                "f_rad",
                lambda: calculate_f_rad(jst, tran_file),
            )
    else:
        print(f"Warning [{run_id}] TRAN file not found")

    return row


def validate_simulations():
    # Validate the manually edited list before reading run data.
    if not SIMULATIONS:
        raise ValueError("Add at least one {'run_id': '...'} to SIMULATIONS")

    run_ids = []
    for simulation in SIMULATIONS:
        if set(simulation) != {"run_id"}:
            raise ValueError(
                "Each SIMULATIONS entry must contain only the 'run_id' field"
            )
        run_ids.append(simulation["run_id"])

    if len(run_ids) != len(set(run_ids)):
        raise ValueError("SIMULATIONS contains duplicate run_id values")


def write_csv(rows, output_path):
    # Write one row per simulation.
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return output_path.resolve()


def main():
    output_path = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(__file__).with_name("simulation_database.csv")
    )

    try:
        validate_simulations()
        rows = []
        for simulation in SIMULATIONS:
            run_id = simulation["run_id"]
            print(f"Processing {run_id}...")
            try:
                rows.append(build_run_row(simulation))
            except Exception as error:
                print(f"Error [{run_id}]: {error}")
                empty_row = {column: None for column in CSV_COLUMNS}
                empty_row["run_id"] = run_id
                rows.append(empty_row)

        database_path = write_csv(rows, output_path)
    except (OSError, ValueError) as error:
        print(f"Error: {error}")
        return 1

    print(f"CSV database written to: {database_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
R. Cicioni, 20-07-2026

Inspect and extract data from JETTO/EDGE2D simulation runs.

JST and JSP files are required. Every SST file present in the run is detected
automatically, while TRAN is optional because it is absent from JETTO-only runs.

Command-line usage:
    $ python available_data.py <run_name>

Module usage:
    >>> from available_data import available_data, get_data
    >>> outputs = available_data("run_sa_nclass2")
    >>> x, y = get_data("JST", "PRAD", "run_sa_nclass2")
"""

import os
import re
import sys

sys.path.append('/home/lgarzot/Work/Python/JET') #python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import numpy as np
import jetto_binary_tools
import eproc as ep


def ask_run_name(prompt="Run name: "):
    """Ask for a non-empty run name."""
    while True:
        run_name = input(prompt).strip()
        if run_name:
            return run_name
        print("The run name cannot be empty.")


def find_run_directory(run_dir):
    """Find a run directory in the standard JDC and EFGW locations."""
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


def find_sst_files(run_path):
    """Return all numbered jetto.sst files, sorted by SST number."""
    sst_files = []
    for filename in os.listdir(run_path):
        match = re.fullmatch(r"jetto\.sst(\d+)", filename, re.IGNORECASE)
        if match:
            sst_files.append((int(match.group(1)), os.path.join(run_path, filename)))
    return sorted(sst_files, key=lambda item: item[0])


def available_data(run_dir):
    """Print all available variables and return the loaded run outputs."""
    run_path = find_run_directory(run_dir)
    outputs = {}

    # JST and JSP are required for every supported run.
    jst_path = os.path.join(run_path, "jetto.jst")
    jsp_path = os.path.join(run_path, "jetto.jsp")
    outputs["JST"] = jetto_binary_tools.read_binary_file(jst_path)

    # SST files are optional and their number is not fixed.
    for sst_number, sst_path in find_sst_files(run_path):
        outputs[f"SST{sst_number}"] = jetto_binary_tools.read_binary_file(sst_path)

    outputs["JSP"] = jetto_binary_tools.read_binary_file(jsp_path)

    # TRAN is optional and normally absent from JETTO-only runs.
    tran_path = os.path.join(run_path, "tran")
    tran_exists = os.path.exists(tran_path)
    if tran_exists:
        try:
            outputs["TRAN"] = ep.names(tran_path)
        except Exception as exc:
            print(f"Warning: could not read TRAN file: {exc}")

    print("\nAvailable data:")
    for data_type, data in outputs.items():
        if data_type == "TRAN":
            names = data
        else:
            names = list(data.keys())
        print(f"{data_type}: {names}")

    if not tran_exists:
        print("TRAN: not available (this may be a JETTO-only run)")

    return outputs


def get_data(data_type, variable_name, run_dir, sst_number=None, timestep=-1):
    """Extract x and y arrays for one variable from a run output file."""
    run_path = find_run_directory(run_dir)
    data_type = data_type.upper()

    if data_type == "JST":
        data = jetto_binary_tools.read_binary_file(
            os.path.join(run_path, "jetto.jst")
        )
        x = np.asarray(data["TVEC1"]).flatten()
        y = np.asarray(data[variable_name]).flatten()

    elif data_type == "JSP":
        data = jetto_binary_tools.read_binary_file(
            os.path.join(run_path, "jetto.jsp")
        )
        x = np.asarray(data["XRHO"])[timestep, :].flatten()
        y = np.asarray(data[variable_name])[timestep, :].flatten()

    elif data_type.startswith("SST"):
        suffix = data_type[3:]
        if suffix:
            sst_number = int(suffix)
        elif sst_number is None:
            raise ValueError("Specify the SST number, for example 'SST1'.")

        data = jetto_binary_tools.read_binary_file(
            os.path.join(run_path, f"jetto.sst{sst_number}")
        )
        x = np.asarray(data["TVEC1"]).flatten()
        y = np.asarray(data[variable_name]).flatten()

    elif data_type == "TRAN":
        tran_path = os.path.join(run_path, "tran")
        if not os.path.exists(tran_path):
            raise FileNotFoundError(
                f"TRAN file not found for run '{run_dir}' (JETTO-only run)."
            )
        result = ep.time(tran_path, variable_name)
        x = np.asarray(result.xData)
        y = np.asarray(result.yData)

    else:
        raise ValueError(f"Unsupported data type: {data_type}")

    return x, y


def list_tran_signals(run_dir):
    """Return TRAN signal names, or None when the run has no TRAN file."""
    run_path = find_run_directory(run_dir)
    tran_path = os.path.join(run_path, "tran")
    if not os.path.exists(tran_path):
        print("TRAN: not available (this may be a JETTO-only run)")
        return None

    signals = ep.names(tran_path)
    print(f"TRAN: {signals}")
    return signals


def main():
    run_dir = sys.argv[1].strip() if len(sys.argv) > 1 else ask_run_name()
    try:
        available_data(run_dir)
    except FileNotFoundError as exc:
        print(exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

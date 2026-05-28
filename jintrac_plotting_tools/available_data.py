"""
R. Cicioni, 25-05-2025

Utility to load and inspect data from JETTO/EDGE2D simulation runs.
Reads binary files (JST, SST1, SST2, SST3, JSP) and TRAN files, prints available data keys/signals.

Example usage:
    $ python available_data.py <run_directory>
Or:
    >>> from available_data import available_data
    >>> available_data("run_sa_nclass2")
"""

import sys
import os

sys.path.append('/home/lgarzot/Work/Python/JET') #python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import numpy as np
import jetto_binary_tools
import eproc as ep

def find_run_directory(run_dir):
    """Find run directory in edge2d and jetto run directories (JDC or EFGW server)."""
    potential_paths = [
        f"/common/cmg/jnv7243/jetto/runs/{run_dir}",
        f"/pfs/work/g2rcicio/jetto/runs/{run_dir}",
        f"/common/cmg/jnv7243/edge2d/runs/{run_dir}",
        f"/pfs/work/g2rcicio/edge2d/runs/{run_dir}",
    ]
    
    for path in potential_paths:
        if os.path.exists(path):
            print(f"Run found at: {path}")
            return path
    
    raise FileNotFoundError(
        f"Run '{run_dir}' not found in any of the expected locations:\n" + 
        "\n".join(potential_paths)
    )

def available_data(run_dir):
    # --- FILE PATHS ---
    base_path = find_run_directory(run_dir)

    fnamet = f"{base_path}/jetto.jst"   # file JST
    fnamep = f"{base_path}/jetto.jsp"   # file JSP
    ftran = f"{base_path}/tran"   # file TRAN

    # --- READING FILE ---
    jst = jetto_binary_tools.read_binary_file(fnamet)
    
    sst1 = None
    fnamet_imp1 = f"{base_path}/jetto.sst1"
    if os.path.exists(fnamet_imp1):
        sst1 = jetto_binary_tools.read_binary_file(fnamet_imp1)
    
    sst2 = None
    fnamet_imp2 = f"{base_path}/jetto.sst2"
    if os.path.exists(fnamet_imp2):
        sst2 = jetto_binary_tools.read_binary_file(fnamet_imp2)
    
    sst3 = None
    fname_imp3 = f"{base_path}/jetto.sst3"
    if os.path.exists(fname_imp3):
        sst3 = jetto_binary_tools.read_binary_file(fname_imp3)
    
    jsp  = jetto_binary_tools.read_binary_file(fnamep)
    
    tran = None
    if os.path.exists(ftran):
        try:
            tran = ep.names(ftran)
        except Exception as e:
            print(f"Warning: Could not read TRAN file. Error: {e}")

    # --- PRINT AVAILABLE DATA ---
    print("JST available data:", list(jst.keys()))
    
    if sst1 is not None:
        print("SST1 available data:", list(sst1.keys()))
    
    if sst2 is not None:
        print("SST2 available data:", list(sst2.keys()))
    
    if sst3 is not None:
        print("SST3 available data:", list(sst3.keys()))
    
    print("JSP available data:", list(jsp.keys()))
    
    if tran is not None:
        print("TRAN available signals:", tran)
    
    return jst, sst1, sst2, sst3, jsp, tran

def main():
    if len(sys.argv) > 1:
        run_dir = sys.argv[1] # with this line you can run the script from command line with: python available_data.py <run_directory>
    else:
        run_dir = input("Run directory: ").strip()
    
    try:
        jst, sst1, sst2, sst3, jsp, tran = available_data(run_dir)
        print("Data successfully loaded.")
    except FileNotFoundError as e:
        print(e)

if __name__ == "__main__":
    main()
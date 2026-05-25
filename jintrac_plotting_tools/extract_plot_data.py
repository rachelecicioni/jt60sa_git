"""
R. Cicioni, 25-05-2025

Utility to extract data from JETTO/EDGE2D simulation runs.
Reads binary files (JST, SST1, SST2, SST3, JSP) and TRAN files, return arrays for plotting.

Example usage:
    $ python extract_plot_data.py <run_directory>
Or:
    >>> from extract_plot_data import extract_plot_data
    >>> extract_plot_data("run_sa_nclass2")
    >>> x, y = extract_plot_data('JST', 'PRAD', 'run_sa_nclass2')
"""

import sys
import os
import numpy as np
import jetto_binary_tools
import eproc as ep
from data_available import find_run_directory


def extract_plot_data(tipo_dato, nome_variabile, run_dir, sst_number=None, timestep=-1):
    
    base_path = find_run_directory(run_dir)
    tipo_dato = tipo_dato.upper()
    
    if tipo_dato == 'JST':
        data = jetto_binary_tools.read_binary_file(f"{base_path}/jetto.jst")
        y = np.array(data[nome_variabile]).flatten()
        x = np.array(data['TVEC1']).flatten()
        
    elif tipo_dato == 'JSP':
        data = jetto_binary_tools.read_binary_file(f"{base_path}/jetto.jsp")
        y = np.array(data[nome_variabile])[timestep, :].flatten()
        x = np.array(data['XRHO'])[timestep, :].flatten()
        
    elif tipo_dato.startswith('SST'):
        if sst_number is None:
            sst_number = int(tipo_dato[3:])
        data = jetto_binary_tools.read_binary_file(f"{base_path}/jetto.sst{sst_number}")
        y = np.array(data[nome_variabile]).flatten()
        x = np.array(data['TVEC1']).flatten()
        
    elif tipo_dato == 'TRAN':
        result = ep.time(f"{base_path}/tran", nome_variabile)
        x = np.array(result.xData)
        y = np.array(result.yData)
        
    else:
        raise ValueError(f"Tipo non supportato: {tipo_dato}")
    
    return x, y

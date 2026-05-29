"""
author: R. Cicioni
last updated: 29-05-2025

This code lists the signals available in the TRAN file of a JT-60SA edge2d simulation.
This code can be run in both JDC and EFGW servers, but the path to the TRAN file (=run path) must be updated accordingly (see RUN_DIR variable).
"""
import sys
import os
import matplotlib
matplotlib.use("TkAgg") 
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sys.path.append('/home/lgarzot/Work/Python/JET') #jintrac python tools on JDC server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/python') #jintrac python tools on EFGW server
sys.path.append('/afs/eufus.eu/user/g/g2fjc/jintrac/37.0.0/libs') #eproc tools on EFGW server

import numpy as np
import jetto_binary_tools
import eproc as ep
import matplotlib.pyplot as plt

#RUN_DIR = "/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/dec0525/seq#1" #example run on JDC server
RUN_DIR = "/pfs/work/g2rcicio/edge2d/runs/run_sa_fk_c_WAr_p20_s19_a3_50ms_std" #example run on EFGW server
TRAN_FILE = os.path.join(RUN_DIR, "tran")

#print signals in the TRAN file
print("Available signals in the TRAN file:")
signals = ep.names(TRAN_FILE)

'''
answer=input("Do you want the desciption of the signal? (y/n): ")
if answer.lower() == "y":
    signal_name=input("Insert the name of the signal you want to know more about: ")
    data = ep.data(TRAN_FILE, signal_name)
    units = data.units if hasattr(data, 'units') else 'N/A'
    print(f"Description: {data.desc}")
    print(f"Units: {units}")
'''

import sys
sys.path.append('/home/lgarzot/Work/Python/JET')

from ppf import *
import jetto_binary_tools
import eproc as ep
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.grid'] = True

'''
jdir_list = [
    "run_sa_sg_sett_WAr_p21_s20",
    "run_sa_sg_sett_WAr_p22_s20",
    "run_sa_sg_sett_WAr_p5e22_s20"
]

'''
jdir_list = [
    "run_sa_sg_sett_WAr_p5e22_s19",
    "run_sa_sg_sett_WAr_p5e22_s20",
    "run_sa_sg_sett_WAr_p5e22_s21"
]

base_root = "/common/cmg/jnv7243/edge2d/runs"

# ============================================================
# 1. Potenza radiata totale core
# ============================================================
plt.figure(figsize=(5,3))
for jdir in jdir_list:
    base_path = f"{base_root}/{jdir}"
    jst = jetto_binary_tools.read_binary_file(f"{base_path}/jetto.jst")
    time = np.array(jst["TVEC1"]).flatten()
    prad = np.array(jst["PRAD"]).flatten()
    plt.plot(time, -prad * 1e-6, label=jdir)
plt.ylabel('[MW]')
plt.title(r'$P_{rad,core}$ (total)')
plt.xlabel('time [s]')
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.legend()
plt.show()

# ============================================================
# 2. Potenza radiata Argon
# ============================================================
plt.figure(figsize=(5,3))
for jdir in jdir_list:
    base_path = f"{base_root}/{jdir}"
    sst1 = jetto_binary_tools.read_binary_file(f"{base_path}/jetto.sst1")
    time = np.array(sst1["TVEC1"]).flatten()
    prad = np.array(sst1["PT"]).flatten()
    plt.plot(time, prad * 1e-6)
plt.ylabel('[MW]')
plt.title(r'$P_{rad}$ Argon')
plt.xlabel('time [s]')
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

# ============================================================
# 3. Potenza radiata Tungsten
# ============================================================
plt.figure(figsize=(5,3))
for jdir in jdir_list:
    base_path = f"{base_root}/{jdir}"
    sst2 = jetto_binary_tools.read_binary_file(f"{base_path}/jetto.sst2")
    time = np.array(sst2["TVEC1"]).flatten()
    prad = np.array(sst2["PT"]).flatten()
    plt.plot(time, prad * 1e-6)
plt.ylabel('[MW]')
plt.title(r'$P_{rad}$ Tungsten')
plt.xlabel('time [s]')
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

# ============================================================
# 4. QMAXIT
# ============================================================
plt.figure(figsize=(5,3))
for jdir in jdir_list:
    base_path = f"{base_root}/{jdir}"
    result = ep.time(f"{base_path}/tran", "QMAXIT")
    plt.plot(result.xData, np.array(result.yData) * 1e-6)
plt.xlabel('time (s)')
plt.ylabel('[MW/m²]')
plt.title(r'$Q_{max,IT}$')
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

# ============================================================
# 5. QMAXOT
# ============================================================
plt.figure(figsize=(5,3))
for jdir in jdir_list:
    base_path = f"{base_root}/{jdir}"
    result = ep.time(f"{base_path}/tran", "QMAXOT")
    plt.plot(result.xData, np.array(result.yData) * 1e-6)
plt.xlabel('time (s)')
plt.ylabel('[MW/m²]')
plt.title(r'$Q_{max,OT}$')
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
import jetto_binary_tools
import eproc as ep
import matplotlib.pyplot as plt
import numpy as np
import os

# ============================================================
# PLOT CONFIGURATION
# ============================================================
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.grid'] = True

# ============================================================
# REQUEST RUN NAME FROM TERMINAL
# ============================================================
jdir = input("Insert the run name: ")

base_root = "/pfs/work/g2rcicio/edge2d/runs"
base_path = f"{base_root}/{jdir}"

print(f"Loading run: {jdir}")

# ============================================================
# FILE PATHS
# ============================================================
fname_jst = f"{base_path}/jetto.jst"
fname_sst1 = f"{base_path}/jetto.sst1"
fname_sst2 = f"{base_path}/jetto.sst2"
fname_tran = f"{base_path}/tran"

# ============================================================
# READ FILES
# ============================================================
jst = jetto_binary_tools.read_binary_file(fname_jst)
sst1 = jetto_binary_tools.read_binary_file(fname_sst1)
sst2 = jetto_binary_tools.read_binary_file(fname_sst2)

# ============================================================
# MULTI-PANEL FIGURE (1x3)
# ============================================================
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5), sharex=False)

# ============================================================
# PANEL A) P_rad,TOT
# ============================================================
time_jst = np.array(jst["TVEC1"]).flatten()
prad_jst = np.array(jst["PRAD"]).flatten()
ax1.plot(time_jst, -prad_jst * 1e-6, label="P_rad, TOT", color='black')

time_sst1 = np.array(sst1["TVEC1"]).flatten()
prad_sst1 = np.array(sst1["PT"]).flatten()
ax1.plot(time_sst1, prad_sst1 * 1e-6, label="P_rad, Ar", color='blue')

time_sst2 = np.array(sst2["TVEC1"]).flatten()
prad_sst2 = np.array(sst2["PT"]).flatten()
ax1.plot(time_sst2, prad_sst2 * 1e-6, label="Prad, W", color='red')

ax1.set_xlabel('time [s]')
ax1.set_ylabel('[MW]')
ax1.set_title(r'$P_{rad,core}$')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.4)

# ============================================================
# PANEL B) P_rad,SOL
# ============================================================
result_powsol = ep.time(fname_tran, "POWSOL")
result_zrad1 = ep.time(fname_tran, "ZRAD_1")
result_zrad2 = ep.time(fname_tran, "ZRAD_2")

ax2.plot(result_powsol.xData, np.array(result_powsol.yData) * 1e-6, label="Prad, tot", color='black')
ax2.plot(result_zrad1.xData, np.array(result_zrad1.yData) * 1e-6, label="Prad, Ar", color='blue')
ax2.plot(result_zrad2.xData, np.array(result_zrad2.yData) * 1e-6, label="Prad, W", color='red')

ax2.set_xlabel('time [s]')
ax2.set_ylabel('[MW]')
ax2.set_title(r'$P_{rad,SOL}$')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.4)

# ============================================================
# PANEL C) Power load targets
# ============================================================
result_qmaxot = ep.time(fname_tran, "QMAXOT")
result_qmaxit = ep.time(fname_tran, "QMAXIT")

ax3.plot(result_qmaxot.xData, np.array(result_qmaxot.yData) * 1e-6, label="outer target", color='green')
ax3.plot(result_qmaxit.xData, np.array(result_qmaxit.yData) * 1e-6, label="inner target", color='orange')

ax3.set_xlabel('time [s]')
ax3.set_ylabel('[MW/m²]')
ax3.set_title('Power load divertor targets')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.4)

# ============================================================
# LEGEND AND FORMATTING
# ============================================================
for ax in [ax1, ax2, ax3]:
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)

plt.show()

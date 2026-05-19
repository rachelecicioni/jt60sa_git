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
# MULTI-PANEL FIGURE (1x2)
# ============================================================
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 5), sharex=False)

# ============================================================
# PANEL A) Ar content
# ============================================================
result_ztotal_1 = ep.time(fname_tran, "ZTOTAL_1")
result_zlevsn_1_ar = ep.time(fname_tran, "ZLEVSN_1")

time_sst1 = np.array(sst1["TVEC1"]).flatten()
prad_sst1 = np.array(sst1["NP"]).flatten()

# SOL content: ZTOTAL_1 - ZLEVSN_1
ar_sol_content = np.array(result_ztotal_1.yData) - np.array(result_zlevsn_1_ar.yData)

ax1.plot(result_ztotal_1.xData, np.array(result_ztotal_1.yData), label="Ar total content", color='black')
ax1.plot(result_ztotal_1.xData, ar_sol_content, label="Ar SOL content", color='blue')
ax1.plot(time_sst1, prad_sst1, label="Ar core content", color='red')

ax1.set_xlabel('time [s]')
ax1.set_ylabel('content')
ax1.set_title('Ar content')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.4)

# ============================================================
# PANEL B) W content
# ============================================================
result_ztotal_2 = ep.time(fname_tran, "ZTOTAL_2")
result_zlevsn_2_w = ep.time(fname_tran, "ZLEVSN_2")

time_sst2 = np.array(sst2["TVEC1"]).flatten()
prad_sst2 = np.array(sst2["NP"]).flatten()

# SOL content: ZTOTAL_2 - ZLEVSN_2
w_sol_content = np.array(result_ztotal_2.yData) - np.array(result_zlevsn_2_w.yData)

ax2.plot(result_ztotal_2.xData, np.array(result_ztotal_2.yData), label="W total content", color='black')
ax2.plot(result_ztotal_2.xData, w_sol_content, label="W SOL content", color='green')
ax2.plot(time_sst2, prad_sst2, label="W core content", color='orange')

ax2.set_xlabel('time [s]')
ax2.set_ylabel('content')
ax2.set_title('W content')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.4)

# ============================================================
# LEGEND AND FORMATTING
# ============================================================
for ax in [ax1, ax2]:
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)

plt.show()

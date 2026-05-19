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
# LIST OF SIMULATIONS TO PLOT
# ============================================================
#c-wall ped
'''
jdir_list = [
    "run_sa_fk_c_p20_a3_ArW_fc_chi3e3_d8e2_noelm",
    "run_sa_fk_c_p20_a3_ArW_fc_chi4e3_d8e2_noelm",
    "run_sa_fk_c_p20_a3_ArW_fc_chi6e3_d8e2_noelm",
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_noelm"
]
'''

#puff scan
'''
jdir_list = [
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_noelm",
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_p19",
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_p21"
]
'''

#seeding scan
'''
jdir_list = [
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_noelm",
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_p20_s19",
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_p20_s21"

]
'''

#original pedestal
'''
jdir_list = [
    "run_sa_fk_c_p20_a3_ArW_fc",
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_noelm"
]
'''

#original pedestal D puff scan
jdir_list = ["run_sa_fk_c_p20_a3_ArW_fc_p18",
"run_sa_fk_c_p20_a3_ArW_fc_p19",
"run_sa_fk_c_p20_a3_ArW_fc",
"run_sa_fk_c_p20_a3_ArW_fc_p21",
"run_sa_fk_c_p20_a3_ArW_fc_p22" ]

#original pedestal Ar scan
'''
jdir_list = ["run_sa_fk_c_p20_a3_ArW_fc_p20_s21",
"run_sa_fk_c_p20_a3_ArW_fc_p20_s22"]
'''


base_root = "/common/cmg/jnv7243/edge2d/runs"

# ============================================================
# GLOBAL TIME WINDOW FOR X-LIMITS
# ============================================================
# Set the same x-axis limits for all plots that use an xlim
t_min = 61.0
t_max = 61.1
# ============================================================
# MULTI-PANEL FIGURE (3x3)
# ============================================================
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(18, 12), sharex=False)
(
    ax1, ax2, ax3,
    ax4, ax5, ax6,
    ax7, ax8, ax9
) = axs.flatten()

# ============================================================
# LOOP OVER SIMULATIONS
# ============================================================
for jdir in jdir_list:
    print(f"Upload run: {jdir}")
    base_path = f"{base_root}/{jdir}"

    fname_jst = f"{base_path}/jetto.jst"
    fname_sst1 = f"{base_path}/jetto.sst1"
    fname_sst2 = f"{base_path}/jetto.sst2"
    fname_jsp = f"{base_path}/jetto.jsp"
    fname_tran = f"{base_path}/tran"

    # ============================================================
    # READ FILES
    # ============================================================
    jst = jetto_binary_tools.read_binary_file(fname_jst)
    sst1 = jetto_binary_tools.read_binary_file(fname_sst1)
    sst2 = jetto_binary_tools.read_binary_file(fname_sst2)
    jsp = jetto_binary_tools.read_binary_file(fname_jsp)

    # ============================================================
    # 1. Total radiated power (core)
    # ============================================================
    time = np.array(jst["TVEC1"]).flatten()
    prad = np.array(jst["PRAD"]).flatten()
    ax1.plot(time, -prad * 1e-6, label=jdir)
    ax1.set_ylabel('[MW]')
    ax1.set_title(r'$P_{rad,core}$ (total)')
    ax1.set_xlabel('time [s]')

    # ============================================================
    # 2. Impurity 1 power (Argon)
    # ============================================================
    time_imp1 = np.array(sst1["TVEC1"]).flatten()
    prad_imp1 = np.array(sst1["PT"]).flatten()
    ax4.plot(time_imp1, prad_imp1 * 1e-6, label=jdir)
    ax4.set_ylabel('[MW]')
    ax4.set_title(r'$P_{rad}$ Argon')
    ax4.set_xlabel('time [s]')

    # ============================================================
    # 3. Impurity 2 power (Tungsten)
    # ============================================================
    time_imp2 = np.array(sst2["TVEC1"]).flatten()
    prad_imp2 = np.array(sst2["PT"]).flatten()
    ax7.plot(time_imp2, prad_imp2 * 1e-6, label=jdir)
    ax7.set_ylabel('[MW]')
    ax7.set_title(r'$P_{rad}$ Tungsten')
    ax7.set_xlabel('time [s]')

    # (Removed TOB/SEP time-series panels - now using 3x3 layout)

    # ============================================================
    # 4. Radial profile of Zeff (last timestep)
    # ============================================================
    idx_last = -1
    rho = np.array(jsp["XRHO"][idx_last, :]).flatten()
    zeff = np.array(jsp["ZEFF"][idx_last, :]).flatten()
    ax2.plot(rho, zeff, label=jdir)
    ax2.set_xlabel(r'$\rho$')
    ax2.set_ylabel(r'$Z_{\mathrm{eff}}$')
    ax2.set_title(r'$Z_{\mathrm{eff}}$ last timestep')

    # ============================================================
    # 5. QMAXIT
    # ============================================================
    result = ep.time(fname_tran, "QMAXIT")
    ax5.plot(result.xData, np.array(result.yData) * 1e-6, label=jdir)
    ax5.set_xlabel('time (s)')
    ax5.set_ylabel('[MW/m²]')
    ax5.set_xlim(t_min, t_max)
    ax5.axhline(y=10, color="red", linewidth=1.0, linestyle='-')
    ax5.set_title(r'$Q_{max,IT}$')

    # ============================================================
    # 6. QMAXOT
    # ============================================================
    result = ep.time(fname_tran, "QMAXOT")
    ax8.plot(result.xData, np.array(result.yData) * 1e-6, label=jdir)
    ax8.set_xlabel('time (s)')
    ax8.set_ylabel('[MW/m²]')
    ax8.set_xlim(t_min, t_max)
    ax8.axhline(y=10, color="red", linewidth=1.0, linestyle='-')
    ax8.set_title(r'$Q_{max,OT}$')

    # ============================================================
    # 7. Radial profile of Te (last timestep)
    # ============================================================
    te = np.array(jsp["TE"][idx_last, :]).flatten()
    ax3.plot(rho, te*1e-3, label=jdir)
    ax3.axhline(y=0.150, color="red", linewidth=1.0, linestyle='--')
    ax3.set_xlabel(r'$\rho$')
    ax3.set_ylabel('[keV]')
    ax3.set_title(r'$T_e$ last timestep')

    # ============================================================
    # 8. Radial profile of Ti (last timestep)
    # ============================================================
    ti = np.array(jsp["TI"][idx_last, :]).flatten()
    ax6.plot(rho, ti*1e-3, label=jdir)
    ax6.axhline(y=0.230, color="red", linewidth=1.0, linestyle='--')
    ax6.set_xlabel(r'$\rho$')
    ax6.set_ylabel('[keV]')
    ax6.set_title(r'$T_i$ last timestep')

    # ============================================================
    # 9. Radial profile of ne (last timestep)
    # ============================================================
    ne = np.array(jsp["NE"][idx_last, :]).flatten()
    ax9.plot(rho, ne, label=jdir)
    ax9.axhline(y=3e19, color="red", linewidth=1.0, linestyle='--')
    ax9.set_xlabel(r'$\rho$')
    ax9.set_ylabel('[m-3]')
    ax9.set_title(r'$n_e$ last timestep')

# ============================================================
# LEGEND AND FORMATTING
# ============================================================
# ============================================================
# OVERLAY C-WALL SCENARIO PROFILES
# ============================================================
cwall_jsp_path = "/home/lgarzot/cmg/catalog/jetto/jet/70000/oct2115/seq#1/jetto.jsp"
if os.path.exists(cwall_jsp_path):
    jsp_c = jetto_binary_tools.read_binary_file(cwall_jsp_path)
    idx_last = -1
    rho_c = np.array(jsp_c["XRHO"][idx_last, :]).flatten()
    te_c = np.array(jsp_c["TE"][idx_last, :]).flatten()
    ti_c = np.array(jsp_c["TI"][idx_last, :]).flatten()
    ne_c = np.array(jsp_c["NE"][idx_last, :]).flatten()
    # Replace last-column plots with C-wall profiles (place into new 3x3 axes)
    ax3.plot(rho_c, te_c*1e-3, color='k', linewidth=1.5, label='C-wall scenario')
    ax3.axhline(y=0.150, color="red", linewidth=1.0, linestyle='--')
    ax3.set_xlabel(r'$\rho$')
    ax3.set_ylabel('[keV]')
    ax3.set_title(r'$T_e$ last timestep (C-wall)')

    ax6.plot(rho_c, ti_c*1e-3, color='k', linewidth=1.5, label='C-wall scenario')
    ax6.axhline(y=0.230, color="red", linewidth=1.0, linestyle='--')
    ax6.set_xlabel(r'$\rho$')
    ax6.set_ylabel('[keV]')
    ax6.set_title(r'$T_i$ last timestep (C-wall)')

    ax9.plot(rho_c, ne_c, color='k', linewidth=1.5, label='C-wall scenario')
    ax9.axhline(y=3e19, color="red", linewidth=1.0, linestyle='--')
    ax9.set_xlabel(r'$\rho$')
    ax9.set_ylabel('[m-3]')
    ax9.set_title(r'$n_e$ last timestep (C-wall)')
else:
    print(f"C-wall jsp not found: {cwall_jsp_path}")
for ax in [
    ax1, ax2, ax3,
    ax4, ax5, ax6,
    ax7, ax8, ax9
]:
    ax1.legend(fontsize=7)
    ax.grid(True, alpha=0.4)

plt.show()

plot_dir = "/home/jnv7243/JT-60SA/plot"
os.makedirs(plot_dir, exist_ok=True)
out_path = os.path.join(plot_dir, "puff_seeding_scan.png")
plt.tight_layout()
fig.savefig(out_path, dpi=200)
print(f"Saved figure to {out_path}")
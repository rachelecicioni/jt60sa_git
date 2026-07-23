import sys
sys.path.append('/home/lgarzot/Work/Python/JET')

from ppf import *
import jetto_binary_tools
import eproc as ep
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# CONFIGURAZIONE GRAFICI
# ============================================================
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.grid'] = True

# ============================================================
# LISTA DELLE SIMULAZIONI DA PLOTTARE
# ============================================================
'''
jdir_list = [
    "run_sa_fk_c_p20_a3_ArW_fc"
]
'''
jdir_list = [
    "run_sa_fk_c_p20_a3_ArW_fc_chi8e3_d8e2_noelm",
    "run_sa_fk_c_p20_a3_ArW_fc"
]

base_root = "/common/cmg/jnv7243/edge2d/runs"

# ============================================================
# FIGURA MULTI-PANNELLO
# ============================================================
fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(18, 10), sharex=False)
(
    ax1, ax2, ax3, ax4,
    ax5, ax6, ax7, ax8,
    ax9, ax10, ax11, ax12
) = axs.flatten()

# ============================================================
# CICLO SULLE SIMULAZIONI
# ============================================================
for jdir in jdir_list:
    print(f"Caricando simulazione: {jdir}")
    base_path = f"{base_root}/{jdir}"

    fname_jst = f"{base_path}/jetto.jst"
    fname_sst1 = f"{base_path}/jetto.sst1"
    fname_sst2 = f"{base_path}/jetto.sst2"
    fname_jsp = f"{base_path}/jetto.jsp"
    fname_tran = f"{base_path}/tran"

    #print("\n--- Variabili disponibili nel file TRAN ---")
    #tran_vars = ep.names(fname_tran)

    #for v in tran_vars:
        #print(v)

    # ============================================================
    # LETTURA FILES
    # ============================================================
    jst = jetto_binary_tools.read_binary_file(fname_jst)
    sst1 = jetto_binary_tools.read_binary_file(fname_sst1)
    sst2 = jetto_binary_tools.read_binary_file(fname_sst2)
    jsp = jetto_binary_tools.read_binary_file(fname_jsp)

    # ============================================================
    # 1. Potenza radiata totale core
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
    ax5.plot(time_imp1, prad_imp1 * 1e-6, label=jdir)
    ax5.set_ylabel('[MW]')
    ax5.set_title(r'$P_{rad}$ Argon')
    ax5.set_xlabel('time [s]')

    # ============================================================
    # 3. Impurity 2 power (Tungsten)
    # ============================================================
    time_imp2 = np.array(sst2["TVEC1"]).flatten()
    prad_imp2 = np.array(sst2["PT"]).flatten()
    ax9.plot(time_imp2, prad_imp2 * 1e-6, label=jdir)
    ax9.set_ylabel('[MW]')
    ax9.set_title(r'$P_{rad}$ Tungsten')
    ax9.set_xlabel('time [s]')

    # ============================================================
    # 4. Densità elettronica TOB / SEP
    # ============================================================
    ne_tob = np.array(jst["NEBA"]).flatten()
    ne_sep = np.array(jst["NEBO"]).flatten()
    ax2.plot(time, ne_tob, label=fr'$n_{{e,TOB}}$ {jdir}')
    ax2.plot(time, ne_sep, linestyle="--", label=fr'$n_{{e,SEP}}$ {jdir}')
    ax2.axhline(y=3e19, color="red", linewidth=1.0, linestyle='--')
    ax2.set_ylabel(r'$n_e$ [m$^{-3}$]')
    ax2.set_title(r'$n_{e,TOB}$ / $n_{e,SEP}$')
    ax2.set_xlabel('time [s]')

    # ============================================================
    # 5. Temperatura elettronica TOB / SEP
    # ============================================================
    te_tob = np.array(jst["TEBA"]).flatten() / 1e3
    te_sep = np.array(jst["TEBO"]).flatten() / 1e3
    ax6.plot(time, te_tob, label=fr'$T_{{e,TOB}}$ {jdir}')
    ax6.plot(time, te_sep, linestyle="--", label=fr'$T_{{e,SEP}}$ {jdir}')
    ax6.axhline(y=0.150, color="red", linewidth=1.0, linestyle='--')
    ax6.set_ylabel('[keV]')
    ax6.set_title(r'$T_{e,TOB}$ / $T_{e,SEP}$')
    ax6.set_xlabel('time [s]')

    # ============================================================
    # 6. Temperatura ionica TOB / SEP
    # ============================================================
    ti_tob = np.array(jst["TIBA"]).flatten() / 1e3
    ti_sep = np.array(jst["TIBO"]).flatten() / 1e3
    ax10.plot(time, ti_tob, label=fr'$T_{{i,TOB}}$ {jdir}')
    ax10.plot(time, ti_sep, linestyle="--", label=fr'$T_{{i,SEP}}$ {jdir}')
    ax10.axhline(y=0.230, color="red", linewidth=1.0, linestyle='--')
    ax10.set_ylabel('[keV]')
    ax10.set_title(r'$T_{i,TOB}$ / $T_{i,SEP}$')
    ax10.set_xlabel('time [s]')

    # ============================================================
    # 7. Profilo radiale di Zeff (ultimo timestep)
    # ============================================================
    idx_last = -1
    rho = np.array(jsp["XRHO"][idx_last, :]).flatten()
    zeff = np.array(jsp["ZEFF"][idx_last, :]).flatten()
    ax3.plot(rho, zeff, label=jdir)
    ax3.set_xlabel(r'$\rho$')
    ax3.set_ylabel(r'$Z_{\mathrm{eff}}$')
    ax3.set_title(r'$Z_{\mathrm{eff}}$ last timestep')

    # ============================================================
    # 8. QMAXIT
    # ============================================================
    result = ep.time(fname_tran, "QMAXIT")
    ax7.plot(result.xData, np.array(result.yData) * 1e-6, label=jdir)
    ax7.set_xlabel('time (s)')
    ax7.set_ylabel('[MW/m²]')
    ax7.set_xlim(61.0, 61.1)
    ax7.set_title(r'$Q_{max,IT}$')

    # ============================================================
    # 9. QMAXOT
    # ============================================================
    result = ep.time(fname_tran, "QMAXOT")
    ax11.plot(result.xData, np.array(result.yData) * 1e-6, label=jdir)
    ax11.set_xlabel('time (s)')
    ax11.set_ylabel('[MW/m²]')
    ax11.set_xlim(61.0, 61.1)
    ax11.set_title(r'$Q_{max,OT}$')

    # ============================================================
    # 10. Profilo radiale di Te (ultimo timestep)
    # ============================================================
    te = np.array(jsp["TE"][idx_last, :]).flatten()
    ax4.plot(rho, te*1e-3, label=jdir)
    ax4.set_xlabel(r'$\rho$')
    ax4.set_ylabel('[keV]')
    ax4.set_title(r'$T_e$ last timestep')

    # ============================================================
    # 11. Profilo radiale di Ti (ultimo timestep)
    # ============================================================
    ti = np.array(jsp["TI"][idx_last, :]).flatten()
    ax8.plot(rho, ti*1e-3, label=jdir)
    ax8.set_xlabel(r'$\rho$')
    ax8.set_ylabel('[keV]')
    ax8.set_title(r'$T_i$ last timestep')

    # ============================================================
    # 12. Profilo radiale di ne (ultimo timestep)
    # ============================================================
    ne = np.array(jsp["NE"][idx_last, :]).flatten()
    ax12.plot(rho, ne, label=jdir)
    ax12.set_xlabel(r'$\rho$')
    ax12.set_ylabel('[m-3]')
    ax12.set_title(r'$n_e$ last timestep')

# ============================================================
# LEGENDA E FORMATTING
# ============================================================
for ax in [
    ax1, ax2, ax3, ax4,
    ax5, ax6, ax7, ax8,
    ax9, ax10, ax11, ax12
]:
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

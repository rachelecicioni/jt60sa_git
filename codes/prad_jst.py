#!/usr/bin/env python

import sys
sys.path.append('/home/lgarzot/Work/Python/JET')

from ppf import *
import jetto_binary_tools

import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib import cm 
import numpy as np

plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['black'])
plt.rcParams['axes.grid'] = True  # griglia attiva di default

# --- INPUT ---
code_type = input("Simulation code (jetto / edge2d): ").strip().lower()
jdir = input("Run directory: ").strip()

# --- Controllo input ---
if code_type not in ["jetto", "edge2d"]:
    raise ValueError(" Inserisci 'jetto' oppure 'edge2d' come tipo di simulazione.")

# --- FILE PATHS ---
base_path = f"/common/cmg/jnv7243/{code_type}/runs/{jdir}"

fnamet       = f"{base_path}/jetto.jst"   # file JST
fnamet_imp1  = f"{base_path}/jetto.sst1"  # impurity 1
fnamet_imp2  = f"{base_path}/jetto.sst2"  # impurity 2
fnamep       = f"{base_path}/jetto.jsp_savs2"   # file JSP


# --- LETTURA FILE ---
jst = jetto_binary_tools.read_binary_file(fnamet)
sst1 = jetto_binary_tools.read_binary_file(fnamet_imp1)
sst2 = jetto_binary_tools.read_binary_file(fnamet_imp2)
jsp  = jetto_binary_tools.read_binary_file(fnamep)



# Tipo di file
print("Tipo JST:", type(jst))
print("Tipo JSP:", type(jsp))
print("Tipo SST1:", type(sst1))
print("Tipo SST2:", type(sst2))

# Chiavi disponibili
print("\nChiavi in JST (time traces):")
print(list(jst.keys()))

print("\nChiavi in JSP (profiles):")
print(list(jsp.keys()))

print("\nChiavi in SST1 (time traces):")
print(list(sst1.keys()))

print("\nChiavi in SST2 (time traces):")
print(list(sst2.keys()))


# --- FIGURA ---
fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(10, 12), sharex=False)
(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10) = axs.flatten()

# ============================================================
# 1. Potenza radiata totale (core)
# ============================================================
prad = jst['PRAD'].flatten()
time = jst['TVEC1'].flatten()
ax1.plot(time, -prad * 1e-6)
ax1.set_ylabel('[MW]')
ax1.set_title('JST/PRAD — P_rad core (total)')
ax1.set_xlabel('time [s]')

# ============================================================
# 2. Bremsstrahlung
# ============================================================
prad = jst['PBRE'].flatten()
ax2.plot(time, -prad * 1e-6)
ax2.set_ylabel('[MW]')
ax2.set_title('JST/PBRE — Bremsstrahlung')
ax2.set_xlabel('time [s]')

# ============================================================
# 3. Line emission
# ============================================================
prad = jst['PLIN'].flatten()
ax3.plot(time, -prad * 1e-6)
ax3.set_ylabel('[MW]')
ax3.set_title('JST/PLIN — Line emission (impurities)')
ax3.set_xlabel('time [s]')

# ============================================================
# 4. Synchrotron radiation
# ============================================================
prad = jst['PSYR'].flatten()
ax4.plot(time, -prad * 1e-6)
ax4.set_ylabel('[MW]')
ax4.set_title('JST/PSYR — Synchrotron radiation')
ax4.set_xlabel('time [s]')

# ============================================================
# 5. Impurity 1 power
# ============================================================
prad = sst1['PT'].flatten()
time_imp1 = sst1['TVEC1'].flatten()
ax5.plot(time_imp1, prad * 1e-6)
ax5.set_ylabel('[MW]')
ax5.set_title('SST1/PT — Impurity 1 P_rad')
ax5.set_xlabel('time [s]')

# ============================================================
# 6. Impurity 2 power
# ============================================================
prad = sst2['PT'].flatten()
time_imp2 = sst2['TVEC1'].flatten()
ax6.plot(time_imp2, prad * 1e-6)
ax6.set_ylabel('[MW]')
ax6.set_title('SST2/PT — Impurity 2 P_rad')
ax6.set_xlabel('time [s]')

# ============================================================
# 7. Densità elettronica TOB / SEP
# ============================================================
y_t = jst['TVEC1'].flatten()
ne_tob = jst['NEBA'].flatten()
ne_sep = jst['NEBO'].flatten()
ax7.plot(y_t, ne_tob, label=r'$n_{e,TOB}$', color='blue')
ax7.plot(y_t, ne_sep, label=r'$n_{e,SEP}$', color='black')
ax7.set_ylabel(r'$n_e$ [m$^{-3}$]')
ax7.set_title(r'$n_{e,TOB}$ / $n_{e,SEP}$')
ax7.legend()
ax7.set_xlabel('time [s]')

# ============================================================
# 8. Temperatura elettronica TOB / SEP (in keV)
# ============================================================
te_tob = jst['TEBA'].flatten() / 1e3
te_sep = jst['TEBO'].flatten() / 1e3
ax8.plot(y_t, te_tob, label=r'$T_{e,TOB}$', color='blue')
ax8.plot(y_t, te_sep, label=r'$T_{e,SEP}$', color='black')
ax8.set_ylabel('[keV]')
ax8.set_title(r'$T_{e,TOB}$ / $T_{e,SEP}$')
ax8.legend()
ax8.set_xlabel('time [s]')

# ============================================================
# 9. Temperatura ionica TOB / SEP (in keV)
# ============================================================
ti_tob = jst['TIBA'].flatten() / 1e3
ti_sep = jst['TIBO'].flatten() / 1e3
ax9.plot(y_t, ti_tob, label=r'$T_{i,TOB}$', color='blue')
ax9.plot(y_t, ti_sep, label=r'$T_{i,SEP}$', color='black')
ax9.set_ylabel('[keV]')
ax9.set_title(r'$T_{i,TOB}$ / $T_{i,SEP}$')
ax9.legend()
ax9.set_xlabel('time [s]')

# ============================================================
# 10. Profilo radiale di Zeff (ultimo istante)
# ===========================================================
t_jsp = jsp['TIME'].flatten()
idx_last = -1  # ultimo tempo
rho = jsp['XRHO'][idx_last, :].flatten()
zeff = jsp['ZEFF'][idx_last, :].flatten()
ax10.plot(rho, zeff, color='purple', linewidth=1.5)
ax10.set_xlabel(r'$\rho$')
ax10.set_ylabel(r'$Z_{\mathrm{eff}}$')
ax10.set_title(rf'$Z_{{\mathrm{{eff}}}}$ profile — last timestep (t = {t_jsp[idx_last]:.3f} s)')

# ============================================================
# FORMATTAZIONE GENERALE
# ============================================================
for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]:
    ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

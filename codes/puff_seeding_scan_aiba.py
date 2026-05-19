import sys
import shutil
sys.path.append('/home/lgarzot/Work/Python/JET')

from ppf import *
import jetto_binary_tools
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.grid'] = True

# ============================================================
# Inserimento da terminale del run da analizzare
# ============================================================
run_name = input("Inserisci il nome della cartella del run: ")

base_root = "/common/cmg/jnv7243/edge2d/runs"
fname_jsp = f"{base_root}/{run_name}/jetto.jsp_savs"

# ============================================================
# Lettura del file .jsp_savs
# ============================================================
jsp = jetto_binary_tools.read_binary_file(fname_jsp)

# ============================================================
# Array dei tempi disponibili
# ============================================================
time_array = np.array(jsp["TIME"]).flatten()
print("\nElenco degli istanti temporali con indice:")
for idx, t in enumerate(time_array):
    print(f"{idx:3d}: {t:.6f} s")

# ============================================================
# Inserimento da terminale dell'indice da plottare
# ============================================================
selected_idx = int(input("\nInserisci l'indice del timestep da plottare: "))
print(f"Selezionato timestep: index {selected_idx}, time {time_array[selected_idx]:.6f} s")

# ============================================================
# Preparazione figura per i profili
# ============================================================
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
ax_zeff, ax_te, ax_ti, ax_ne = axs.flatten()

# Coordinate radiali
psi_norm = np.array(jsp["XRHO"][selected_idx, :]).flatten()

# ============================================================
# ZEFF
# ============================================================
zeff = np.array(jsp["ZEFF"][selected_idx, :]).flatten()
ax_zeff.plot(psi_norm, zeff, label=run_name)
ax_zeff.set_xlabel('psi_norm')
ax_zeff.set_ylabel('Z_eff')
ax_zeff.set_title('Z_eff')

# ============================================================
# TE
# ============================================================
te = np.array(jsp["TE"][selected_idx, :]).flatten() * 1e-3  # keV
ax_te.plot(psi_norm, te, label=run_name)
ax_te.set_xlabel('psi_norm')
ax_te.set_ylabel('T_e [keV]')
ax_te.set_title('T_e')

# ============================================================
# TI
# ============================================================
ti = np.array(jsp["TI"][selected_idx, :]).flatten() * 1e-3  # keV
ax_ti.plot(psi_norm, ti, label=run_name)
ax_ti.set_xlabel('psi_norm')
ax_ti.set_ylabel('T_i [keV]')
ax_ti.set_title('T_i')

# ============================================================
# NE
# ============================================================
ne = np.array(jsp["NE"][selected_idx, :]).flatten()  # m^-3
ax_ne.plot(psi_norm, ne, label=run_name)
ax_ne.set_xlabel('psi_norm')
ax_ne.set_ylabel('n_e [m^-3]')
ax_ne.set_title('n_e')

# ============================================================
# Legende e layout
# ============================================================
for ax in [ax_zeff, ax_te, ax_ti, ax_ne]:
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

# ============================================================
# Creazione sottocartella per il run
# ============================================================
profiles_dir = "/home/jnv7243/JT-60SA/codes/profiles"
run_save_dir = os.path.join(profiles_dir, run_name)
os.makedirs(run_save_dir, exist_ok=True)

# ============================================================
# Salvataggio dei profili in file .txt
# ============================================================
outname = os.path.join(run_save_dir, f"profile_{run_name}.txt")

with open(outname, "w") as f:
    f.write("# psi_norm     Z_eff     n_e[m^-3]     T_e[keV]     T_i[keV]\n")
    for ps, z, n, te_i, ti_i in zip(psi_norm, zeff, ne, te, ti):
        f.write(f"{ps:12.6e}  {z:10.6f}  {n:12.6e}  {te_i:10.6e}  {ti_i:10.6e}\n")

print(f"\nProfilo salvato in: {outname}")

# ============================================================
# Copia del file "jetto.eqdsk_out" nella sottocartella del run
# ============================================================
src_eqdsk = os.path.join(base_root, run_name, "jetto.eqdsk_out")
dst_eqdsk = os.path.join(run_save_dir, "jetto.eqdsk_out")

if os.path.exists(src_eqdsk):
    shutil.copy(src_eqdsk, dst_eqdsk)
    print(f"File 'jetto.eqdsk_out' copiato in: {dst_eqdsk}")
else:
    print(f"Attenzione: file 'jetto.eqdsk_out' non trovato in {src_eqdsk}")
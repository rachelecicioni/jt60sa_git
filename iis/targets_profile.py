import matplotlib.pyplot as plt
import numpy as np
from pytran import Tran

# === Percorso al file caricato ===
filename = '/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/jun2325/seq#2/tran'

# === Caricamento dati ===
t = Tran(filename)

# === Tempo desiderato (in secondi) ===
t_target = 64055  # <-- cambia a piacere

# === Trova l'indice temporale più vicino ===
t_array = np.array(t.tvec)
j = np.abs(t_array - t_target).argmin()
print(f"Tempo richiesto: {t_target:.5f} s → più vicino: {t_array[j]:.5f} s (j = {j})")

# === Outer target ===
jtarget = 2  # di solito 2 = outer target, 1 = inner target
cell_ids = t.korxy[1:-1, jtarget]
valid = cell_ids > 0
cell_ids = cell_ids[valid] - 1  # Fortran → Python

# === Coordinate lungo target (opzionale: rmesh[cell_ids]) ===
s = np.arange(len(cell_ids))

# === Campo fisico: potenza irradiata ===
power = t.den


if power.ndim == 2:
    power = power[j]  # usa il frame temporale selezionato

pvals = power[cell_ids]

# === Plot ===
plt.figure(figsize=(8, 5))
plt.plot(s, pvals, '-o', color='darkblue')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r'cell grid[]', fontsize=14)
plt.ylabel(r'$n_i (m^{-3})$', fontsize=14)
plt.title(r'$n_i$ along outer target')
plt.grid(True)
plt.tight_layout()
plt.show()


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from jetto_tools.binary import *
import numpy as np
from pytran import Tran #questa è la libreria ncessaria per leggere i file tran

# === Percorso al file caricato ===
filename = '/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/jun2325/seq#2/tran'

# === Caricamento dati ===
t = Tran(filename)

# === Tempo desiderato (secondi) ===
t_target = 63850 # <-- cambia questo valore a piacere

# === Trova l'indice temporale più vicino ===
t_array = np.array(t.tvec)
j = np.abs(t_array - t_target).argmin()
print(f"Tempo richiesto: {t_target:.5f} s → tempo più vicino disponibile: {t_array[j]:.5f} s (indice j = {j})")

korxy = t.korxy[1:-1, 1:-1]
cell_indices = korxy.flatten()
valid = cell_indices > 0 #filtro boolenao per selezionare solo le celle con indice >0

# Correggi per 1-based indexing (Fortran style) questo perchè in Forstran gli indici partono da 1 mentre in python da 0
cell_indices = cell_indices[valid] - 1

r_all = t.rmesh
z_all = t.zmesh

r = r_all[cell_indices] #tutte le posizioni r dei centri delle celle
z = z_all[cell_indices] #tutte le posizioni z dei centri delle celle


# === Scegli campo fisico da plottare ===

values = t.denz01  # potenza irradiata per cella


if values.ndim == 2:
    values = values[j]  # seleziona il frame temporale più vicino a t_target

# === Valori da colorare ===
data = values[cell_indices]

# === Plot ===
plt.figure(figsize=(8, 6))
tri = plt.tricontourf(r, z, data, levels=100, cmap='plasma')

plt.xlabel("R [m]")
plt.ylabel("Z [m]")
plt.title(f"Potenza Irradiata – t = {t_array[j]:.5f} s")
plt.colorbar(tri, label='[W/m³]')
plt.axis('equal')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()
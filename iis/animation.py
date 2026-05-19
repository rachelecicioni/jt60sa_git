import matplotlib.pyplot as plt
import numpy as np
import time  # per rallentare tra i frame
from pytran import Tran

# === Percorso al file caricato ===
filename = '/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/jun2325/seq#2/tran'

# === Caricamento dati ===
t = Tran(filename)

# === Outer target (di solito j=2) ===
jtarget = 2
korxy_col = t.korxy[1:-1, jtarget]
valid = korxy_col > 0
cell_ids = korxy_col[valid] - 1  # Fortran → Python

# === Coordinate lungo il target (indice o fisica)
s = np.arange(len(cell_ids))  # oppure t.rmesh[cell_ids]

# === Campo fisico: potenza irradiata ===
try:
    power = t.sqehrad
except AttributeError:
    power = t.sqezr_1

# === Verifica che abbia dipendenza temporale ===
if power.ndim != 2 or power.shape[0] != len(t.tvec):
    raise ValueError("Il campo selezionato non ha evoluzione temporale compatibile.")

# === Loop sugli istanti temporali ===
for j, tval in enumerate(t.tvec):
    pvals = power[j][cell_ids]

    plt.figure(figsize=(8, 5))
    plt.plot(s, pvals, '-o')
    plt.xlabel('Indice cella lungo target')
    plt.ylabel('P rad [W/m³]')
    plt.title(f'Potenza Irradiata – Outer Target\nTempo = {tval:.5f} s')
    plt.grid(True)
    plt.tight_layout()
    plt.pause(0.5)  # mostra il plot per mezzo secondo
    plt.clf()       # pulisce per il prossimo frame

plt.close()

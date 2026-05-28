#!/usr/bin/env python3
"""
Estrai profilo LOS dal file tran JINTRAC CORRETTAMENTE
Usa ring() e row() per ricostruire la mesh non-rettangolare
"""

import eproc as ep
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# ============================================================
# CONFIGURAZIONE
# ============================================================
RUN_DIR = "/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/dec0525/seq#1"
TRAN_FILE = os.path.join(RUN_DIR, "tran")

SIGNAL_NAME = "TEVE"

# Parametri LOS (da estremi a estremi)
LOS_X1, LOS_Y1 = 4.5, 0.0    # Punto 1
LOS_X2, LOS_Y2 = 3.0, 2.0    # Punto 2

LOS_ORIGIN = np.array([LOS_X1, LOS_Y1])
LOS_DIRECTION = np.array([LOS_X2 - LOS_X1, LOS_Y2 - LOS_Y1])

print("=" * 80)
print(f"LOS PROFILE EXTRACTION")
print("=" * 80)
print(f"Punto 1: ({LOS_X1}, {LOS_Y1})")
print(f"Punto 2: ({LOS_X2}, {LOS_Y2})")
print(f"Origine: {LOS_ORIGIN}")
print(f"Direzione: {LOS_DIRECTION}")

# ============================================================
# CARICA GEOMETRIA
# ============================================================

geom = ep.geom(TRAN_FILE)
nRings = int(geom.nRings)
nRows = int(geom.nRows)

print(f"\nGeometria:")
print(f"  nRings: {nRings}")
print(f"  nRows: {nRows}")
print(f"  firstCoreRow: {geom.firstCoreRow}")
print(f"  lastCoreRow: {geom.lastCoreRow}")
print(f"  firstOpenRing: {geom.firstOpenRing}")
print(f"  wallRing: {geom.wallRing}")

# ============================================================
# RICOSTRUISCI MESH USANDO ring() E row()
# ============================================================

print(f"\nRicostruzione mesh usando ring() e row()...")

# Estrai tutti i ring
nodes = []
quads = []
values = []

ring_names = []
for i in range(nRings):
    # Nomi dei ring: S01, S02, ..., oppure S1, S2, ... a seconda del numero
    if i < geom.firstOpenRing:
        ring_name = f"S{i:02d}" if i >= 10 else f"S{i}"
    else:
        ring_name = f"S{i+1}"
    
    ring_names.append(ring_name)

# Prova a caricare rings e capire la struttura
print(f"\nCaricamento ring...")

ring_data_list = []
for i, ring_name in enumerate(ring_names):
    try:
        ring = ep.ring(TRAN_FILE, SIGNAL_NAME, ring_name)
        ring_data = {
            'name': ring_name,
            'ring_idx': i,
            'nPts': ring.nPts,
            'xData': np.array(ring.xData),  # coord poloidale
            'yData': np.array(ring.yData),  # valori del segnale
        }
        ring_data_list.append(ring_data)
        if i < 5:  # Stampa primi 5
            print(f"  Ring {ring_name} (idx {i}): {ring.nPts} punti")
    except Exception as e:
        if i < 5:
            print(f"  Ring {ring_name} (idx {i}): ERROR - {e}")
        continue

print(f"\nRing caricati: {len(ring_data_list)}")

# Estrai Row per ottenere coordinate R
print(f"\nCaricamento row...")

row_names = ["IT", "OT", "1", "2", "10", "20", "30", "40", "50"]
row_data_list = []

for row_name in row_names:
    try:
        row = ep.row(TRAN_FILE, SIGNAL_NAME, row_name)
        row_data = {
            'name': row_name,
            'nPts': row.nPts,
            'xData': np.array(row.xData),  # coord radiale (offset poloidale)
            'yData': np.array(row.yData),  # valori del segnale
        }
        row_data_list.append(row_data)
        print(f"  Row {row_name}: {row.nPts} punti, range coord={row.xData.min():.4f} to {row.xData.max():.4f}")
    except:
        continue

# ============================================================
# RICOSTRUISCI LA MESH POLOIDALE (da un ring)
# ============================================================

if len(ring_data_list) > 0:
    # Usa il primo ring come riferimento
    ref_ring = ring_data_list[0]
    
    print(f"\nRicostruzione mesh da Ring {ref_ring['name']}...")
    
    # Le coordinate X di un ring sono coordinate poloidali
    # Abbiamo bisogno di R e Z di questi punti
    
    # Per ora, usiamo direttamente i dati di ring e row
    # Ring = sezione radiale (costante in poloidale, varia in radiale)
    # Row = sezione poloidale (costante in radiale, varia in poloidale)
    
    print(f"  {ref_ring['name']}: {ref_ring['nPts']} punti poloidali")
    
    if len(row_data_list) > 0:
        ref_row = row_data_list[0]
        print(f"  {ref_row['name']}: {ref_row['nPts']} punti radiali")

# ============================================================
# APPROCCIO ALTERNATIVO: Usa RMESH E ZMESH DIRETTAMENTE
# ============================================================

print(f"\nUsando RMESH e ZMESH direttamente...")

rmesh = ep.data(TRAN_FILE, "RMESH")
zmesh = ep.data(TRAN_FILE, "ZMESH")
data = ep.data(TRAN_FILE, SIGNAL_NAME)

R_all = np.array(rmesh.data, dtype=float)
Z_all = np.array(zmesh.data, dtype=float)
signal_all = np.array(data.data, dtype=float)

# Scarta i punti nulli
valid_idx = (R_all != 0) | (Z_all != 0) | (signal_all != 0)
R = R_all[valid_idx]
Z = Z_all[valid_idx]
signal = signal_all[valid_idx]

print(f"Punti validi (non nulli): {len(R)} / {len(R_all)}")
print(f"R range: {R.min():.4f} - {R.max():.4f}")
print(f"Z range: {Z.min():.4f} - {Z.max():.4f}")
print(f"{SIGNAL_NAME} range: {signal.min():.4f} - {signal.max():.4f}")

# ============================================================
# ALGORITMO RAY-QUAD
# ============================================================

def cross2d(a, b):
    """2D cross product"""
    return a[0] * b[1] - a[1] * b[0]

def line_segment_intersect_t(origin, direction, p1, p2):
    """Intersect infinite line with segment"""
    d = direction
    v = p2 - p1
    denom = cross2d(d, v)
    if abs(denom) < 1e-12:
        return None, None
    w = p1 - origin
    t = cross2d(w, v) / denom
    s = cross2d(w, d) / denom
    return t, s

def point_in_los(point, origin, direction, tolerance=0.1):
    """Check if point is close to the ray"""
    # Distance from point to line
    # d = |((p - origin) × direction)| / |direction|
    p_rel = point - origin
    cross = cross2d(p_rel, direction)
    dist = abs(cross) / (np.linalg.norm(direction) + 1e-10)
    return dist < tolerance

# ============================================================
# ESTRAI PROFILO SEMPLICE (basato su distanza dal ray)
# ============================================================

print(f"\nEstrazione profilo LOS (distanza dal ray)...")

# Normalizza la direzione
dir_norm = LOS_DIRECTION / (np.linalg.norm(LOS_DIRECTION) + 1e-10)

# Calcola t e distanza da ray per ogni punto
t_params = []
distances = []

for i, (r, z) in enumerate(zip(R, Z)):
    point = np.array([r, z])
    p_rel = point - LOS_ORIGIN
    
    # t = proiezione del punto sulla linea
    t = np.dot(p_rel, dir_norm)
    
    # distanza perpendiculare
    p_perp = p_rel - t * dir_norm
    dist = np.linalg.norm(p_perp)
    
    t_params.append(t)
    distances.append(dist)

t_params = np.array(t_params)
distances = np.array(distances)

# Seleziona punti vicini al ray (entro 0.1 m)
threshold = 0.15
close_idx = distances < threshold

print(f"Punti vicini al ray (distanza < {threshold} m): {np.sum(close_idx)}")

if np.sum(close_idx) > 0:
    t_close = t_params[close_idx]
    signal_close = signal[close_idx]
    
    # Ordina per t
    sort_idx = np.argsort(t_close)
    t_sorted = t_close[sort_idx]
    signal_sorted = signal_close[sort_idx]
    
    print(f"\nProfilo LOS (primi 50 punti):")
    print(f"{'t':>10} {SIGNAL_NAME:>10} {'distanza':>10}")
    print("-" * 35)
    for i in range(min(50, len(t_sorted))):
        dist_val = distances[close_idx][sort_idx[i]]
        print(f"{t_sorted[i]:10.4f} {signal_sorted[i]:10.4f} {dist_val:10.6f}")
    
    # ============================================================
    # VISUALIZZAZIONE
    # ============================================================
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot sinistra: mesh con LOS
    ax = axes[0]
    scatter = ax.scatter(R, Z, c=signal, cmap='viridis', s=10, alpha=0.6)
    
    # Disegna LOS
    t_range = np.array([t_sorted.min() - 0.2, t_sorted.max() + 0.2])
    los_pts = LOS_ORIGIN[:, np.newaxis] + dir_norm[:, np.newaxis] * t_range
    ax.plot(los_pts[0], los_pts[1], 'r-', linewidth=2, label='LOS')
    
    # Evidenzia punti sulla LOS
    scatter2 = ax.scatter(R[close_idx], Z[close_idx], c='red', s=30, 
                          alpha=0.8, edgecolors='darkred', linewidths=0.5,
                          label='Punti sulla LOS')
    
    ax.set_xlabel('R [m]', fontsize=12)
    ax.set_ylabel('Z [m]', fontsize=12)
    ax.set_title(f'Mesh 2D con LOS')
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label=SIGNAL_NAME)
    
    # Plot destra: profilo 1D
    ax2 = axes[1]
    ax2.plot(t_sorted, signal_sorted, 'o-', linewidth=1.5, markersize=4, color='darkblue')
    ax2.fill_between(t_sorted, signal_sorted, alpha=0.3, color='lightblue')
    ax2.set_xlabel('t (parametro lungo LOS)', fontsize=12)
    ax2.set_ylabel(f'{SIGNAL_NAME}', fontsize=12)
    ax2.set_title(f'Profilo 1D lungo LOS')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = f"los_profile_{SIGNAL_NAME}_final.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot salvato: {output_file}")
    plt.show()
    
    # Salva dati profilo
    profile_file = f"los_profile_{SIGNAL_NAME}_final.txt"
    with open(profile_file, 'w') as f:
        f.write(f"# Profilo LOS per {SIGNAL_NAME}\n")
        f.write(f"# Punto 1: ({LOS_X1}, {LOS_Y1})\n")
        f.write(f"# Punto 2: ({LOS_X2}, {LOS_Y2})\n")
        f.write(f"{'t':>12} {SIGNAL_NAME:>12} {'distanza_ray':>12}\n")
        f.write("-" * 40 + "\n")
        for i in range(len(t_sorted)):
            dist_val = distances[close_idx][sort_idx[i]]
            f.write(f"{t_sorted[i]:12.6f} {signal_sorted[i]:12.6f} {dist_val:12.6f}\n")
    
    print(f"Profilo salvato: {profile_file}")

else:
    print("ERROR: Nessun punto trovato vicino al ray!")

print("\n" + "=" * 80)
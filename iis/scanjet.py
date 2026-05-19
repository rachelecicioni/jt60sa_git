import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# === 1. Dati direttamente nel codice ===

alpha_c_vals = [1.9, 2.0, 2.1, 2.2, 2.3]
chi_elm_vals = [0.5, 0.7, 0.8, 0.9, 1.0]

# Matrice di DeltaEp (righe = alpha_c, colonne = chi_elm)
# Usa 'None' per i valori mancanti
delta_Ep_matrix = [
    [None, 52667, 63176, 74901, 88449],    # alpha_c = 1.9
    [33643, 58650, 70817, 80921, 96400],   # alpha_c = 2.0
    [None, 65862, 77854, 87533, 105583],   # alpha_c = 2.1
    [40376, 70098, None, None, None],      # alpha_c = 2.2
    [None, 71214, 84528, 95239, 114690],   # alpha_c = 2.3
]

# === 2. Converti in lista di punti validi (x, y, z) ===
chi_list = []
alpha_list = []
delta_Ep_list = []

for i, alpha in enumerate(alpha_c_vals):
    for j, chi in enumerate(chi_elm_vals):
        val = delta_Ep_matrix[i][j]
        if val is not None:
            alpha_list.append(alpha)
            chi_list.append(chi)
            delta_Ep_list.append(val)

# === 3. Griglia regolare per interpolazione ===
xi = np.linspace(min(chi_elm_vals), max(chi_elm_vals), 100)
yi = np.linspace(min(alpha_c_vals), max(alpha_c_vals), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpola su griglia
zi = griddata((chi_list, alpha_list), delta_Ep_list, (xi, yi), method='linear')

# === 4. Plot ===
plt.figure(figsize=(8, 6))
cp = plt.contourf(xi, yi, zi, levels=20, cmap='viridis')
contour_lines = plt.contour(xi, yi, zi, levels=[74609], colors='red', linewidths=2)
plt.clabel(contour_lines, inline=True, fmt={74609: 'Target'}, fontsize=10)

# Etichette e layout
plt.colorbar(cp, label=r'$\Delta E_p$ (J)')
plt.xlabel(r'$\chi_{\mathrm{ELM}}$ (m$^2$/s)', fontsize=14)
plt.ylabel(r'$\alpha_c$ (a.u.)', fontsize=14)
plt.title(r'Parametric scan for $\Delta E_{p,target}=74609$ J',fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# === 1. Parametri alpha_c e chi_elm ===
alpha_c_vals = [1.9, 2.0, 2.1, 2.2, 2.3]
chi_elm_vals = [0.5, 0.7, 0.8, 0.9, 1.0]

# === 2. Matrice delle frequenze f (Hz) ===
freq_matrix = [
    [None, 105.90,  92.97, 85.36,  69.36],   # alpha_c = 1.9
    [123.12, 86.72, 76.15, 69.14,  56.09],   # alpha_c = 2.0
    [None,   67.49, 61.49, 55.92,  42.92],   # alpha_c = 2.1
    [88.54,  50.11, None,  None,   None],    # alpha_c = 2.2
    [None,   48.87, 43.70, 39.69,  30.52],   # alpha_c = 2.3
]

# === 3. Prepara i dati validi ===
chi_list = []
alpha_list = []
f_list = []

for i, alpha in enumerate(alpha_c_vals):
    for j, chi in enumerate(chi_elm_vals):
        val = freq_matrix[i][j]
        if val is not None:
            alpha_list.append(alpha)
            chi_list.append(chi)
            f_list.append(val)

# === 4. Griglia e interpolazione ===
xi = np.linspace(min(chi_elm_vals), max(chi_elm_vals), 100)
yi = np.linspace(min(alpha_c_vals), max(alpha_c_vals), 100)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((chi_list, alpha_list), f_list, (xi, yi), method='linear')

# === 5. Plot ===
plt.figure(figsize=(8, 6))
cp = plt.contourf(xi, yi, zi, levels=20, cmap='viridis')
contour_lines = plt.contour(xi, yi, zi, levels=[38.7], colors='red', linewidths=2)
plt.clabel(contour_lines, inline=True, fmt={38.7: 'Target'}, fontsize=10)

# Etichette e layout
plt.colorbar(cp, label=r'$f$ (Hz)')
plt.xlabel(r'$\chi_{\mathrm{ELM}}$ (m$^2$/s)', fontsize=14)
plt.ylabel(r'$\alpha_c$ (a.u.)', fontsize=14)
plt.title(r'Parametric scan for $f_{target}=38.7$ Hz', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()

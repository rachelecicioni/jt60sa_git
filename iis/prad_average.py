import ppf
import numpy as np 
import time  
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt  
from getdat import getdat 
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['black'])

pulse = 70000
user='jnv7243'
ppf.ppfgo()
'''
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
sequenza = input("Inserisci la sequenza: ")


prad_core = ppf.ppfdata(pulse, 'JST', 'PRAD', seq=sequenza, uid=user)
ax1.plot(prad_core[2], prad_core[0], label=sequenza)
ax1.set_ylabel(r'$P_{rad,core}$')
ax1.legend()

zeff = ppf.ppfdata(pulse, 'JST', 'ZEFF', seq=sequenza, uid=user)
ax2.plot(zeff[2], zeff[0], label=sequenza)
ax2.set_ylabel(r'$Z_{eff}$')
ax2.set_xlabel("time (s)")

plt.show()

t_min = float(input("t_min: "))
t_max = float(input("t_max: "))
i_min = np.abs(prad_core[2]-t_min).argmin() #nearest index to t_min
i_max = np.abs(prad_core[2]-t_max).argmin() #nearest index to t_max
prad_core_mean = np.mean(prad_core[0][i_min:i_max])
zeff_mean = np.mean(zeff[0][i_min:i_max])
print(prad_core_mean)
print(zeff_mean)
'''
#Array completati a mano facendo girare questo codice per le simulazioni di interesse
num_seq=[580, 564, 575, 581, 579, 578, 576, 577]
zeff_mean=[1.0611914, 1.0608886, 1.1965913, 1.3612049, 1.1964725, 1.3692366, 1.1957258, 1.3603616]
prad_mean=[-816176.75, -816533.8, -1146021.4, -1573088.1, -1146066.1, -1514653.1, -1149043.5, -1575834.1 ]
puffD_mean=[5e18, 5e19, 5e18, 5e18, 1e19, 1e19, 5e19, 5e19]

# Conversione in numpy array
x = np.array(zeff_mean)
y = np.array(puffD_mean)
z = -np.array(prad_mean)*1e-6

# Creazione scatter 2D con colorbar
plt.figure(figsize=(10, 5))
sc = plt.scatter(x, y, c=z, cmap='plasma', s=100, edgecolors='k')

# Aggiunta della colorbar
cb = plt.colorbar(sc)
#cb.set_label('P_rad (W)')

# Etichette
plt.xlabel(r'$Z_{eff}$', fontsize=14)
plt.ylabel(r'Puff D $(s^{-1})$', fontsize=14)
plt.title(r'$P_{rad,core} (MW)$', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()
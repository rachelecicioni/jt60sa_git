import ppf
import numpy as np 
import time  
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt  
from getdat import getdat 
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['black'])

pulse = 96745
user='jnv7243'
dt=40

ppf.ppfgo()
fig, (ax1, ax2)=plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8,6))
fig.suptitle(r'Simulated ELMs dynamics $D\#96745$ ', fontsize=14)
next_value=0

sequenza=975

wmhd= ppf.ppfdata(pulse, 'JST', 'WTH',seq=sequenza, uid=user)[0]
wmhd_t = ppf.ppfdata(pulse, 'JST', 'WTH',seq=sequenza, uid=user)[2]-dt
ax2.plot(wmhd_t, wmhd*1e-6, label=sequenza, color='blue') #Plasma energy in MJ
ax2.set_ylabel(r'$E_p(MJ)$', fontsize=14)
ax2.set_xlabel("time (s)")

elm = ppf.ppfdata(pulse, 'JST', 'XEBA',seq=sequenza, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XEBA',seq=sequenza, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='blue')
ax1.set_ylabel(r'$\chi_{ELM} (m^2/s)$', fontsize=14)
'''
peaks = np.where(np.diff(np.sign(np.diff(elm))) == -2)[0] + 1 #trovo gli indici corrispondenti ai massimi di elm
tempi_counts=[]
counts=[]
count=0
for i in peaks:
    tempi_counts.append(elm_t[i])
    count=count+1
    counts.append(count)

if len(tempi_counts)>1: #if there are more than one elm
    coefficients = np.polyfit(tempi_counts, counts, 1)
    m_r = coefficients[0]  
    q_r = coefficients[1]

    ax2.scatter(tempi_counts, counts, label='count ELMs')
    ax2.set_ylabel('ELMs (a.u.)')

    x_interpolation = np.linspace(min(tempi_counts), max(tempi_counts), 100)
    y_interpolation = m_r * x_interpolation + q_r
    ax2.plot(x_interpolation, y_interpolation, 'r-', label='Retta di interpolazione')

    peaks_wmhd = np.where(np.diff(np.sign(np.diff(wmhd))) == -2)[0] + 1 #indici di wmhd corrispondenti ai massimi locali

    drops_abs=[]
    for i in range(len(peaks_wmhd)-1):
        minimo=np.min(wmhd[peaks_wmhd[i]:peaks_wmhd[i+1]])
        drops_abs.append(wmhd[peaks_wmhd[i]]-minimo)

    tempi_drops=[]
    for i in range(len(peaks_wmhd)-1):
        tempi_drops.append(wmhd_t[peaks_wmhd[i]])     
    ax4.scatter(tempi_drops, drops_abs, label=r'$\Delta E_p (J)')
    ax4.set_ylabel(r'$\Delta E_p abs (J)$')
    ax4.set_xlabel('time (s)')  
    dropabs_medio=np.mean(drops_abs)   
    print(str(sequenza) + "\t f: " + str(m_r) + "\t" + str(q_r) + "\t" + "\t DE (J): " + str(dropabs_medio))
    
elif len(tempi_counts)==1: #if there is one elm
    print(str(sequenza) + "\t" + "one elm" + "\t" + "one elms" + "\t" + "one elms")
else: #if no elm
    print(str(sequenza) + "\t" + "no elms" + "\t" + "no elms" + "\t" + "no elms")
'''
plt.show()


import ppf
import numpy as np 
import time  
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt  
from getdat import getdat 
import h5py
plt.rcParams['lines.linewidth'] = 1.0
dt=40
tmin=47.95-dt
tmax=50.1-dt
tstart=48.0-dt
tstop=49.986-dt
t2D=48.36-dt
t2T=49.27-dt
t2DT=48.52-dt


pulse = 96745
user_color = 'blue'
isotopo= 'D'
user='JETPPF'
title = f"{isotopo} {pulse}"

ppf.ppfgo()
fig, (ax1, ax2)=plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8,6))
fig.suptitle(r'Experimental ELMs dynamics $D\#96745$ ', fontsize=14)
times_imp = np.loadtxt('times_{}.txt'.format(pulse))-dt #qui prendo i tempi salvati con codice Garzotti
times=times_imp[(times_imp>tmin) & (times_imp<tmax)]

#elm
elm = ppf.ppfdata(pulse, 'EDG7', 'BE2H', uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'EDG7', 'BE2H', uid=user)[2]-dt

ax1.plot(elm_t,elm, color=user_color)
ax1.set_xlim(tmin,t2D+0.05)
ax1.set_ylabel(r'$Be_{II} (ph/s/cm^2/sr)$', fontsize=14)
ax1.set_ylim(0,0.6*1e13)
for i in times:
    ax1.axvline(i, color='dimgrey', linewidth=0.8)
ax1.axvline(tstart, color='green', linestyle='dashed')
ax1.axvline(t2D, color='red', linestyle='dashed')

#wmhd
wmhd = getdat('xg/rtss/wmhd', pulse)[0]
wmhd_t = getdat('xg/rtss/wmhd', pulse)[1]-dt

ax2.plot(wmhd_t,wmhd*1e-6,color=user_color)
ax2.set_xlim(tmin,t2D+0.05)
ax2.set_ylabel(r'$E_p$ (MJ)', fontsize=14)
ax2.set_xlabel('time (s)', fontsize=14)
for i in times:
    ax2.axvline(i, color='dimgrey', linewidth=0.8)
ax2.axvline(tstart, color='green', linestyle='dashed')
ax2.axvline(t2D, color='red', linestyle='dashed')
ax2.set_ylim(3,7)


plt.show()

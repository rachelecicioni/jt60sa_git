import matplotlib.pyplot as plt
import eproc as ep
from jetto_tools.binary import *
from scipy import signal 
import ppf
plt.rcParams['lines.linewidth'] = 1.0
dt=40
tmin=47.9-dt
tmax=50.0-dt
tstart=48.0-dt
tstop=49.986-dt
n=16 #font size
n2=14

pulse = 96745
user='jnv7243'
ppf.ppfgo()

trand='/home/jnv7243/cmg/catalog/edge2d/jet/96745/sep2423/seq#5/tran'
trandts='/home/jnv7243/cmg/catalog/edge2d/jet/96745/sep2223/seq#4/tran'
trant='/home/jnv7243/cmg/catalog/edge2d/jet/96745/sep1523/seq#2/tran'

fig, (ax1, ax2, ax3)=plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(9,6))

elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=886, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=886, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='blue')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=896, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=896, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='blue')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=947, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=947, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='blue', label = 'D')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=847, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=847, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='orangered')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=859, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=859, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='orangered')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=925, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=925, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='orangered', label= 'DT(s.)')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=834, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=834, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='darkmagenta')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=851, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=851, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='darkmagenta')
elm = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=912, uid=user)[0]
elm_t = ppf.ppfdata(pulse, 'JST', 'XIBA',seq=912, uid=user)[2]-dt
ax1.plot(elm_t, elm, color='darkmagenta', label='T')
#ax1.set_ylabel(r'$\chi_{ELM}  (m^2/s)$')
ax1.tick_params(axis='both')

ax1.set_title(r'$\chi_{ELM}  (m^2/s)$', loc='center', fontsize=14)
ax1.set_xlim(9.80,9.85)
ax1.set_ylim(0,1.5)
#ax1.legend()

dato2="YLDT_2"

moltd=ep.time(trand, 'GAMIT')
myd=np.array(moltd.yData)

resultd=ep.time(trand, dato2)
td=np.array(resultd.xData)-dt
yd=np.array(resultd.yData)

moltt=ep.time(trant, 'GAMIT')
myt=np.array(moltt.yData)

resultt=ep.time(trant, dato2)
tt=np.array(resultt.xData)-dt
yt=np.array(resultt.yData)

moltdts=ep.time(trandts, 'GAMIT')
mydts=np.array(moltdts.yData)

resultdts=ep.time(trandts, dato2)
tdts=np.array(resultdts.xData)-dt
ydts=np.array(resultdts.yData)

ax2.plot(td,-yd*myd, color='blue')
ax2.plot(tdts,-ydts*mydts, color='orangered')
ax2.plot(tt,-yt*myt, color='darkmagenta')
ax2.set_title("Tungsten sputtering source (a.u.)", loc='center', fontsize=14)
#ax2.set_ylabel("Tungsten \n sputtering source \n (a.u.)")
#ax2.tick_params(axis='both')
#ax2.fill_between(elm_t, 2.6*1e21,  where=(elm_t >= 9.80) & (elm_t <= 9.85), color='yellow', alpha=0.3)
ax2.set_ylim(0,2.5*1e21)


dato3="FZBND_2"

resultd=ep.time(trand, dato3)
td=np.array(resultd.xData)-dt
yd=np.array(resultd.yData)

resultt=ep.time(trant, dato3)
tt=np.array(resultt.xData)-dt
yt=np.array(resultt.yData)

resultdts=ep.time(trandts, dato3)
tdts=np.array(resultdts.xData)-dt
ydts=np.array(resultdts.yData)

ax3.plot(td,yd, color='blue')   
ax3.plot(tdts,ydts, color='orangered')
ax3.plot(tt,yt, color='darkmagenta')
ax3.tick_params(axis='both')
ax3.set_title("Tungsten core boundary flux ($s^{-1}$)", loc='center', fontsize=14)
ax3.set_ylim(-0.5*1e19, 1.5*1e19)
#ax3.set_ylabel("Tungsten core \n boundary flux \n ($s^{-1}$)")


plt.show()

plt.plot(td,yd, color='blue')   
plt.plot(tdts,ydts, color='orangered')
plt.plot(tt,yt, color='darkmagenta')
plt.xlabel("time (s)", fontsize=14)
#plt.xlabel("Tungsten core \n boundary flux \n ($s^{-1}$)", fontsize=n)
plt.tick_params(axis='both')
#plt.set_title('(c)', loc='left', pad=-14)
plt.ylim(-0.5*1e19, 1.5*1e19)
plt.xlim(9.8168,9.83052)


plt.show()
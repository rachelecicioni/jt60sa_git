import matplotlib.pyplot as plt
import eproc as ep
from jetto_tools.binary import *
from scipy import signal

plt.figure(figsize=(6, 6))

prefix=''


dirs=['/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/nov1824/seq#1/']
#'run_JT-60SA_S2_fullyW_5e19_50x_WAr_30',
#'run_JT-60SA_S2_fullyW_5e19_100x_WAr_30']

  # Array di stringhe per la legenda

#dirs=['run_JT-60SA_S2_fullyW_5e18_50x_WAr_30']

cases=[prefix + directory for directory in dirs]

ep.names(cases[0]+'/tran')

#extracase='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/nov1824/seq#1/'
#extracase='/common/cmg/gcor/edge2d/runs/run_ITER_ramp_up_Q10_W_D_ISS_100_RJ8_1_AYIZ2_05_Wconrelax_nG_30a/'

#cases.append(extracase)

labels = ['IC Ar x10']

#EDGE2D signals
#Time
#sig='YLDT_2'
#sig='ZCONTE_2'
#sig='ZRADSN_2'
#sig='ZLEVSN_2'
#sig='ZRAD_2
#sig='HCONTJ'
#sig='DENSIT'
#sig='TESIT'
#sig='ZTOTAL_2'
#sig='TEBNDJ'
#sig='HBND'
sig='ZRAD_1'

#Rows
#sig='TEVE'
#sig='DENEL'
#sig='DPERP'
#sig='CHII'
#sig='SQEZR_2'
#sig='DENZ02'
#sig='DENZT_2'
#sig='VPINCH'
#sig='ZEFF'

#row='OT'
row='OMP'

#JETTO signals
#sig='NGFR'
#sig='PECE'
#sig='S0'

code='E'
Etype='time'
#Etype='row'
start=12




#Plotting section

for case in cases:
#Plot EDGE2D
  if code == 'E':
# Timetraces
    if Etype == 'time':
      result=ep.time(case+'/tran',sig)
      t=np.array(result.xData)
      y=np.array(result.yData)
# Profiles along rows
    if Etype == 'row':
      result=ep.row(case+'/tran',sig,row)
      t=np.array(result.xData)
      y=np.array(result.yData)
#Plot JETTO
  if code == 'J':
    result=read_binary_file(case+'/jetto.jst')
    t=result['TVEC1'].transpose()[0]
    y=result[sig][0]

  plt.plot(t,y)


#plt.yscale('log')

if code == 'E':
  plt.xlabel(result.xDesc.strip().decode('utf-8')+' ['+result.xUnits.strip().decode('utf-8')+']')
  plt.ylabel(result.yDesc.strip().decode('utf-8')+' ['+result.yUnits.strip().decode('utf-8')+']')

if code == 'J':
  plt.xlabel('t [s]')
  plt.ylabel(sig)


plt.legend(labels, loc='upper left', framealpha=0.25)
plt.show()

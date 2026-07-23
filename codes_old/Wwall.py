import matplotlib.pyplot as plt
import eproc as ep
from jetto_tools.binary import *
from scipy import signal 
import ppf
plt.rcParams['lines.linewidth'] = 1.0

pulse = 70000
user='jnv7243'
ppf.ppfgo()

num_seq=int(input("Inserisci il numero di seq. da plottare: "))
sequenze=[]
trans=[]
nomi=[]
for i in range(num_seq):
    elemento=input("Inserisci l'elemento {}: ".format(i+1))
    sequenze.append(elemento)
    tran=input("Inserisci data/seq#: ")
    trans.append(tran)
    nome=input("Nome legenda: ")
    nomi.append(nome)
print(sequenze)
print(trans)
"""
#plot JETTO
#densità elettronica core
for i in range(num_seq):
    y= ppf.ppfdata(pulse, 'JST', 'NEAV',seq=sequenze[i], uid=user)[0]
    x= ppf.ppfdata(pulse, 'JST', 'NEAV',seq=sequenze[i], uid=user)[2]
    plt.plot(x, y, label=nomi[i]) 
    plt.ylabel('electronic density average in the core (m-3)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_NEAV.pdf', format='pdf')
#profilo densità elettronica core
for i in range(num_seq):
    y = ppf.ppfdata(pulse, 'JSP', 'NE',seq=sequenze[i], uid=user)[0]
    x = ppf.ppfdata(pulse, 'JSP', 'XRHO',seq=sequenze[i], uid=user)[0]
    plt.plot(x, y, label=nomi[i]) 
    plt.ylabel('electronic density in the core (m-3)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JSP_NE.pdf', format='pdf')







for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'ZEFF',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'ZEFF',seq=sequenze[i], uid=user)[2]
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/ZEFF(a.u.)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_ZEFF.pdf', format='pdf')


for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'NEBA',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'NEBA',seq=sequenze[i], uid=user)[2]
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/NEBA (m-3)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_NEBA.pdf', format='pdf')

for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'NEBO',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'NEBO',seq=sequenze[i], uid=user)[2]
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/NEBO (m-3)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_NEBO.pdf', format='pdf')

for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'PRAD',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'PRAD',seq=sequenze[i], uid=user)[2]
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/PRAD (W)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_PRAD.pdf', format='pdf')

for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'WTH',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'WTH',seq=sequenze[i], uid=user)[2
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/WTH(J)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_WTH.pdf', format='pdf')

for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'TEAV',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'TEAV',seq=sequenze[i], uid=user)[2]
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/TEAV(eV)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_TEAV.pdf', format='pdf')

for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'TEBA',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'TEBA',seq=sequenze[i], uid=user)[2]
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/TEBA(eV)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_TEBA.pdf', format='pdf')

for i in range(num_seq):
    zeff= ppf.ppfdata(pulse, 'JST', 'IMC3',seq=sequenze[i], uid=user)[0]
    zeff_t = ppf.ppfdata(pulse, 'JST', 'IMC3',seq=sequenze[i], uid=user)[2]
    plt.plot(zeff_t, zeff, label=nomi[i]) 
    plt.ylabel('JST/IMC3(%)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/JST_IMC3.pdf', format='pdf')

#plot di EDGE2D
for i in range(num_seq):
    name='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/{}/tran'.format(trans[i])
    values=ep.time(name, 'DENSOT')
    t=np.array(values.xData)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('ion density OT (m-3)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/EDGE2D_DENSOT.pdf', format='pdf')

for i in range(num_seq):
    name='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/{}/tran'.format(trans[i])
    values=ep.time(name, 'DEESOT')
    t=np.array(values.xData)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('electron density OT (m-3)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/EDGE2D_DEESOT.pdf', format='pdf')

for i in range(num_seq):
    name='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/{}/tran'.format(trans[i])
    values=ep.time(name, 'TESOT')
    t=np.array(values.xData)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('electron temperature OT (m-3)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/EDGE2D_TESOT.pdf', format='pdf')
"""
for i in range(num_seq):
    name='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/{}/tran'.format(trans[i])
    values=ep.time(name, 'ZCONTE_1')
    t=np.array(values.xData)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('Ar content SOL (?) (a.u.)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/EDGE2D_ZCONTE_3.pdf', format='pdf')

for i in range(num_seq):
    name='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/{}/tran'.format(trans[i])
    values=ep.time(name, 'ZLEVSN_1')
    t=np.array(values.xData)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('Ar content core (a.u.)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/EDGE2D_ZLEVSN_3.pdf', format='pdf')

for i in range(num_seq):
    name='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/{}/tran'.format(trans[i])
    values=ep.time(name, 'ZRAD_1')
    t=np.array(values.xData)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('Ar radiation (a.u.)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/EDGE2D_ZRAD_3.pdf', format='pdf')


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
    tran=input("Inserisci nome cartella: ")
    trans.append(tran)
    nome=input("Nome legenda: ")
    nomi.append(nome)

for i in range(num_seq):
    name='common/cmg/jnv7243/edge2d/runs/{}/tran/pltres'.format(trans[i])
    print(f"Percorso generato per il file: {name}")
    values=ep.time(name, 'ZRAD_1')
    t=np.array(values.xData)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('Ar radiation (a.u.)')
    plt.legend()
plt.show()
plt.savefig('/home/jnv7243/JT-60SA/plot/EDGE2D_ZRAD_3.pdf', format='pdf')
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

def indice_elemento_piu_vicino(array, valore):
    indice = (np.abs(array - valore)).argmin()  # Trova l'indice con la distanza minima
    return indice

for i in range(num_seq):
    name='/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/{}/tran'.format(trans[i])
    values=ep.time(name, 'EMID')
    t=np.array(values.xData)
    print(t)
    y=np.array(values.yData)
    plt.plot(t,y, label=nomi[i])
    plt.ylabel('El. mid plane separatrix density (m-3)')
    plt.legend()
    plt.show()
    t_min=float(input("t_min for the average: "))
    t_max=float(input("t_max for the average: "))
    indice1=indice_elemento_piu_vicino(t,t_min)
    indice2=indice_elemento_piu_vicino(t,t_max)
    mean=np.mean(y[indice1:indice2])
    print(str(num_seq)+"$<n_e>_sep$: {}: ".format(mean))
    plt.plot(t,y, label=nomi[i])
    plt.axvline(x=t_min, linestyle='--')
    plt.axvline(x=t_max, linestyle='--')
    plt.ylabel('El. mid plane separatrix density (m-3)')
    plt.xlabel('time (s) ')
    plt.legend()

    plt.show

plt.show()

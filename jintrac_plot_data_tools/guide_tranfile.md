# Guide for tran file structure and reading

In order to read the tran file it is necessary to import the `eproc` library in python.  

In order to know which signals are present in the tran file use `available_data_run.py <run_path>`

```
import eproc as ep

tran_file = <tran file path>
signal_file = ep.data(tran_file,signal)

#Example on JDC: 
tran_file="/common/cmg/jnv7243/edge2d/runs/run_sa_fk_c_p20_a3_ArW_fc_farsol/tran"
signal = "QPARTOT"

print(dir(signal_file)) #Different field of signal_file
    # nPts, name, time, data, units
print(signal_file) #Informations about the signal
print(type(signal_file)) #Data type of signal_data, usually a struct

```


## Grid of EDGE2D
The tran file`s data necessary in order to reconstruct the mesh of EDGE2D are:  
- **RMESH**: R coordinates of the centroids of the cells
- **ZMESH**: Z coordinates of the centroids of the cells
- **RVERTP**: R coordinates of the centroids *and* vertices of the cellsR coordinates of the centroids *and* vertices of the cells
- **ZVERTP**: Z coordinates of the centroids *and* vertices of the cells

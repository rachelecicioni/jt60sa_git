import eproc as ep
import numpy as np
import os
import matplotlib.pyplot as plt

print("=" * 80)
print("AVAILABLE EPROC FUNCTIONS")
print("=" * 80)

for name in dir(ep):
    if not name.startswith("_"): ## In Python, names starting with "_" are conventionally treated as internal/private
        print(name)

# ============================================================
# RUN DIRECTORY
# ============================================================
RUN_DIR = "/home/jnv7243/cmg/catalog/edge2d/jt60sa/70000/dec0525/seq#1"
TRAN_FILE = os.path.join(RUN_DIR, "tran")

def inspect_datastruct(data, name=""):
    print("\n" + "="*80)
    print(f"DATASTRUCT INSPECTION: {name}")
    print("="*80)

    print("Type:", type(data))
    print("-"*80)

    for attr in dir(data):
        if attr.startswith("_"):
            continue

        try:
            val = getattr(data, attr)

            # try to convert to numpy array for shape inspection
            arr = np.array(val)

            print(f"{attr:<20} shape={arr.shape} dtype={arr.dtype}")

        except Exception:
            # fallback for scalars or weird objects
            try:
                val = getattr(data, attr)
                print(f"{attr:<20} type={type(val)}")
            except:
                pass

def inspect_tran():
    
    print("=" * 80)
    print("READING TRAN FILE \n")
    print(TRAN_FILE)
    print("=" * 80)

    if not os.path.exists(TRAN_FILE):
        print("ERROR: tran file not found")
        return
    
    #List of signals in TRAN_FILE
    signals = ep.names(TRAN_FILE) 
    geom = ep.geom(TRAN_FILE)
    print(geom)
    #geom.keys()
    
    '''
    #Type and dimension of selected signal in TRAN_FILE
    #var_sel = input("Name of signal: ")
    data = ep.data(TRAN_FILE, "TEVE")
    time = ep.time(TRAN_FILE, "TEVE")

    print("signal:", data.name)
    print("desc:", data.desc)
    print("units:", data.units)

    print("nPts =", data.nPts)
    print(data.nPts)
    print("time =", data.time)
    print("timestep =", data.timestep)

    print("len(data) =", len(data.data))

    arr = np.array(data.data)

    print("shape:", arr.shape)
    print("min/max:", arr.min(), arr.max())
    '''

if __name__ == "__main__":
    inspect_tran()


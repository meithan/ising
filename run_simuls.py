# Run a series of ising models
from math import log, sqrt
import os
import numpy as np

Tc = 2/log(1+sqrt(2))

temperatures = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, "Tc", 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5]

base_dir = "simuls"

for T in temperatures:
  
  if T == "Tc":
    out_dir = os.path.join(base_dir, "Tc")
    T = Tc
  else:
    out_dir = os.path.join(base_dir, f"T{T:.1f}")
  os.makedirs(out_dir, exist_ok=True)
  
  cmd = f"./ising {T}"
  print(cmd)
  os.system(cmd)

  cmd = f"mv T{T:.3f}_*.dat {out_dir}"
  print(cmd)
  os.system(cmd)

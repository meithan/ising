import os
import sys
import matplotlib.pyplot as plt
import numpy as np

fname = sys.argv[1]

magn = []; energy = []
with open(fname) as f:
  
  # Detrermine NGRID
  while True:
    line = f.readline()
    if line == "":
      print("Coudln't determine grid size")
      sys.exit()
    if "Grid N =" in line:
      NGRID = int(line.strip().split()[-1])
      break
  
  print(f"{NGRID} x {NGRID} grid")

  # Read and plot grids
  while True:
    line = f.readline()
    if line == "": break
    if "GEN" in line:
      
      gen = int(line.strip().split()[-1])
      
      rows = []
      for i in range(NGRID):
        rows.append([int(x) for x in f.readline().strip()])
      grid = np.array(rows, dtype=int)
      
      plt.figure(figsize=(8,8))
      grid[grid == 0] = -1
      
      plt.imshow(grid, origin="lower", cmap="Greys_r")

      plt.axis("off")
      
      plt.title(f"Generation {gen}")
      plt.tight_layout()

      if "--save" in sys.argv:
        out_fname = os.path.splitext(fname)[0] + f"_gen{gen}" + ".png"
        plt.savefig(out_fname)
        print("Wrote", out_fname)
      else:
        plt.show()

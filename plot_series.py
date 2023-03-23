import os
import sys
import matplotlib.pyplot as plt

fname = sys.argv[1]

magn = []; energy = []
with open(fname) as f:
  for line in f:
    if line.startswith("#"): continue
    magn.append(float(line.split()[1]))
    energy.append(float(line.split()[2]))

# plt.figure(figsize=(10,4))
plt.figure(figsize=(10,7))

plt.subplot(2,1,1)
plt.plot(magn, color="C0")
plt.xlabel("Generation")
plt.title("Magnetization")
plt.grid(ls=":")

plt.subplot(2,1,2)
plt.plot(energy, color="C1")
plt.xlabel("Generation")
plt.title("Energy")
plt.grid(ls=":")

plt.tight_layout()

if "--save" in sys.argv:
  out_fname = os.path.splitext(fname)[0] + ".png"
  plt.savefig(out_fname)
  print("Wrote", out_fname)
else:
  plt.show()
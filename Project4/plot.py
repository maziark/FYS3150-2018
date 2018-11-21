import numpy as np
import math
import sys
import matplotlib.pyplot as plt


f1 = "configsRanT1.txt"
f2 = "configsHomoT1.txt"
x = 0
y1 = 0
y2 = 0
x_label = "Monto Carlo spins"
y_label = "Accepted Configurations"

#if (len(sys.argv) > 1):
#	y1 = int(sys.argv[2])
	#y2 = int(sys.argv[3])


file = open(f1)
lines1 = [line.rstrip('\n') for line in file]
file.close()

file = open(f2)
lines2 = [line.rstrip('\n') for line in file]
file.close()


rows1 = [l.split() for l in lines1]
rows1 = [[float(x) for x in y] for y in rows1]

rows2 = [l.split() for l in lines2]
rows2 = [[float(x) for x in y] for y in rows2]
#print (rows)

#xs = [r[x] for r in rows]
ys1 = [r[y1] for r in rows1]
ys2 = [r[y2] for r in rows2]
xs = list(range(len(ys1)))
plt.xlabel(x_label)
plt.ylabel(y_label)

plt.plot (xs, ys1, 'g', label="Homogeneos")
plt.plot (xs, ys2, 'b', label="Random")
plt.legend()
plt.show()

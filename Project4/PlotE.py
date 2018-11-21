import numpy as np
import math
import sys
import matplotlib.pyplot as plt


f = "result"
x = 0
y1 = 1
y2 = 2
x_label = "Monto Carlo Spins Spins"
y_label = "Numerical value of E, |M|"

if (len(sys.argv) > 1):
	f = sys.argv[1]

if (len(sys.argv) > 3):
	y1 = int(sys.argv[2])
	y2 = int(sys.argv[3])


file = open(f)
lines = [line.rstrip('\n') for line in file]
file.close()


rows = [l.split() for l in lines]
rows = [[float(x) for x in y] for y in rows]
#print (rows)

#xs = [r[x] for r in rows]
ys1 = [r[y1] for r in rows]
ys2 = [r[y2] for r in rows]
xs = list(range(len(ys1)))
plt.xlabel(x_label)
plt.ylabel(y_label)

plt.plot (xs, ys1, 'g', label="E")
plt.plot (xs, ys2, 'b', label="|M|")
plt.legend()
plt.show()

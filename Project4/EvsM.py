import numpy as np
import math
import sys
import matplotlib.pyplot as plt


f = "result"
x = 0
y1 = 3
y2 = 4

x_label = "Monto Carlo Spins"
y_label = "Energy - Magnetization "

if (len(sys.argv) > 1):
	f = sys.argv[1]



file = open(f)
lines = [line.rstrip('\n') for line in file]
file.close()


rows = [l.split() for l in lines]
rows = [[float(x) for x in y] for y in rows]
#print (rows)

xs = [r[x] for r in rows]
ys1 = [r[y1] for r in rows]
ys2 = [r[y2] for r in rows]
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.plot (xs, ys1, label='Energy')
plt.plot (xs, ys2, label='Magnetization')
plt.show()

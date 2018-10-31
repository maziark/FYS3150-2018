import numpy as np
import math
import sys
import matplotlib.pyplot as plt

n = 3
d = 3
f = "body2.txt"
if (len(sys.argv) > 2):
	n = int(sys.argv[1])
	d = int(sys.argv[2])

if (len(sys.argv) > 3):
	f = sys.argv[3]

fig = plt.figure()
planets = []
# number of planets;

file = open(f)

lines = [line.rstrip('\n') for line in file]
file.close()	
planets = []
for r in range(n) : 
	planets.append([lines[i] for i in range(len(lines)) if i%n==r])

for r in range (d):
	name  = [a.split()[0] for a in planets[r]]
	pos_x = [a.split()[2] for a in planets[r]]
	pos_y = [a.split()[3] for a in planets[r]]
	plt.plot (pos_x, pos_y, label = name[0])

plt.xlabel('x(AU)')
plt.ylabel('y(AU)')
plt.legend()

plt.show()



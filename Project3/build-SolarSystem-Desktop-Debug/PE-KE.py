import numpy as np
import math
import matplotlib.pyplot as plt
p = 1

def plotMe (s, o_0, o_1, o_2, o_3, filename):
	ax1 = plt.subplot(2, 1, 1)
	ax1.set_title("Euler")
	plt.semilogy(s, o_0, 'b', label="Potential Energy")
	plt.semilogy(s, o_1, 'r', label="Kinetic Energy")
	#plt.plot(s, o_0, 'b', label="Potential Energy")
	#plt.plot(s, o_1, 'r', label="Kinetic Energy")

	plt.legend()
	plt.xlabel('Time [year]')
	plt.ylabel('Energy [J]')


	ax2 = plt.subplot(2, 1, 2)
	ax2.set_title("Velocity Verlet")
	#plt.plot(s, o_2, 'g', label="Potential Energy")
	#plt.plot(s, o_3, 'k', label="Kinetic Energy")
	plt.semilogy(s, o_2, 'g', label="Potential Energy")
	plt.semilogy(s, o_3, 'k', label="Kinetic Energy")

	plt.legend()
	plt.xlabel('Time [year]')
	plt.ylabel('Energy [J]')
	#plt.show()
	plt.savefig(filename+ ".jpg")
	plt.close()

for i in range (4):
	p = 10**i

	filename = "Report/twoBodyEnergy_"+str(p)
	file = open(filename+ ".txt", "r")
	obj = file.read() [:-1]
	obj = obj.replace ('\n', ' ')

	obj = obj.split(' ')

	n = [0.01, 0.5, 1.0, 5.0]
	o_0 = []
	o_1 = []
	o_2 = []
	o_3 = []

	i = 0
	for x in obj:
		print (x)
		if (i%4 == 0):
			o_0.append(float(x))
			print (n)
		elif (i%4 == 1):
			o_1.append(float(x))
		elif (i%4 == 2):
			o_2.append(float(x))
		else:
			o_3.append(float(x))
			print (n)
		i+=1

	s = list(range(len(o_0)))
	plotMe (s, o_0, o_1, o_2, o_3, filename)
	file.close()

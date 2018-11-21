from __future__ import division
import plotly.plotly as py
import plotly.tools as tls

import matplotlib.pyplot as plt
import numpy as np

from matplotlib.pyplot import *
from numpy import *

E_file = open("output.txt")

E = []
temp = []

#For T=1:
possible_E = linspace(-720.,-800., 21) 

#For T = 2:
#possible_E = linspace(-540., -840., 76)

#For T = 2.4:
#possible_E = linspace(-340., -640., 76) #Gir DE = 4

size_E = len(possible_E)  
E_count_vec = zeros(size_E) 
prob_vec = zeros(size_E) 


count = 0
n_spins = 20

#E_file.readline()

#Leser fra E_file:
for line in E_file:
    numbers = line.split()
    temp.append(float(numbers[0]))
    E.append(float(numbers[1]))
    count+=1
E_file.close()


#Counts number of times each of the energies occur in order to calculate the probability for a given energy
for i in range(len(E)):
    #E_count_vec[possible_E.index(E[i])] += 1
    for j in range(size_E):    
        if E[i] == possible_E[j]: 
            E_count_vec[j] += 1


    
prob_vec = E_count_vec / sum(E_count_vec) #Vector with all probabilities
prob_list = list(prob_vec)











gaussian_numbers = np.random.randn(1000)
plt.hist(gaussian_numbers)
plt.title("Gaussian Histogram")
plt.xlabel("Value")
plt.ylabel("Frequency")

fig = plt.gcf()
plotly_fig = tls.mpl_to_plotly( fig )
py.iplot(plotly_fig, filename='mpl-basic-histogram')

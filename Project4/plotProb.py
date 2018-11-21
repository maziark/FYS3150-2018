from __future__ import division
from matplotlib.pyplot import *
from numpy import *

E_file = open("output.txt")

E = []
temp = []

#For T=1:
#possible_E = linspace(-720.,-800., 21) 

#For T = 2:
#possible_E = linspace(-540., -840., 76)

#For T = 2.4:
possible_E = linspace(-340., -640., 76)

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



for i in range(len(E)):
    #E_count_vec[possible_E.index(E[i])] += 1
    for j in range(size_E):    
        if E[i] == possible_E[j]: 
            E_count_vec[j] += 1


    
prob_vec = E_count_vec / sum(E_count_vec) 
prob_list = list(prob_vec)

#hist2d(possible_E, prob_vec, range=((-800, -720), (0, 1)))
plot(possible_E, prob_vec,'*') 
title("Energy-values for 20x20 spins and T=2.4")
xlabel("Energy")
ylabel("Probability Distribution")

show()
    

#mplt.show()

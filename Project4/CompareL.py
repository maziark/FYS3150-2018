# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
from numpy import *

E = []
E_ = []

E2 = []
E2_ = []

M = []
M_ = []

M2 = []
M2_ = []

T = []
T_ = []

S = []
S_ = []

C_v = []
C_v_ = []
file1 = open("outputL40.txt", "r")
file2 = open("output60.txt", "r")


for line in file1:
    numbers = line.split()
    T.append(float(numbers[0]))
    E.append(float(numbers[1]))
    E2.append(float(numbers[2]))
    M.append(float(numbers[3]))
    M2.append(float(numbers[4]))
    C_v.append(float(numbers[5]))
    S.append(float(numbers[6]))

for line in file2:
    numbers = line.split()
    T_.append(float(numbers[0]))
    E_.append(float(numbers[1]))
    E2_.append(float(numbers[2]))
    M_.append(float(numbers[3]))
    M2_.append(float(numbers[4]))
    C_v_.append(float(numbers[5]))
    S_.append(float(numbers[6]))


file1.close()
file2.close()

plot(T, E, label = 'L = 40')
plot(T, E_, label = 'L = 60')
xlabel ('Temperature')
ylabel ('Energy')

"""plot(T, M, label = 'L = 40')
plot(T, M_, label = 'L = 60')
xlabel ('Temperature')
ylabel ('|M|')
"""

"""plot(T, C_v, label = 'L = 40')
plot(T, C_v_, label = 'L = 60')
xlabel ('Temperature')
ylabel ('C_v')"""


"""plot(T, S, label = 'L = 40')
plot(T, S_, label = 'L = 60')
xlabel ('Temperature')
ylabel ('Susceptibility')
"""

legend()
show()


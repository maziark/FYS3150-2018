# -*- coding: utf-8 -*-

from __future__ import division
from matplotlib.pyplot import *
from numpy import *

E = []
E2 = []
M = []
M2 = []
T = []
S = []
C_v = []

file = open("output.txt", "r")

for line in file:
    numbers = line.split()
    T.append(float(numbers[0]))
    E.append(float(numbers[1]))
    E2.append(float(numbers[2]))
    M.append(float(numbers[3]))
    M2.append(float(numbers[4]))
    C_v.append(float(numbers[5]))
    S.append(float(numbers[6]))

file.close()

"""plot(T, E)
xlabel ('Temperature')
ylabel ('Energy')"""

"""plot(T, M)
xlabel ('Temperature')
ylabel ('|M|')"""


"""plot(T, C_v)
xlabel ('Temperature')
ylabel ('C_v')"""


plot(T, S)
xlabel ('Temperature')
ylabel ('Susceptibility')


"""plot(T, E)
xlabel ('Temperature')
ylabel ('Energy')"""


show()


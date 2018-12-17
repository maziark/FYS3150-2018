import numpy as np
import matplotlib.pyplot as plt
from math import *

def create_lambda(A, w, C):
    return lambda t: A * sin(w * t) + C

class SIRS:
    N = 1.0e5
    a = 1.0
    b = 1.0
    c = 1.0
    # For more advanced model
    e = 0.0
    d = 0.0
    d_i = 0.0
    dt = 1e-5

    def __init__ (self, N, a, b, c,
                d = create_lambda(0, 0, 0),
                d_i = create_lambda(0, 0, 0),
                e = create_lambda(0, 0, 0)) :
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.d_i = d_i
        self.N = N

        c1 = 4.0/(self.a(0) * self.N)
        c2 = 1.0/(self.b(0) * self.N)
        c3 = 1.0/(self.c(0) * self.N)

        self.dt = min (c1, min(c2, c3))




    def dS (self, R , I, S):
        result = (self.c(0)* R - self.a(0)* S * I / self.N - self.d(0)* S + self.e(0)* self.N)
        return result

    def dI (self, S, I):
        result = (self.a(0)* S * I / self.N - self.b(0)* I - self.d(0)* I - self.d_i(0) * I)
        return result

    def dR (self, I, R):
        result = (self.b(0)* I - self.c(0)* R - self.d(0)* R)
        return result

    def seasonal_a(self, t):
        #return self.a * sin(self.omega * t) + self.h
        return self.a

    def RK4 (self, S, I , time_frame):
        S0 = S
        I0 = I
        R0 = self.N - S0 - I0

        t_array = np.arange(0, time_frame, step = self.dt)
        S_array = np.zeros(len(t_array))
        I_array = np.zeros(len(t_array))
        R_array = np.zeros(len(t_array))
        counter = 0
        t = 0.0

        #while (t < time_frame):
        for t in t_array:
            S_array[counter] = S0
            I_array[counter] = I0
            R_array[counter] = R0
            counter += 1
            # RK step 1

            S1 = self.dt * self.dS(R0, I0, S0)
            I1 = self.dt * self.dI(S0, I0)
            R1 = self.dt * self.dR(I0, R0)

            # RK step 2

            S2 = self.dt * self.dS(R0 + R1 * 0.5 , I0 + I1 * 0.5, S0 + S1 * 0.5)
            I2 = self.dt * self.dI(S0 + S1 * 0.5 , I0 + I1 * 0.5)
            R2 = self.dt * self.dR(I0 + I1 * 0.5 , R0 + R1 * 0.5)

            # RK step 3
            S3 = self.dt * self.dS(R0 + R2 * 0.5 , I0 + I2 * 0.5, S0 + S2 * 0.5)
            I3 = self.dt * self.dI(S0 + S2 * 0.5 , I0 + I2 * 0.5)
            R3 = self.dt * self.dR(I0 + I2 * 0.5 , R0 + R2 * 0.5)

            # RK step 4
            S4 = self.dt * self.dS(R0 + R3, I0 + I3, S0 + S3)
            I4 = self.dt * self.dI(S0 + S3, I0 + I3)
            R4 = self.dt * self.dR(I0 + I3, R0 + R3)


            S_new = 1.0 * (S1 + 2.0*S2 + 2.0*S3 + S4) / 6.0
            I_new = 1.0 * (I1 + 2.0*I2 + 2.0*I3 + I4) / 6.0
            R_new = 1.0 * (R1 + 2.0*R2 + 2.0*R3 + R4) / 6.0

            S0 += S_new
            I0 += I_new
            R0 += R_new


        return (t_array, S_array, I_array, R_array)


    def MC_ (self, S_init, I_init, time_frame, n_samples = 10):
        S = S_init
        I = I_init
        R = self.N - S - I

        t_array = np.arange(0, time_frame, step = self.dt)
        S_array = np.zeros((len(t_array), n_samples))
        I_array = np.zeros((len(t_array), n_samples))
        R_array = np.zeros((len(t_array), n_samples))

        S_avg = np.zeros (len(t_array))
        I_avg = np.zeros (len(t_array))
        R_avg = np.zeros (len(t_array))

        for i in range(n_samples):
            S = S_init
            I = I_init
            R = self.N - S_init - I_init
            for j in range (len(t_array)):

                # Generate Random Number :
                r = np.random.random(3)

                # P (S -> I) : S-- I++
                if (r[0] < (self.dt * self.a(0)* S * I / self.N) and (S > 0)):
                    S -= 1
                    I += 1

                # P (I -> R) : I-- R++
                if (r[1] < (self.dt * self.b(0)* I ) and (I > 0)) :
                    I -= 1
                    R += 1


                # P (R -> S) : R-- S++
                if (r[2] < (self.dt * self.c(0)* R ) and (R > 0)) :
                    R -= 1
                    S += 1

                if (S + I + R > self.N) : print (S, I, R)

                S_array[j][i] = S
                I_array[j][i] = I
                R_array[j][i] = R

        print (S_array)

        for i in range(len(S_avg)) : S_avg[i] = np.mean (S_array[i])
        for i in range(len(I_avg)) : I_avg[i] = np.mean (I_array[i])
        for i in range(len(R_avg)) : R_avg[i] = np.mean (R_array[i])

        return (t_array, S_avg, S_array, I_avg, I_array, R_avg, R_array)




    def MC_vital (self, S_init, I_init, time_frame, n_samples = 10):
        S = S_init
        I = I_init
        R = self.N - S - I



        t_array = np.arange(0, time_frame, step = self.dt)
        S_array = np.zeros((len(t_array), n_samples))
        I_array = np.zeros((len(t_array), n_samples))
        R_array = np.zeros((len(t_array), n_samples))

        S_avg = np.zeros (len(t_array))
        I_avg = np.zeros (len(t_array))
        R_avg = np.zeros (len(t_array))

        for i in range(n_samples):
            S = S_init
            I = I_init
            R = self.N - S_init - I_init
            t = 0
            for j in range (len(t_array)):


                while (t + self.dt < t_array[j]):
                    t += self.dt
                    # Generate Random Number :
                    r = np.random.random(8)
                    N = S + I + R
                    # Here we have more actions
                    # P(newBorn) [0]
                    if (r[0] < (self.e(t) * N * self.dt)) :
                        S += 1

                    # P(S -> Death) [1]
                    if (S > 0 and r[1] < (self.d(t) * S * self.dt)) :
                        S -= 1

                    # P(S -> I) [2]
                    if (S > 0 and r[2] < (self.a(t) * S * I / N * self.dt)):
                        I += 1
                        S -= 1

                    # P(I -> Death_I) [3]
                    if (I > 0 and r[3] < (self.dt * self.d_i(t) * I)):
                        I -= 1

                    # P(I -> Death) [4]
                    if (I > 0 and r[4] < (self.dt * self.d(t) * I)):
                        I -= 1

                    # P(I -> R) [5]
                    if (I > 0 and r[5] < (self.dt * self.b(t) * I)):
                        I -= 1
                        R += 1

                    # P(R -> Death) [6]
                    if (R > 0 and r[6] < (self.dt * self.d(t) * R)) :
                        R -= 1

                    # P(R -> S) [7]
                    if (R > 0 and r[7] < (self.dt * self.c(t) * R)) :
                        R -= 1
                        S += 1

                    N = S + I + R
                    delta_t  = []
                    if (N > 0) :
                        if (self.a(t) > 0) : delta_t.append(4.0 / (self.a(t) * N))
                        if (self.b(t) > 0) : delta_t.append(1.0 / (self.b(t) * N))
                        if (self.c(t) > 0) : delta_t.append(1.0 / (self.c(t) * N))
                        if (self.e(t) > 0) : delta_t.append(1.0 / (self.e(t) * N))
                        if (self.d(t) > 0) : delta_t.append(1.0 / (self.d(t) * N))
                        if (self.d(t) + self.d_i(t) > 0) : delta_t.append(1.0 / ((self.d(t) + self.d_i(t)) * N))

                    if (len(delta_t) > 0) : self.dt = min(delta_t)

                S_array[j][i] = S
                I_array[j][i] = I
                R_array[j][i] = R

        print (S_array)

        for i in range(len(S_avg)) : S_avg[i] = np.mean (S_array[i])
        for i in range(len(I_avg)) : I_avg[i] = np.mean (I_array[i])
        for i in range(len(R_avg)) : R_avg[i] = np.mean (R_array[i])

        return (t_array, S_avg, S_array, I_avg, I_array, R_avg, R_array)

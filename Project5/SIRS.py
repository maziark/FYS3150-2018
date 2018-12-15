import numpy as np
import matplotlib.pyplot as plt
from math import *
class SIRS:
    N = 1.0e5
    a = 1.0
    b = 1.0
    c = 1.0
    # For more advanced model
    e = 0.0
    d = 0.0
    d_l = 0.0
    dt = 1e-5

    def __init__ (self, N, a, b, c, d = 0, d_l = 0, e = 0) :
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.d_l = d_l
        self.N = N

        c1 = 4.0/(self.a * self.N)
        c2 = 1.0/(self.b * self.N)
        c3 = 1.0/(self.c * self.N)

        self.dt = min (c1, min(c2, c3))




    def dS (self, R , I, S):
        result = (self.c * R - self.a * S * I / self.N - self.d * S + self.e * self.N)
        return result

    def dI (self, S, I):
        result = (self.a * S * I / self.N - self.b * I - self.d * I - self.d_l * I)
        return result

    def dR (self, I, R):
        result = (self.b * I - self.c * R - self.d * R)
        return result


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
                if (r[0] < (self.dt * self.a * S * I / self.N) and (S > 0)):
                    S -= 1
                    I += 1

                # P (I -> R) : I-- R++
                if (r[1] < (self.dt * self.b * I ) and (I > 0)) :
                    I -= 1
                    R += 1


                # P (R -> S) : R-- S++
                if (r[2] < (self.dt * self.c * R ) and (R > 0)) :
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

        return (t_array, S_array, I_avg, R_avg)

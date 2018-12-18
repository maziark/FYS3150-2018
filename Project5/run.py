from SIRS import *
from math import *

def MC_deadly ():
    a = create_lambda(0, 0, 4.0)
    b = create_lambda(0, 0, 1.0)
    c = create_lambda(0,0, 0.5)
    d = create_lambda(0,0, 0.1)

    d_i_A = create_lambda(0,0, 0.0)
    d_i_B = create_lambda(0,0, 1.0)
    d_i_C = create_lambda(0,0, 1.5)
    d_i_D = create_lambda(0,0, 2.0)

    e = create_lambda(0,0, 0.1)

    sirs_A = SIRS(400.0, a, b, c, d, d_i_A, e)
    sirs_B = SIRS(400.0, a, b, c, d, d_i_B, e)
    sirs_C = SIRS(400.0, a, b, c, d, d_i_C, e)
    sirs_D = SIRS(400.0, a, b, c, d, d_i_D, e)


    (t_array_A, S_avg_A, S_array_A, I_avg_A, I_array_A, R_avg_A, R_array_A) = sirs_A.MC_vital(300.0, 100.0, 100.0)
    (t_array_B, S_avg_B, S_array_B, I_avg_B, I_array_B, R_avg_B, R_array_B) = sirs_B.MC_vital(300.0, 100.0, 100.0)
    (t_array_C, S_avg_C, S_array_C, I_avg_C, I_array_C, R_avg_C, R_array_C) = sirs_C.MC_vital(300.0, 100.0, 100.0)
    (t_array_D, S_avg_D, S_array_D, I_avg_D, I_array_D, R_avg_D, R_array_D) = sirs_D.MC_vital(300.0, 100.0, 100.0)


    plt.subplot(221)
    plt.plot (t_array_A, S_avg_A, color = '#F4C430', label="S")
    plt.plot (t_array_A, I_avg_A, color = 'indigo', label="I")
    plt.plot (t_array_A, R_avg_A, color = 'red', label="R")

    plt.plot (t_array_A, S_array_A, '-', color = '#F4C430', alpha= 0.3) # #F4C430 = saffron
    plt.plot (t_array_A, I_array_A, '-', color = 'indigo', alpha= 0.3) # indigo for I
    plt.plot (t_array_A, R_array_A, '-', color = 'red', alpha= 0.3) # red for R

    plt.title("a = " + str(a(0)) + ", b = " + str(b(0)) + ", c = " + str(c(0)) + ", d = " + str(d(0))
        + ", d_i = " + str(d_i_A(0)) + ", e = " + str(e(0)))
    plt.xlabel ('Time')
    plt.ylabel ('Population')
    N = S_avg_A[-1] + I_avg_A[-1] + R_avg_A[-1]
    print ("population A : S/N = ", S_avg_A[-1]/N)
    print ("population A : I/N = ", I_avg_A[-1]/N)
    print ("population A : R/N = ", R_avg_A[-1]/N)

    plt.legend()

    plt.subplot(222)
    plt.plot (t_array_B, S_avg_B, color = '#F4C430', label="S")
    plt.plot (t_array_B, I_avg_B, color = 'indigo', label="I")
    plt.plot (t_array_B, R_avg_B, color = 'red', label="R")

    plt.plot (t_array_B, S_array_B, '-', color = '#F4C430', alpha= 0.3) # #F4C430 = saffron
    plt.plot (t_array_B, I_array_B, '-', color = 'indigo', alpha= 0.3) # indigo for I
    plt.plot (t_array_B, R_array_B, '-', color = 'red', alpha= 0.3) # red for R

    plt.title("a = " + str(a(0)) + ", b = " + str(b(0)) + ", c = " + str(c(0)) + ", d = " + str(d(0))
        + ", d_i = " + str(d_i_B(0)) + ", e = " + str(e(0)))
    plt.xlabel ('Time')
    plt.ylabel ('Population')
    N = S_avg_A[-1] + I_avg_A[-1] + R_avg_A[-1]
    print ("population B : S/N = ", S_avg_A[-1]/N)
    print ("population B : I/N = ", I_avg_A[-1]/N)
    print ("population B : R/N = ", R_avg_A[-1]/N)

    plt.legend()

    plt.subplot(223)
    plt.plot (t_array_C, S_avg_C, color = '#F4C430', label="S")
    plt.plot (t_array_C, I_avg_C, color = 'indigo', label="I")
    plt.plot (t_array_C, R_avg_C, color = 'red', label="R")

    plt.plot (t_array_C, S_array_C, '-', color = '#F4C430', alpha= 0.3) # #F4C430 = saffron
    plt.plot (t_array_C, I_array_C, '-', color = 'indigo', alpha= 0.3) # indigo for I
    plt.plot (t_array_C, R_array_C, '-', color = 'red', alpha= 0.3) # red for R

    plt.title("a = " + str(a(0)) + ", b = " + str(b(0)) + ", c = " + str(c(0)) + ", d = " + str(d(0))
        + ", d_i = " + str(d_i_C(0)) + ", e = " + str(e(0)))
    plt.xlabel ('Time')
    plt.ylabel ('Population')
    N = S_avg_A[-1] + I_avg_A[-1] + R_avg_A[-1]
    print ("population C : S/N = ", S_avg_A[-1]/N)
    print ("population C : I/N = ", I_avg_A[-1]/N)
    print ("population C : R/N = ", R_avg_A[-1]/N)

    plt.legend()

    plt.subplot(224)
    plt.plot (t_array_D, S_avg_D, color = '#F4C430', label="S")
    plt.plot (t_array_D, I_avg_D, color = 'indigo', label="I")
    plt.plot (t_array_D, R_avg_D, color = 'red', label="R")

    plt.plot (t_array_D, S_array_D, '-', color = '#F4C430', alpha= 0.3) # #F4C430 = saffron
    plt.plot (t_array_D, I_array_D, '-', color = 'indigo', alpha= 0.3) # indigo for I
    plt.plot (t_array_D, R_array_D, '-', color = 'red', alpha= 0.3) # red for R

    plt.title("a = " + str(a(0)) + ", b = " + str(b(0)) + ", c = " + str(c(0)) + ", d = " + str(d(0))
        + ", d_i = " + str(d_i_D(0)) + ", e = " + str(e(0)))
    plt.xlabel ('Time')
    plt.ylabel ('Population')
    N = S_avg_A[-1] + I_avg_A[-1] + R_avg_A[-1]
    print ("population A : S/N = ", S_avg_A[-1]/N)
    print ("population A : I/N = ", I_avg_A[-1]/N)
    print ("population A : R/N = ", R_avg_A[-1]/N)

    plt.legend()


    plt.show()

def main ():
    MC_deadly()
    #MC_deadly()

main ()

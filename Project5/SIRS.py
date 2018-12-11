import numpy as np

class SIRS:
    N = 1e5;
    a = 1.0;
    b = 1.0;
    c = 1.0;
    dt = 1e-5;

    def __init__ (self, N, a, b, c) :
        self.a = a;
        self.b = b;
        self.c = c;
        self.N = N;

        c1 = 4.0/(self.a * self.N);
        c2 = 1.0/(self.b * self.N);
        c3 = 1.0/(self.c * self.N);

        self.dt = min (c1, min(c2, c3));

    

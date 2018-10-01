#include <cmath>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <string>
#include <time.h>
#include <armadillo>
#include "time.h"


# define TYPE double
# define M_PI           3.14159265358979323846  /* pi */
# define EPSILON	10e-8;
# define RHO_MIN	0.0;
# define RHO_MAX	1.0;


using namespace std;
using namespace arma;

struct Eig {
	vec	eigvalue;
	mat	eigvector;
	mat	S;
	mat	D;
};


TYPE find_a (int N);
TYPE find_d (int N, int i = 0, int task = 0, TYPE w_r = 0);
TYPE **allocateMemory (int N, int M, bool I = false);
void deallocateMemory (int N, int M, TYPE **A);
TYPE max_off_diagonal (int N, int *x, int *y, TYPE **A);
void setA(TYPE **A, int N);
void jacobi_rotate (int N, int l, int k, TYPE **A, TYPE **R);
int jacobi_method(int N, TYPE **A, TYPE **R);



/* ************************************ */
/* ARMADILLOOOOOOOOOOOOOOOOOOOOOOOOOOOO */
/* ************************************ */


mat generate_mat (int N, TYPE h, TYPE a, TYPE d);
void arma_diag (Eig *eig);
Eig arma_eig (int N, TYPE a, TYPE d);




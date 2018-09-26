#include <iostream>
#include <cmath>

#define TYPE double

using namespace std;

TYPE max_off_diagonal (int N, int *x, int *y, TYPE **A);
void setA(TYPE **A, int N);
void jacobi_rotate (int N, int l, int k, TYPE **A, TYPE **R);
int jacobi_method(int N, TYPE **A, TYPE **R);


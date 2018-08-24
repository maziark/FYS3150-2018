#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <string>
#include <time.h>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

#define TYPE double

void solve_LU				(int N, TYPE *v, TYPE *u);
void solve_matrix_equation_general	(int N, TYPE *v, TYPE *f, TYPE *a, TYPE *b, TYPE *c);
void solve_matrix_equation		(int N, TYPE *v, TYPE *f);
void initiate_abc			(int N, TYPE *a, TYPE *b, TYPE *c, TYPE *f);

TYPE get_u(TYPE x){
	return (1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x));
}

void initiate_abc			(int N, TYPE *a, TYPE *b, TYPE *c, TYPE *f){
	TYPE h = 1.0/(N+1.0);
	cout << "H : " << h << endl;
	for (int i = 0; i < N; i++){
		if (i > 0 && i < N-1) {
			a[i] = -1.0;
			c[i] = -1.0;
		}
		b[i] = 2.0;	
		f[i] = pow(h, 2) * 100 * exp(-10 * (i * h));
	}
}


/*void solve_matrix_equation		(int N, TYPE *v, TYPE *f){
	TYPE	btemp;
	TYPE 	*temp = new TYPE [N+2];
	
	btemp = 2;
	for (int i = 2; i <= N; i++){
		temp[i] = c[i - 1]/btemp;
		btemp = (i + 1.0)/(1.0 * i);
		
	}
}*/

void solve_matrix_equation_general	(int N, TYPE *v, TYPE *f, TYPE *a, TYPE *b, TYPE *c){
	TYPE 	*temp = new TYPE [N+2];
	TYPE 	btemp;
	
	btemp = b[0];
	f[0] = f[0]/btemp;
	for (int i = 1; i < N; i++){
		temp[i] = c[i - 1] / btemp;
		f[i] = f[i] - (a[i] * f[i - 1])/btemp;
		btemp = b[i] - a[i] * temp[i];
		
	}
	v[N - 1] =  f[N - 1] / b[N - 1];
	for (int i = N-2; i>= 1; i--){
		v[i] -= temp[i + 1] * f[i + 1];
	}
	for (int i = 1; i < N; i++){
		cout << a[i] << " " << b[i] << " " << c[i] << " " << v[i] << " " << f[i] << endl;
	}
	delete []temp;
}

void solve_LU				(int N, TYPE *v, TYPE *u){
	mat A = zeros <mat>(N, N);
	vec f(N);
	mat L, U;
	vec X;
	vec temp;
	// initializing mat A:
	for (int i = 0; i < N; i++){
		A(i, i) = 2;
		if (i < N-1) {
			A(i, i+1) = -1;
			A(i+1, i) = -1;
		}
	}
	cout << A << endl;
	
	// initializing x : 
	
	TYPE h = 1.0/(N+1);
	for (int i = 1; i < N+1; i++) {
		v[i-1] = h*i;
		cout << v[i-1] << endl;
	}
	
	//initializing the value for function f:
	for (int i = 0; i < N; i++) f[i] = h * h * 100.0 * exp(-10.0 * v[i]);
	
	lu (L, U, A);
	X = solve (L, f);
	temp = solve (U, X);
	for (int i = 0; i < N; i++) {
		u[i] = temp(i);
		cout << u[i]<< endl;
	}
}


int main (){
	
	int 	N 	= 10;
	TYPE 	*a 	= new TYPE [N-1];
	TYPE 	*b 	= new TYPE [N];
	TYPE 	*c 	= new TYPE [N-1];
	TYPE 	*v	= new TYPE [N];
	TYPE 	*f	= new TYPE [N];
	TYPE	*f_LU	= new TYPE [N];
	initiate_abc (N, a, b, c, f);
	/*for (int i = 1; i <= N; i++){
		cout << a[i] << " " << b[i] << " " << c[i] << " " << x[i] << " " << f[i] << endl;
	}*/
	solve_matrix_equation_general	(N, v, f, a, b, c);
	solve_LU (N, v, f_LU);
	
	for (int i = 0; i < N; i++){
		cout << f[i] << " " << f_LU[i] << " " << get_u((TYPE)i) <<endl;
	}
	// Freeing memory :
	delete [] a; delete [] b; delete [] c; delete [] v; delete []f; delete []f_LU;
}


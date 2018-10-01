#include "JacobiMethod.h"
#include <stdlib.h>
#include <algorithm>

const TYPE epsilon = 10e-8;
const TYPE rho_min = 0;
const TYPE rho_max = 1;


/*
	Function that would find the initial value of a
*/
TYPE find_a (int N){
	TYPE rho_0, rho_N;

	rho_0 = rho_min;
	rho_N = rho_max;
	
	TYPE h = (rho_N - rho_0) / (1.0 * N);
	
	TYPE h2 = h*h;	

	TYPE a;
	
	a = -1.0/h2;
	return a;
}
/*
	Function that would find the initial value of d
*/
TYPE find_d (int N){
	TYPE rho_0, rho_N;

	rho_0 = rho_min;
	rho_N = rho_max;
	
	TYPE h = (rho_N - rho_0) / (1.0 * N);
	
	TYPE h2 = h*h;	

	TYPE d;
	
	d = 2.0/h2;
	return d;
}

void analytical_solution (int N, TYPE a, TYPE d, TYPE *lambda){
	for (int i = 1; i < N+1; i++){
		lambda [i-1] = d + 2*a*cos(i * M_PI / ((TYPE)N + 1));
	}
}

void jacobi_solution (int N, TYPE a, TYPE d, TYPE *lambda){

	TYPE **A = new TYPE* [N];
	TYPE **R = new TYPE* [N];
	for (int i = 0; i < N; i++){
		A[i] = new TYPE [N];
		R[i] = new TYPE [N];
	}
	
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (i == j){
				R[i][j] = 1;
				A[i][j] = d;
			}else if (abs(i - j) == 1){
				A[i][j] = a;
			}else {
				A[i][j] = 0;
				R[i][j] = 0;
			}
		}
	}
	

	int iter = jacobi_method(N, A, R);
	cout << iter << endl;
	
	for (int i = 0; i < N; i++)	lambda [i] = A[i][i];
	sort (lambda, lambda + N);
	
	delete [] A;delete [] R;
	//return 1;
}



int main (){
	int  N = 100;
	cout << "Enter value for N" << endl;
	cin >> N;
	TYPE *lambda_anal = new TYPE [N];
	TYPE *lambda_jaco = new TYPE [N];
	

	analytical_solution (N, find_a(N), find_d(N), lambda_anal);
	jacobi_solution (N, find_a(N), find_d(N), lambda_jaco);

	for (int i = 0; i < N; i++){
		printf (" %4.2f \t %4.2f \n", lambda_anal[i], lambda_jaco[i]);

	}
	
	return 0;
}


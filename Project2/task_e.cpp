#include "JacobiMethod.h"
#include <stdlib.h>
#include <algorithm>

const TYPE epsilon = 10e-8;
const TYPE rho_min = 0;
const TYPE rho_max = 5;

void printA(int N, TYPE **A){
	for (int i = 0; i < N; i++){
		cout << endl;
		for (int j = 0; j < N; j++){
			cout << A[i][j] << " ";
		}
	}
	cout << endl;
}

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
TYPE find_d (int N, int i, TYPE w_r){
	TYPE rho_0, rho_N, rho_i;
	TYPE V_i;
	TYPE d;
	rho_0 = rho_min;
	rho_N = rho_max;
	
	TYPE h = (rho_N - rho_0) / (1.0 * N);
	
	TYPE h2 = h*h;	

	rho_i = rho_0 + (i+1) * h;

	V_i = rho_i * rho_i;
	
	
	d = 2.0/h2 + w_r*w_r*V_i + 1.0/rho_i;
	return d;
}

void analytical_solution (int N, TYPE a, TYPE d, TYPE *lambda){
	for (int i = 1; i < N+1; i++){
		lambda [i-1] = d + 2*a*cos(i * M_PI / ((TYPE)N + 1));
	}
}

void jacobi_solution (int N, TYPE a, TYPE w_r, TYPE *lambda){

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
				A[i][j] = find_d(N, i, w_r);
			}else if (abs(i - j) == 1){
				A[i][j] = a;
			}else {
				A[i][j] = 0;
				R[i][j] = 0;
			}
		}
	}

	//printA (N, A);
	cout << "Hello!";

	int iter = jacobi_method(N, A, R);
	cout << iter << endl;
	
	for (int i = 0; i < N; i++)	lambda [i] = A[i][i];
	sort (lambda, lambda + N);
	
	delete [] A;delete [] R;
	//return 1;
}



int run_d (){
	int  N = 100;
	int  c = 4;
	TYPE omega_r[] = {0.01, 0.5, 1.0, 5.0};
	cout << "Enter value for N" << endl;
	
	cin >> N;
	
	TYPE **lambda_jaco = new TYPE* [c];
	for (int i = 0; i < 4; i++) lambda_jaco[i] = new TYPE [N];


	int p;
	for (int i = 0; i < 4; i++){
		//analytical_solution (N, find_a(N), find_d(N, 0), lambda_anal);
		
		cout << endl << rho_max << endl<<endl;
		jacobi_solution (N, find_a(N), omega_r[i], lambda_jaco[i]);

		
		
	}
	for (int j = 0; j < 10; j++){
			printf (" %4.2f \t %4.2f \t %4.2f \t %4.2f \t  \n", 
				lambda_jaco[0][j], lambda_jaco[1][j], lambda_jaco[2][j], lambda_jaco[3][j]);

	}
	return 0;
}

int main (){
	return 0;
}


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
void initiate				(int N, TYPE *a, TYPE *b , TYPE *c);
void solve_LU				(int N, TYPE v[], TYPE *u);
void solve_general			(int N, TYPE *v);
void solve_specific			(int N, TYPE *v);
TYPE find_max_error			(int N, TYPE v[], TYPE u[]);
string CPU_runtime			(void (*calculate) (int, TYPE *), int N, TYPE *v, TYPE u[], string method_name);


string CPU_runtime			(void (*calculate) (int, TYPE *), int N, TYPE *v, TYPE u[], string method_name){
	clock_t start, finish;
	start = clock();
	calculate (N, v);
	finish = clock ();
	TYPE e = find_max_error			(N, v, u);

	cout << "CPU runtime of method for N = "<< N << " "  << method_name << " " << ((finish - start)*1.0/CLOCKS_PER_SEC) << " " << e <<endl;
	TYPE took = ((finish - start)*1.0/CLOCKS_PER_SEC);
	string output =to_string (took) + " " + to_string(e);
	return output;
}


void initiate				(int N, TYPE *a, TYPE *b , TYPE *c){
	for (int j = 0; j < N; j++){
		a[j] = -1;
		c[j] = -1;
		b[j] = 2;
	}
}
void initiate_ft			(int N, TYPE *ft){
	TYPE h = 1.0 /(N + 1.0);

	TYPE	*x = new TYPE [N];
	TYPE	*f = new TYPE [N];

	for (int j = 0; j < N; j++){
		x[j] = h * (j + 1);
		f[j] = 100 * exp(-10 *x[j]);
		ft[j] = f[j] * h * h;
	}
	delete []x; delete[] f;
}

void solve_LU				(int N, TYPE *u){
	TYPE *v = new TYPE [N];
	initiate_ft (N, v);
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
	//cout << A << endl;
	
	// initializing x : 
	
	TYPE h = 1.0/(N+1);
	for (int i = 1; i < N+1; i++) {
		v[i-1] = h*i;
		//cout << v[i-1] << endl;
	}
	
	//initializing the value for function f:
	for (int i = 0; i < N; i++) f[i] = h * h * 100.0 * exp(-10.0 * v[i]);
	
	lu (L, U, A);
	X = solve (L, f);
	temp = solve (U, X);
	for (int i = 0; i < N; i++) {
		u[i] = temp(i);
		//cout << u[i]<< endl;
	}
	delete []v;
}


void solve_general 			(int N, TYPE *v){
	TYPE *f = new TYPE [N];
	
	TYPE *a = new TYPE[N];
	TYPE *b = new TYPE[N];
	TYPE *c = new TYPE[N];
	initiate (N, a, b, c);
	initiate_ft (N, f);
	for(int i=1; i<N; i++){
		b[i] -= c[i-1]*a[i-1]/b[i-1];
		f[i] -= f[i-1]*a[i-1]/b[i-1];
	}
	v[N-1] = f[N-1]/b[N-1];
	for(int i=N-2; i>=0; i--)
		v[i] = (f[i]-v[i+1]*c[i])/b[i];
		
	delete []a; delete []b; delete []c; delete []f;
}

void solve_specific 			(int N, TYPE *v){
	TYPE *f = new TYPE [N];

	TYPE *b = new TYPE [N];
	b [0] = 2.0;
	initiate_ft (N, f);
	// Forward

	for(int i = 1; i < N + 1; i++){
		b[i] = (i + 1.0) / (i * 1.0);
		f[i] +=  (i-1)*1.0 * f[i - 1] / (1.0 * i) ;
	}
	//b[N-1] = 2.0;
	v[N-1] = f[N - 1] / b[N - 1];
	// Backward
	for(int i = N - 2; i >= 0; i--){
		//v[i] = (f[i] + v[i+1]) / b[i];
		v[i] = i*1.0/(1.0*(i+1.0)) * (f[i] + v[i+1]);
	}
	delete []b; delete []f;
}


TYPE find_max_error 			 (int N, TYPE v[], TYPE u[]){
    TYPE max_error = -99999;
    for (int i = 1; i < N; i++){
	TYPE temp = log10 (abs(v[i] - u[i])/ u[i]);
	//cout << temp << endl;
	if (max_error < temp) max_error = temp;
    }
    return max_error;
}

void store (string file_name, string result){
	ofstream output (file_name);
	if (output.is_open())
		output << result;
	output.close();
}

void taskb (){
	int	N = 1;
	TYPE 	h;
	string r = "";
	for (int i = 1; i <= 3; i++){
		N *= 10;
		h = 1.0 /(N + 1.0);
		TYPE	*x = new TYPE [N];
		TYPE	*u = new TYPE [N];

		for (int j = 0; j < N; j++){
			x[j] = h * (j + 1);
			u[j] = 1 - (1-exp(-10))*x[j] - exp(-10*x[j]);
		}

		TYPE *v_general_solution = new TYPE[N];
		
		r +=  to_string(N) + " " + CPU_runtime(
			solve_general, N, v_general_solution, u, "General Solution"
		) + "\n";
		
		delete []x; delete []u;
		delete []v_general_solution;
	}
	string file_name = "taskb.txt";
	store (file_name, r);
}


void taskc (){
	int	N = 1;
	TYPE 	h;
	string r = "";
	for (int i = 1; i <= 6; i++){
		N *= 10;
		h = 1.0 /(N + 1.0);
		TYPE	*x = new TYPE [N];
		TYPE	*u = new TYPE [N];

		for (int j = 0; j < N; j++){
			x[j] = h * (j + 1);
			u[j] = 1 - (1-exp(-10))*x[j] - exp(-10*x[j]);
		}

		TYPE *v_general_solution = new TYPE[N];
		TYPE *v_s		 = new TYPE[N];
		r +=  to_string(N) + " " + CPU_runtime(
			solve_general, N, v_general_solution, u, "General Solution"
		) + " "
		+ CPU_runtime(
			solve_specific, N, v_s, u, "specific Solution"
		) + "\n";
		
		delete []x; delete []u;
		delete []v_general_solution;
		delete []v_s; 
	}
	string file_name = "taskc.txt";
	store (file_name, r);
}

void taskd (){
	int	N = 1;
	TYPE 	h;
	string r = "";
	for (int i = 1; i <= 4; i++){
		N *= 10;
		h = 1.0 /(N + 1.0);
		TYPE	*x = new TYPE [N];
		TYPE	*u = new TYPE [N];

		for (int j = 0; j < N; j++){
			x[j] = h * (j + 1);
			u[j] = 1 - (1-exp(-10))*x[j] - exp(-10*x[j]);
		}

		TYPE *v_general_solution = new TYPE[N];
		TYPE *v_LU		 = new TYPE[N];
		r +=  to_string(N) + " " + CPU_runtime(
			solve_general, N, v_general_solution, u, "General Solution"
		) + " "
		+ CPU_runtime(
			solve_LU, N, v_LU, u, "LU Solution"
		) + "\n";
		
		delete []x; delete []u;
		delete []v_general_solution;
		delete []v_LU; 
	}
	string file_name = "taskd.txt";
	store (file_name, r);
}

int main (){
	//int 	k = 1;
/*	int	N = 1;
	TYPE 	h;
	
	for (int i = 1; i <= 7; i++){
		N *= 10;
		h = 1.0 /(N + 1.0);
		TYPE	*x = new TYPE [N];
		//TYPE	*f = new TYPE [N];
		TYPE	*u = new TYPE [N];
		//TYPE	*ft = new TYPE [N];

		for (int j = 0; j < N; j++){
			x[j] = h * (j + 1);
			//f[j] = 100 * exp(-10 *x[j]);
			u[j] = 1 - (1-exp(-10))*x[j] - exp(-10*x[j]);
			//ft[j] = f[j] * h * h;
		}

		TYPE *v_general_solution = new TYPE[N];
		TYPE *v_LU = new TYPE[N];
		TYPE *v_s = new TYPE[N];
		//solve_specific	(N, v_s);
		cout << CPU_runtime(
			solve_general, N, v_general_solution, u, "General Solution"
		);
		/*cout << CPU_runtime(
			solve_LU, N, v_LU, u, "LU Solution"
		);*/
		
		/*
		cout << "\n\n\n\n\n";
		//for (int j = 0; j < N; j++)
			//cout << v_general_solution[j] << "\t\t" << v_s[j] << "\t\t" << v_LU[j] << "\t\t" << u[j] << '\n'; 
		delete []x; delete []u;
		delete []v_general_solution;
		delete []v_LU;
		delete []v_s;
	}*/
	taskb();
	taskc();
	taskd();
}

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
TYPE RHO_MIN =	0.0;
TYPE RHO_MAX =	1.0;


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



  /*
   * The function                             
   *      TYPE find_a (int N)            
   * Find the initial values of a, for the matrix 
   * and initiates the matrix to be I(N*M)
   * 
   * Input data:                      
   *  int N		- number of steps          
   *
   * Returns the value of a based on N
   */

TYPE find_a (int N){
	TYPE rho_0, rho_N;

	rho_0 = RHO_MIN;
	rho_N = RHO_MAX;
	
	TYPE h = (rho_N - rho_0) / (1.0 * N);
	
	TYPE h2 = h*h;	

	TYPE a;
	
	a = -1.0/h2;
	return a;
}


  /*
   * The function                             
   *      TYPE find_d (int N, int i, int task, TYPE w_r = 0)
   * Find the initial values of d (e_i), for the matrix 
   * and initiates the matrix to be I(N*M)
   * 
   * Input data:                      
   *  int N		- number of steps          
   *  int i		- for task d (to find the potential)  
   *  int task		- to add or not to add!        
   *  TYPE w_r		- for task e to include the omega factor          
   *
   * Returns the value of d based on N, rho_i, and w_r
   */

TYPE find_d (int N, int i= 0, int task = 0, TYPE w_r = 0){
	TYPE rho_0, rho_N, rho_i;
	TYPE V_i;
	TYPE d = 0;
	rho_0 = RHO_MIN;
	rho_N = RHO_MAX;
	
	TYPE h = (rho_N - rho_0) / (1.0 * N);
	
	TYPE h2 = h*h;	

	rho_i = rho_0 + (i+1) * h;

	V_i = rho_i * rho_i;
	d = 2.0/h2;
	if (task == 1) d += V_i;
	if (task == 2) d = w_r*w_r*V_i + 1.0/rho_i;
	
	return d;
}



  /*
   * The function                             
   *      TYPE **allocateMemory (int N, int M, bool I = false)                  
   * allocates the memory for a matrix of size N*M 
   * and initiates the matrix to be I(N*M)
   * 
   * Input data:                      
   *  int N		- number of rows          
   *  int M		- number of cols
   *  bool I		- true : I (N*N)
   * Returns the pointer to the array
   */



TYPE **allocateMemory (int N, int M, bool I = false){
	TYPE **A = new TYPE *[N];
		
	for (int i = 0; i < N; i++){
		A[i] = new TYPE[M];
		for (int j = 0; j < M; j++){ 
			A[i][j] = 0;
			if (I == true &&  i == j) A[i][j] = 1;
		}
	}

	return A;
}

  /*
   * The function                             
   *      void deallocateMemory (int N, int M, TYPE **A)
   * deallocates the memory of the matrix of size N*M 
   * and initiates the matrix to be I(N*M)
   * 
   * Input data:                      
   *  int N		- number of rows          
   *  int M		- number of cols
   *  TYPE **A		- The matrix
   * Returns void!
   */
void deallocateMemory (int N, int M, TYPE **A){
	for (int i = 0; i < N; i++)
		delete []A[i];
}


  /*
   * The function                             
   *      TYPE max_off_diagonal (int N, int *x, int *y, TYPE **A)                    
   * finds the element with the largest |e_ij| in the matrix
   * where i != j, considering that the matrix is symmetric, 
   * just checks for the upper triangle
   * Input data:                      
   *  int N		- number of  rows          
   *  int &x		- returns x index of the element  (God bless pointers!)       
   *  int &y		- returns y index of the element  (God bless pointers!)
   *  TYPE **A 		- pointer to the matrix                 
   * Returns the value of the largest element                                
   */
TYPE max_off_diagonal (int N, int *x, int *y, TYPE **A){
	TYPE element_max= 0.0;
	for (int i = 0; i < N; i++){
		for (int j = i + 1; j < N; j++){
			if (fabs(A[i][j]) > element_max){
				*x = i;
				*y = j;
				element_max = fabs(A[i][j]);
			}
		}
	}
	return element_max;
}

  /*
   * The function                             
   *      void jacobi_rotate (int N, int l, int k, TYPE **A, TYPE **R)                   
   * Executes one jacobi rotation on the matrix A,
   * 
   * Input data:                      
   *  int N		- number of  rows
   *  int l		- x index of the max_off_diagonal element
   *  int k		- y index of the max_off_diagonal element
   *  int &x		- returns x index of the element  (God bless pointers!)       
   *  TYPE **A		- pointer to the matrix A, the rotation will be on matrix A
   *  TYPE **R 		- pointer to the matrix R, initially R is I_{n*n}                
   * Returns void! The Eigenvalues and eigenvectors will be kept on R and A                                
   */
void jacobi_rotate (int N, int l, int k, TYPE **A, TYPE **R){
	TYPE tau;
	TYPE s, c, t;


	if (A[k][l] != 0) {
		tau = (A[l][l] - A[k][k])/(2.0 * A[k][l]);
		t = -tau - sqrt(1 + tau * tau);
		c = 1.0 / sqrt(1 + t * t);
		s = t * c;
	} else {
		c = 1.0;
		s = 0.0;
	}

	TYPE A_ik, A_kk;
	TYPE r_ik, r_il;
	for (int i = 0; i < N; i++) {
		if (i != k && i != l) {
			A_ik = A[i][k];
			A[i][k] = A_ik * c - A[i][l] * s;
			A[k][i] = A[i][k];
			A[i][l] = A[i][l] * c + A_ik * s;
			A[l][i] = A[i][l];
		}
		r_ik = R[i][k];
		r_il = R[i][l];
		R[i][k] = c * r_ik - s * r_il;
		R[i][l] = c * r_il + s * r_ik;
	}

	A_kk = A[k][k];
	A[k][k] = A_kk * c * c - 2 * A[k][l] * c * s + A[l][l] * s * s;
	A[l][l] = A[l][l] * c * c + 2 * A[k][l] * c * s + A_kk * s * s;
	A[k][l] = 0;
	A[l][k] = 0;

}

  /*
   * The function                             
   *      int jacobi_method(int N, TYPE **A, TYPE **R)                
   * Executes one jacobi rotation on the matrix A,
   * 
   * Input data:                      
   *  int N		- number of  rows
   *  int l		- x index of the max_off_diagonal element
   *  int k		- y index of the max_off_diagonal element
   *  int &x		- returns x index of the element  (God bless pointers!)       
   *  TYPE **A		- pointer to the matrix A, the rotation will be on matrix A
   *  TYPE **R 		- pointer to the matrix R, initially R is I_{n*n}                
   * Returns void! The Eigenvalues and eigenvectors will be kept on R and A                                
   */

int jacobi_method(int N, TYPE **A, TYPE **R){
	int iter = 0;
	int i, j;
	int max_iter = N * N * N;

	TYPE tolerance = 1.0e-8;
	TYPE max_off_diag_value = 0;
	max_off_diag_value = max_off_diagonal(N, &i, &j, A);



	while (max_off_diag_value > tolerance && iter < max_iter){
		iter ++;
		jacobi_rotate (N, i, j, A, R);
		max_off_diag_value = max_off_diagonal (N, &i, &j, A);
	}
	// To know how many iteration it took
	return iter;
}


/* ************************************ */
/* ARMADILLOOOOOOOOOOOOOOOOOOOOOOOOOOOO */
/* ************************************ */



mat generate_mat (int N, TYPE h, TYPE a, TYPE d){
	mat  A = zeros <mat>(N, N);
	for (int i = 0; i < N; i++){
		A(i, i) = d;
		if (i < N-1) {
			A(i, i+1) = a;
			A(i+1, i) = a;
		}
	}
	return A; 
}


void arma_diag (Eig *eig){
	eig->S = eig->eigvector;
	eig->D = diagmat(eig->eigvalue);
}

Eig arma_eig (int N, TYPE a, TYPE d){
	TYPE h = 1.0/(1.0+ N);
	mat A = generate_mat (N, h, a, d);
	Eig result;
		
	eig_sym(result.eigvalue, result.eigvector, A); 
	//cout << diagmat(result.eigvalue) << endl << result.eigvector << endl;
	return result;
}



int run_jacob (int N){
	TYPE **A = allocateMemory (N, N);
	TYPE a = find_a (N);
	TYPE d = find_d (N);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (i == j){
				A[i][j] = find_d(N, i);
			}else if (abs(i - j) == 1){
				A[i][j] = a;
			}
		}
	}
	TYPE **R = allocateMemory (N, N, true);
	int iter = jacobi_method(N, A, R);
	return iter;
}


void run_arma (int N){
	TYPE a = find_a (N);
	TYPE d = find_d (N);
	Eig result = arma_eig (N, a, d);
	arma_diag(&result);
}

void jacob_vs_arma_stat (){
	int N = 600;
	//cin >> N;
	ofstream fout ("jacob_vs_arma.txt");
	for (int i = 100; i < N; i+= 100){
		cout << i << endl;
		clock_t start_j, finish_j, start_a, finish_a;
		start_j = clock();
		run_jacob (i);
		finish_j = clock ();

		start_a = clock();
		run_arma (i);
		finish_a = clock ();
		
		TYPE time_jacob = ((finish_j - start_j)*1.0/CLOCKS_PER_SEC);

		TYPE time_arma = ((finish_a - start_a)*1.0/CLOCKS_PER_SEC);
		string output = to_string (i) + " " + to_string (time_jacob) + " " + to_string(time_arma) + "\n";	
		cout << output << endl;
		fout << output;
	}
	fout.close ();
}

void jacob_iter (){
	int N = 600;
	//cin >> N;
	ofstream fout ("jacob_iter.txt");
	for (int i = 100; i < N; i+= 100){
		cout << i << endl;
		
		int iter = run_jacob (i);
		
		string output = to_string (i) + " " + to_string (iter) + "\n";	
		cout << output << endl;
		fout << output;
	}
	fout.close ();
}


bool test_max_off_diagonal() {
	int N = 3;
	TYPE **A = new TYPE* [N];
	for (int i = 0; i < N; i++){
		A[i] = new TYPE [N];
	}
	A[0][0] = 2;
	A[0][1] = 3;
	A[0][2] = 1;

	A[1][0] = 3;
	A[1][1] = 5;
	A[1][2] = 0;

	A[2][0] = 1;
	A[2][1] = 0;
	A[2][2] = -7;

	/*
	A[0] = {2, 3, 1};
	A[1] = {3, 5, 0};
	A[2] = {1, 0, -7};*/
	int x, y;
	TYPE diag = max_off_diagonal(3, &x, &y, A);
	cout << x << " " << y << endl;
	TYPE test_diag = 0.0;
	delete [] A;
	return  fabs(diag - 3) < EPSILON ;
}



bool test_Jacobi_rotation() {
	int N = 3;
	TYPE **A = new TYPE* [3];
	TYPE **R = new TYPE* [3];
	for (int i = 0; i < 3; i++){
		A[i] = new TYPE [N];
		R[i] = new TYPE [N];
	}
	A[0][0] = 2; A[0][1] = 3; A[0][2] = 1;
	A[1][0] = 3; A[1][1] = 5; A[1][2] = 0;
	A[2][0] = 1; A[2][1] = 0; A[2][2] = -7;

	R[0][0] = 1; R[0][1] = 0; R[0][2] = 0;
	R[1][0] = 0; R[1][1] = 1; R[1][2] = 0;
	R[2][0] = 0; R[2][1] = 0; R[2][2] = 1;

	int x, y;
	x = 0; y = 1;
	int iter = jacobi_method(N, A, R);
	cout << iter << endl;
	
	for (int i = 0; i < N; i ++){
		printf (" %4.2f \t %4.2f \t %4.2f \n", A[i][0], A[i][1], A[i][2]);
	}
	
	cout << endl;

	for (int i = 0; i < N; i ++){
		printf (" %4.2f \t %4.2f \t %4.2f \n", R[i][0], R[i][1], R[i][2]);
	}

	TYPE test_diag = 0.0;

	
	delete [] A;
	return 1;
}



void printA(int N, TYPE **A){
	for (int i = 0; i < N; i++){
		cout << endl;
		for (int j = 0; j < N; j++){
			cout << A[i][j] << " ";
		}
	}
	cout << endl;
}


void analytical_solution (int N, TYPE a, TYPE d, TYPE *lambda){
	for (int i = 1; i < N+1; i++){
		lambda [i-1] = d + 2*a*cos(i * M_PI / ((TYPE)N + 1));
	}
}

void jacobi_solution_d (int N, TYPE a, TYPE d, TYPE *lambda){

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
				A[i][j] = find_d(N, i, 1);
			}else if (abs(i - j) == 1){
				A[i][j] = a;
			}else {
				A[i][j] = 0;
				R[i][j] = 0;
			}
		}
	}

	//printA (N, A);

	int iter = jacobi_method(N, A, R);
	cout << iter << endl;
	
	for (int i = 0; i < N; i++)	lambda [i] = A[i][i];
	sort (lambda, lambda + N);
	
	delete [] A;delete [] R;
	//return 1;
}


int run_task_d (){
	int  N = 200;
	TYPE rho_max = 10;
	cout << "Enter value for N" << endl;
	//cin >> N;
	TYPE *lambda_anal = new TYPE [N];
	TYPE *lambda_jaco = new TYPE [N];
	int p;
	for (int i = 3; i < 7; i += 2){
		//analytical_solution (N, find_a(N), find_d(N, 0), lambda_anal);
		RHO_MAX = i*1.0;
		cout << endl << RHO_MAX << endl<<endl;
		jacobi_solution_d (N, find_a(N), find_d(N, i, 1), lambda_jaco);

		for (int j = 0; j < 10; j++){
			printf (" %4.2f  \n", lambda_jaco[j]);

		}
		
	}
	return 0;
}

void jacobi_solution_e (int N, TYPE a, TYPE w_r, TYPE *lambda){

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
				A[i][j] = (N, i, 2, w_r);
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



int run_task_e (){
	int  N = 100;
	int  c = 4;
	RHO_MAX = 5.0;
	TYPE omega_r[] = {0.01, 0.5, 1.0, 5.0};
	cout << "Enter value for N" << endl;
	
	//cin >> N;
	
	TYPE **lambda_jaco = new TYPE* [c];
	for (int i = 0; i < 4; i++) lambda_jaco[i] = new TYPE [N];


	int p;
	for (int i = 0; i < 4; i++){
		cout << endl << RHO_MAX - RHO_MIN << endl<<endl;
		jacobi_solution_e (N, find_a(N), omega_r[i], lambda_jaco[i]);

		
		
	}
	ofstream fout ("jacob_task_e.txt");
	for (int j = 0; j < 10; j++){
		fout << lambda_jaco[0][j] <<" " << lambda_jaco[1][j] << " " << lambda_jaco[2][j]<< " " << lambda_jaco[3][j] << endl;
	}
	fout.close();
	return 0;
}


int main (){
	
	//jacob_vs_arma_stat();
	//jacob_iter();

	//run_task_d();
	run_task_e();
	return 0;
}




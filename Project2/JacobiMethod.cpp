#include "JacobiMethod.h"

/*struct Eig {
	vec	eigvalue;
	mat	eigvector;
	mat	S;
	mat	D;
};*/




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
int main (){
	int n;
	cin >> n;
	
	Eig result = arma_eig (n, -1, 2);
	arma_diag(&result);
	cout << result.D << endl;
	return 0;
}



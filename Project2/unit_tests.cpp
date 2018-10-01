#include "JacobiMethod.h"

const TYPE epsilon = 10e-8;

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
	return (fabs(diag - 3) < epsilon);
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

int main (){
	cout << "result " << test_max_off_diagonal() << endl;
	cout << test_Jacobi_rotation ();	
	return 0;
}


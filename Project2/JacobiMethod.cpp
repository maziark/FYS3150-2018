#include "JacobiMethod.h"

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



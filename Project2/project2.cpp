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
struct Eig {
	vec	eigvalue;
	mat	eigvector;
	mat	S;
	mat	D;
};

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
	eig.S = eig.eigvector;
	eig.D = diagmat(eig.eigenvalue);
}

Eig arma_eig (int N, TYPE a, TYPE d){
	TYPE h = 1.0/(1.0+ N);
	mat A = generate_mat (N, h, a, d);
	Eig result;
		
	eig_sym(result.eigvalue, result.eigvector, A); 
	cout << diagmat(result.eigvalue) << endl << result.eigvector << endl;
	return result;
}
int main (){
	int n;
	cin >> n;
	
	Eig *result = arma_eig (n, -1, 2);
	arma_diag(result);
	cout << result->D << endl;
	return 0;
}


void generate_A (){
	
}


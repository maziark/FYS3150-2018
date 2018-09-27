#include "JacobiMethod.h"

const TYPE epsilon = 10e-8;


<<<<<<< HEAD
bool test_max_off_diagonal() {
  int N = 3;
  TYPE **A = new TYPE* [3];


  for (int i = 0; i < 3; i++){
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

  
  int x, y;
  TYPE diag = max_off_diagonal(3, &x, &y, A);
  cout << x << " " << y << endl;
  TYPE test_diag = 0.0;
  delete [] A;
  return (fabs(diag - 3) < epsilon);
=======
bool test_max_off_diagonal_3() {
	int N = 3;
	TYPE **A = new TYPE* [3];
	for (int i = 0; i < 3; i++){
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
>>>>>>> e688164b9378f24a1376fadfea9c86e92e752ab4

}

int main (){
<<<<<<< HEAD
  cout << "result " << test_max_off_diagonal() << endl;
  return 0;
}



=======
	cout << "result " << test_max_off_diagonal() << endl;
	return 0;
}
>>>>>>> e688164b9378f24a1376fadfea9c86e92e752ab4

#include "time.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>

using namespace std;
// TODO : File I/O initialization here!
ofstream ofile;

random_device rd;
mt19937 engine(rd());
uniform_real_distribution <double> randd(0.0, 1.0);

//Random numbers:
double ran2(){
    return randd(engine);
}


// Function calls
void initialize(int n_spins, int **spin_matrix, double& E, double& M, double tmp);
void metropolis(int n_spins, int **spin_matrix, double& E, double& M, double *w, int& accepted, vector <int>& accepted_vec, int& count);
void **matrix(int row, int col, int num_bytes);
void free_matrix(void **matr);

// inline function for Periodic Boundary Conditions (PBC)
inline int periodic(int i, int limit, int add) { 
  return (i + limit + add) % (limit);
}




int main (int argc, char* argv[]){
	// init variables;
	// MPI Config variables : 	
	int my_rank;
	int numprocs;
	
	char *outfilename, name;
	long idum;
	int **spin_matrix, n_spins, mcs, n;
	double w[17], average[5], total_average[5], 
	initial_temp, final_temp, E, M, temp_step;



	// init MPI
	n_spins = 20; mcs = 1000000;  initial_temp = 2.1; final_temp = 2.4; temp_step =0.02;	
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	int no_intervalls = mcs/numprocs;
	int myloop_begin = my_rank*no_intervalls + 1;
	int myloop_end = (my_rank+1)*no_intervalls;
	if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

	// broadcast to all nodes common variables
	MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*MPI_Get_processor_name(
	    &name,
	    &n);
	cout << name << " " << numprocs << " " << my_rank << endl;*/
	
	spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
	if (my_rank == 0 && argc <= 1) {
		cout << "Bad Usage: " << argv[0] << " read output file" << endl;
		exit(1);
	}
	if (my_rank == 0 && argc > 1) {
		outfilename=argv[1];
		ofile.open(outfilename); 
	}

	// Start monto Carlo :
	double  TimeStart, TimeEnd, TotalTime;
	TimeStart = MPI_Wtime();

	vector <double> E_vec(mcs);
	vector <double> absM_vec(mcs);
	vector <double> E2_vec(mcs);
	vector <double> absM2_vec(mcs);
	vector <double> num_susceptibility_vec(mcs);
	vector <double> num_heatCapacity_vec(mcs);
	int accepted = 0, count = 0;
	vector <int> accepted_vec(mcs);

	for (double temp = initial_temp; temp <= final_temp; temp+= temp_step){
		E = M = 0;
		// Initializing spin matrix
		initialize(n_spins, spin_matrix, E, M, temp);
		// Setting all variables to their initial values
		for( int de = -8; de <= 8; de++) w[de + 8] = 0;
		for( int de = -8; de <= 8; de += 4) w[de + 8] = exp(-de / temp);
		for( int i = 0; i < 5; i++) average[i] = 0.;
		for( int i = 0; i < 5; i++) total_average[i] = 0.;

		// Start Monto Carlo Calculations

		for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
			metropolis(n_spins, spin_matrix, E, M, w, accepted, accepted_vec, count);
			// Updating the average values;
			average[0] += E;    average[1] += E*E;
			average[2] += M;    average[3] += M*M; average[4] += fabs(M);
			
			// Keeps the values for plotting porpuses

			E_vec[cycles] = E;
			E2_vec[cycles] = E*E;
			absM_vec[cycles] = fabs(M);
			absM2_vec[cycles] = fabs(M) * fabs(M);
		}

		// Find total average : 
		for( int i =0; i < 5; i++)
			MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		// print results for each set of cycles
		
		// TODO : output function is not defined yet
		/*if (my_rank == 0)
			output(n_spins, mcs, temp, average, total_average, E_vec, E2_vec, absM_vec, absM2_vec, num_susceptibility_vec, num_heatCapacity_vec);*/
	
	}

	free_matrix((void **) spin_matrix); // free memory

	TimeEnd = MPI_Wtime();
	TotalTime = TimeEnd-TimeStart;
	if ( my_rank == 0) 
		cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
	

	// END MPI
	MPI_Finalize (); 
	// dealocate memory
	
	// save to File


	cout << "Hello World!" << endl; 


	return 0;
}


// Function definisions : 
void initialize(int n_spins, int **spin_matrix, double& E, double& M, double tmp){
	// Initializing the spin matrix
	for (int i = 0; i < n_spins; i++){
		for (int j = 0; j < n_spins; j++){
			double chance= ran2();
			if (tmp>1.5) chance = 0;
			spin_matrix[i][j] = (chance >= 0.5)? -1 : 1;
			M += (double)spin_matrix[i][j];
		}
	}

	// Calculating Magnetism
	for (int i = 0; i < n_spins; i++)
		for (int j = 0; j < n_spins; j++)
			E -= (double) spin_matrix[i][j] * 
				(spin_matrix[periodic(i, n_spins, -1)][j] + spin_matrix[i][periodic(j, n_spins, -1)]);

}

void metropolis(int n_spins, int **spin_matrix, double& E, double& M, double *w, int& accepted, vector <int>& accepted_vec, int& count){
	for (int i = 0; i < n_spins; i++){
		for (int j = 0; j < n_spins; j++){
			// Selecting a random element of the matrix
			int x = (int) (ran2() * (double)n_spins);
			int y = (int) (ran2() * (double)n_spins);

			// Calculating Energy difference : 

			int deltaE =  2 * spin_matrix[y][x] *
				(spin_matrix[y][periodic(x, n_spins, -1)] +
					spin_matrix[periodic(y, n_spins, -1)][x] +
					spin_matrix[y][periodic(x , n_spins, 1)] +
					spin_matrix[periodic(y, n_spins, 1)][x]);

			// Metropolis test
			if (ran2() <= w[8 + deltaE]){
				spin_matrix [y][x] *= -1; // Change the spin
				M += (double)(2 * spin_matrix[y][x]); // (1 - (-1)) = 2
				E += (double)deltaE;
				accepted ++;
			}
		}
	}
	// Measure of efficiency of the algorithm
	accepted_vec[count] = accepted;
	count ++;
}



void **matrix(int row, int col, int num_bytes){
	int      i, num;
	char     **pointer, *ptr;

	pointer = new(nothrow) char* [row];
	if(!pointer) {
		cout << "Exception handling: Memory allocation failed";
		cout << " for "<< row << "row addresses !" << endl;
		return NULL;
	}
	i = (row * col * num_bytes)/sizeof(char);
	pointer[0] = new(nothrow) char [i];
	if(!pointer[0]) {
		cout << "Exception handling: Memory allocation failed";
		cout << " for address to " << i << " characters !" << endl;
		return NULL;
	}
	ptr = pointer[0];
	num = col * num_bytes;
	for(i = 0; i < row; i++, ptr += num )   {
		pointer[i] = ptr;
	}

	return  (void **)pointer;

} // end: function void **matrix()

    /*
     * The function
     *      void free_matrix()
     * releases the memory reserved by the function matrix()
     *for the two-dimensional matrix[][]
     * Input data:
     *  void far **matr - pointer to the matrix
     */

void free_matrix(void **matr)
{

	delete [] (char *) matr[0];
	delete [] matr;

}  // End:  function free_matrix()


//Function that calculates analytic average energy
double analyticEavg(double Z, double beta){
    double analytic_Eavg = (16*exp(-8*beta) - 16*exp(8*beta))/Z;
    return analytic_Eavg; //Skal bli ca -2
}

//Function that calculates analytic average magnetization
double analyticMavg(double Z, double beta){
    double analytic_Mavg = (8*exp(8*beta) + 16)/Z;
    return analytic_Mavg; //Skal bli ca
}

//Function that calculates analytic average energy squared
double analyticE2avg(double Z, double beta){
    double analytic_E2avg = 256*cosh(8*beta)/Z;
    return analytic_E2avg;
}

//Function that calculates analytic average magnetization squared:
double analyticM2avg(double Z, double beta){
    double analytic_M2avg = (32 + 32*exp(8*beta))/Z;
    return analytic_M2avg;
}

double analyticSusceptibility(double analytic_Mavg, double analytic_M2avg, double temp){
    //double analytical_Susceptibility = 1/temp*( (32 + +32*exp(8*beta))/Z - pow( 16 + 8*exp(8*beta), 2 )/(Z*Z) );
    double analytic_Susceptibility = (analytic_M2avg - analytic_Mavg*analytic_Mavg)/temp;
    return analytic_Susceptibility;
}

double analyticHeatCapacity(double analytic_Eavg, double analytic_E2avg, double temp){
    double analytic_Heat_Capacity = (analytic_E2avg - analytic_Eavg*analytic_Eavg)/(temp*temp);
    return analytic_Heat_Capacity;
}

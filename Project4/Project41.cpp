/*
   Program to solve the two-dimensional Ising model
   with zero external field using MPI
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
*/
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <random>
#include <time.h>
using namespace  std;

random_device rd;
mt19937 engine(rd());
uniform_real_distribution <double> randd(0.0, 1.0);

// output file
ofstream ofile;
ofstream confFile;
ofstream cycleOut;



// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i + limit + add) % (limit);
}


double analyticEavg(double, double);
double analyticE2avg(double, double);
double analyticMavg(double, double);
double analyticM2avg(double, double);
double analyticSusceptibility(double, double, double);
double analyticHeatCapacity(double, double, double);

//Random numbers:
double rrandom(){
    return randd(engine);
}




// Function to initialise energy and magnetization
void initialize(int, int **, double&, double&);
void initialize(int, int **, double&, double&, int);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);
void metropolis(int, long&, int **, double& , double& , double *, int& , vector <int>& , int& );
// prints to file the results of the calculations
void output(int, int, double, double *);
//  Matrix memory allocation
//  allocate space for a matrix
void  **matrix(int, int, int);
//  free space for  a matrix
void free_matrix(void **);
// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);

// Main program begins here

int main(int argc, char* argv[])
{
    confFile.open("conf.txt");
    cycleOut.open("cycleStat.txt");
    char *outfilename;
    
    long idum;
    int **spin_matrix, n_spins, mcs, my_rank, numprocs;
    double w[17], average[5], total_average[5], initial_temp, final_temp, E, M, temp_step;

    //  MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    if (my_rank == 0 && argc <= 1) {
        cout << "Bad Usage: " << argv[0] <<
        " read output file" << endl;
        exit(1);
    }

    if (my_rank == 0 && argc > 1) {
        outfilename=argv[1];
        ofile.open(outfilename);
    }
    n_spins = 2; mcs = 100;  initial_temp =1.0; final_temp = 1.0; temp_step =0.0;
    if (my_rank == 0 && argc > 2)
        n_spins=atoi(argv[2]);

    /*
    Determine number of intervall which are used by all processes
    myloop_begin gives the starting point on process my_rank
    myloop_end gives the end point for summation on process my_rank
    */

    int no_intervalls = mcs/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

    vector <double> E_vec(mcs);
    vector <double> absM_vec(mcs);
    vector <double> E2_vec(mcs);
    vector <double> absM2_vec(mcs);
    vector <double> num_susceptibility_vec(mcs);
    vector <double> num_heatCapacity_vec(mcs);

    vector <int> accepted_vec(mcs);
    int accepted = 0;
    int mc_count = 0;
    double k = 1; // We assume that k = 1, for this assignment


    // broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&mc_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&accepted, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //  Allocate memory for spin matrix
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    // every node has its own seed for the random numbers, this is important otherwise
    // if one starts with the same seed, one ends with the same random numbers
    idum = -1-my_rank;  // random starting point
    // Start Monte Carlo sampling by looping over T first

    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
    double temperature= initial_temp;
    double  beta, Z, aEavg, aE2avg, aMavg, aM2avg, aS, aHC;
    //for ( temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){}
	// Produce the analyticValues : 
	beta = 1.0/ (k * temperature);
	if (n_spins == 2) {
		Z = 12 + 2*exp(8*beta) + 2*exp(-8*beta);
		aEavg = analyticEavg(Z, beta);
		aE2avg = analyticE2avg(Z, beta);
		aMavg = analyticMavg(Z, beta);
		aM2avg = analyticM2avg(Z, beta);
		aS = analyticSusceptibility(aMavg, aM2avg, temperature);
		aHC = analyticHeatCapacity(aEavg, aE2avg, temperature);
		double aEvariance = (aE2avg- aEavg * aEavg) /4;
		double aMvariance = (aM2avg - aMavg*aMavg)/4;
		cout << "Temperature : " <<temperature << endl;
		cout << "Average_E : " << aEavg/4 << endl;
		cout << "Variance_E : " << aEvariance/temperature/temperature << endl;
		cout << "|M| : " << aMavg/4 << endl;
		cout << "variance_|M| : " << aMvariance /4 << endl;
		cout << "C_v : " << aHC << endl;
		cout << "Susceptibility : " << aS << endl;
		
	}

		
	//    initialise energy and magnetization
        E = M = 0.;
        // initialise array for expectation values
        initialize(n_spins, spin_matrix, E, M, 0);
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
        for( int i = 0; i < 5; i++) average[i] = 0.;
        for( int i = 0; i < 5; i++) total_average[i] = 0.;
        // start Monte Carlo computation
	cout << myloop_begin << " " << myloop_end << endl;
        for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
	    // For Config statistics
            metropolis(n_spins, idum, spin_matrix, E, M, w, accepted, accepted_vec, mc_count);
	    
	    // For results on range of tempertures
	    Metropolis(n_spins, idum, spin_matrix, E, M, w);
            // update expectation values  for local node
            average[0] += E;    average[1] += E*E;
            average[2] += M;    average[3] += M*M; average[4] += fabs(M);
	    
            E_vec[cycles] = average[0];  
            E2_vec[cycles] = E*E;
            absM_vec[cycles] = fabs(M);
            absM2_vec[cycles] = fabs(M) * fabs(M); 
        }
        // Find total average
        for( int i =0; i < 5; i++)
            MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        // print results
        if ( my_rank == 0)
            output(n_spins, mcs, temperature, total_average);

	if (my_rank == 0){
	    //confFile << "cycle" << " " << "AcceptedConfig" << " " << "E" << " " << "|M|" <<  " " << mcs<< endl;
            for (int i=0; i < mcs; i++){
		//confFile << i << " " << accepted_vec[i] << " " << temperature << "  " << E_vec[i]/(n_spins*n_spins) << "  " << E2_vec[i]/(n_spins*n_spins) << "  " << absM_vec[i]/(n_spins*n_spins) << "  " << absM2_vec[i]/(n_spins*n_spins) << "  " << num_heatCapacity_vec[i]/(n_spins*n_spins) << "  " << num_susceptibility_vec[i]/(n_spins*n_spins) << endl;
                confFile << i << " " << accepted_vec[i] << " " 
			<< temperature << "  " 
			<< E_vec[i] << "  " 
			<< absM_vec[i]/(n_spins*n_spins) << " " 
			<< num_heatCapacity_vec[i]/(n_spins*n_spins) << "  "
			<< num_susceptibility_vec[i]/(n_spins*n_spins) << endl;
	     }
		
        }

    //} // NOTE : Take this one out to close the loop! *** Config statistics ***
    free_matrix((void **) spin_matrix); // free memory
    ofile.close();  // close output file
    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0)
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;


    // End MPI
    MPI_Finalize ();
    return 0;
}

// function to initialize energy, spin matrix and magnetization
void initialize(int n_spins, int **spin_matrix,
        double& E, double& M)
{
// setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            spin_matrix[y][x] = 1; // spin orientation for the ground state
            M +=  (double) spin_matrix[y][x];
        }
    }
  // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*
            (spin_matrix[periodic(y,n_spins,-1)][x] +
            spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise


void initialize(int n_spins, int **spin_matrix,
        double& E, double& M, int randomness)
{
srand(time(0));
// setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            if (randomness > 0) {
		int r = rrandom();//((rand()%2)==0) ? -1:1; 
	    	spin_matrix[y][x] = (r>0.5)? 1:-1;
	    }else spin_matrix[y][x] = 1; // spin orientation for the ground state
            M +=  (double) spin_matrix[y][x];
        }
    }
  // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*
            (spin_matrix[periodic(y,n_spins,-1)][x] +
            spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise


void metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double& M, double *w, int& accepted, vector <int>& accepted_vec, int& mc_count)

{
    // loop over all spins
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            //int ix = (int) (ran2(&idum)*(double)n_spins);
            //int iy = (int) (ran2(&idum)*(double)n_spins);
	    int ix = floor(rrandom() * n_spins); //Random x and y which means a random spin is picked
            int iy = floor(rrandom() * n_spins);

            int deltaE =  2*spin_matrix[iy][ix]*
            (spin_matrix[iy][periodic(ix,n_spins,-1)]+
            spin_matrix[periodic(iy,n_spins,-1)][ix] +
            spin_matrix[iy][periodic(ix,n_spins,1)] +
            spin_matrix[periodic(iy,n_spins,1)][ix]);
            if (rrandom() <= w[deltaE+8] ){
            //if ( ran2(&idum) <= w[deltaE+8] ) {
                spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
                accepted++;
            }
        }
    }
    accepted_vec[mc_count] = accepted;
    mc_count ++;
} // end of Metropolis sampling over spins

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
{
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran2(&idum)*(double)n_spins);
      int iy = (int) (ran2(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
	(spin_matrix[iy][periodic(ix,n_spins,-1)]+
	 spin_matrix[periodic(iy,n_spins,-1)][ix] +
	 spin_matrix[iy][periodic(ix,n_spins,1)] +
	 spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran2(&idum) <= w[deltaE+8] ) {
	spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double *total_average)
{
    double norm = 1/((double) (mcs));  // divided by total number of cycles
    double Etotal_average = total_average[0]*norm;
    double E2total_average = total_average[1]*norm;
    double Mtotal_average = total_average[2]*norm;
    double M2total_average = total_average[3]*norm;
    double Mabstotal_average = total_average[4]*norm;
    double HCapacity = (E2total_average - Etotal_average * Etotal_average)/(temperature * temperature);
    double Susceptibility = (M2total_average - Mabstotal_average*Mabstotal_average)/((double) temperature);
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
    double Mvariance = (M2total_average - Mtotal_average*Mtotal_average)/n_spins/n_spins;
    cout << "temperature " << temperature << endl;
    cout << "Etotal " << Etotal_average << endl;
    cout << "E2total_average " << E2total_average << endl;
    cout << "Mabstotal_average " << Mabstotal_average << endl;
    cout << "M2total_average " << M2total_average << endl;
    cout << "HCapacity " << HCapacity << endl;
    cout << "Susceptibility " << Susceptibility << endl;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
    ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins;
    ofile << setw(15) << setprecision(8) << Mvariance/temperature;
    ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins << endl;
} // end output function

/*
** The function
**         ran2()
** is a long periode (> 2 x 10^18) random number generator of
** L'Ecuyer and Bays-Durham shuffle and added safeguards.
** Call with idum a negative integer to initialize; thereafter,
** do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
    int            j;
    long           k;
    static long    idum2 = 123456789;
    static long    iy=0;
    static long    iv[NTAB];
    double         temp;

    if(*idum <= 0) {
    if(-(*idum) < 1) *idum = 1;
    else             *idum = -(*idum);
    idum2 = (*idum);
    for(j = NTAB + 7; j >= 0; j--) {
    k     = (*idum)/IQ1;
    *idum = IA1*(*idum - k*IQ1) - k*IR1;
    if(*idum < 0) *idum +=  IM1;
    if(j < NTAB)  iv[j]  = *idum;
    }
    iy=iv[0];
    }
    k     = (*idum)/IQ1;
    *idum = IA1*(*idum - k*IQ1) - k*IR1;
    if(*idum < 0) *idum += IM1;
    k     = idum2/IQ2;
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
    if(idum2 < 0) idum2 += IM2;
    j     = iy/NDIV;
    iy    = iv[j] - idum2;
    iv[j] = *idum;
    if(iy < 1) iy += IMM1;
    if((temp = AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran2()


/*
    * The function
    *      void  **matrix()
    * reserves dynamic memory for a two-dimensional matrix
    * using the C++ command new . No initialization of the elements.
    * Input data:
    *  int row      - number of  rows
    *  int col      - number of columns
    *  int num_bytes- number of bytes for each
    *                 element
    * Returns a void  **pointer to the reserved memory location.
*/

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


double analyticEavg(double Z, double beta){
    return (16 * exp(-8 * beta) - 16 * exp(8 * beta)) / Z;
}

double analyticMavg(double Z, double beta){
    return (8 * exp(8 * beta) + 16) / Z;
}

double analyticE2avg(double Z, double beta){
    return 256 * cosh(8 * beta) / Z;
}

double analyticM2avg(double Z, double beta){
    return (32 + 32 * exp(8 * beta)) / Z;
}

double analyticSusceptibility(double analytic_Mavg, double analytic_M2avg, double temperature){
    return (analytic_M2avg - analytic_Mavg * analytic_Mavg) / temperature;
}

double analyticHeatCapacity(double analytic_Eavg, double analytic_E2avg, double temperature){
    return (analytic_E2avg - analytic_Eavg * analytic_Eavg)/(temperature * temperature);
}




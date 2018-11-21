#include "time.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>

using namespace std;


random_device rd;
mt19937 engine(rd());
uniform_real_distribution <double> randd(0.0, 1.0);


ofstream ofile;
ofstream outfile;

//Inline function for periodic boundary conditions:
inline int periodic(int i, int limit, int add){
    return (i + limit + add) % (limit);
}


double analyticEavg(double, double);
double analyticE2avg(double, double);
double analyticMavg(double, double);
double analyticM2avg(double, double);
double analyticSusceptibility(double, double, double);
double analyticHeatCapacity(double, double, double);

// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&); //Bruker ikke denne

void initialize(int, int **, double&, double&, int);

void Metropolis(int, int&, int**, double&, double&, double *, vector <int>&, int&);

// prints to file the results of the calculations
void output(int, int, double, double *, double *, vector <double>, vector <double>, vector <double>, vector <double>, vector <double>, vector <double>);

//  Matrix memory allocation
//  allocate space for a matrix
void  **matrix(int, int, int);
//  free space for  a matrix
void free_matrix(void **);


int accepted_configs = 0;


//Random numbers:
double rrandom(){
    return randd(engine);
}

int my_rank;
int numprocs;

//main program:
int main(int argc, char* argv[]){

    //MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    double TimeStart, TimeEnd, TotalTime;
    TimeStart = clock();

    char *outfilename;
    double k = 1;
    double temp = 1.0;
    double initial_temp = 1.0;
    double final_temp = 2.4;
    double temp_step = 0.02;
    double beta = 1./(k*temp);  
    
    int n_spins = 60; 
    int accepted_configs = 0;
    int mcs = 100000; 
    int mc_counter = 0;
    //int my_rank, numprocs;
    double w[17], average[5], total_average[5], E, M; 
    double Z = 12 + 2*exp(8*beta) + 2*exp(-8*beta); 
    int multiplicity = pow(2, (n_spins*n_spins));

    vector <double> E_vec(mcs);
    vector <double> absM_vec(mcs);
    vector <double> E2_vec(mcs);
    vector <double> absM2_vec(mcs);
    vector <double> num_susceptibility_vec(mcs);
    vector <double> num_heatCapacity_vec(mcs);

    vector <int> accepted_configs_vec(mcs);

    int** spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));

    // Analytical Values
    double analytic_Eavg = analyticEavg(Z, beta); //function call
    double analytic_E2avg = analyticE2avg(Z, beta);
    double analytic_Mavg = analyticMavg(Z, beta);
    double analytic_M2avg = analyticM2avg(Z, beta);
    double analytic_Susceptibility = analyticSusceptibility(analytic_Mavg, analytic_M2avg, temp);
    double analytic_HeatCapacity = analyticHeatCapacity(analytic_Eavg, analytic_E2avg, temp);

    if (my_rank == 0){
        cout << "Antall Monte Carlo cycles: " << mcs << endl;
        cout << "Spins in each direction: " << n_spins << endl;
        cout << "Number of processors: " << numprocs << endl;
        //Skriver ut analytiske verdier:
        cout << "Analytical values:" << endl;
        cout << "Analytical expectation value of energy: " << analytic_Eavg/n_spins/n_spins << endl; //Eavg/n_spins/n_spins
        cout << "Analytical expectation value of energy^2: " << analytic_E2avg/n_spins/n_spins << endl; //Eavg/n_spins/n_spins
        cout << "Analytical expectation value of magnetization: " << analytic_Mavg/n_spins/n_spins << endl; //Mavg/n_spins/n_spins
        cout << "Analytical expectation value of magnetization^2: " << analytic_M2avg/n_spins/n_spins << endl; //Mavg/n_spins/n_spins
        cout << "Analytical susceptibility: " << analytic_Susceptibility << endl;
        cout << "Analytical heat capacity: " << analytic_HeatCapacity << endl;
    }


    //outfilename = "output.txt";
    //Open file for writing:
    if (my_rank == 0){
        ofile.open("output.txt"); //File for energy and magnetization
        outfile.open("configs.txt"); //File for number of accepted configs
    }

    for (double temp = initial_temp; temp <= final_temp; temp += temp_step){ 
        if (my_rank == 0){
            cout << endl;
            cout << "Temperatur= " << temp << " (values beneath):" << endl << endl;
        }
        E = M = 0; //Initialize energy and magnetization

        mc_counter = 0;

        //Set up array for possible energy changes:
        for (int de = -8; de <= 8; de++) {
            w[de+8] = 0;
        }

        for (int de = -8; de <= 8; de += 4) w[de+8] = exp(-de/temp);
        //Initialize array for expectation values:
        for (int i = 0; i < 5; i++) average[i] = total_average[i] = 0.; //Setter alle elementer lik null. [E, E^2, M, M^2, abs(M)]

        initialize(n_spins, spin_matrix, E, M, 1); //Kaller initialize
        //Start Monte Carlo computation:
        for (int cycles = 0; cycles < mcs; cycles++){
            Metropolis(n_spins, accepted_configs, spin_matrix, E, M, w, accepted_configs_vec, mc_counter);
            //Update expectation values:
            average[0] += E; average[1] += E*E;  //average is expectation values. I c) vil jeg plotte E jeg faar etter hver sweep gjennom lattice mot mcs
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
            //Fyller arrays for plotting:
            E_vec[cycles] = E;  //Inneholder E-verdiene jeg skal plotte. Antall E-verdier i E_vec vil vaere lik mcs (per temperatur).
            E2_vec[cycles] = E*E;
            absM_vec[cycles] = fabs(M);
            absM2_vec[cycles] = fabs(M) * fabs(M);

        }

        //num_heatCapacity_vec = ( E2_vec - E_vec*E_vec )/(temp*temp);
        //num_susceptibility_vec = ( absM2_vec - absM_vec*absM_vec )/temp;

        //Find total average
        for( int i =0; i < 5; i++){
            MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        //print results:
        if ( my_rank == 0) {
            output(n_spins, mcs, temp, average, total_average, E_vec, E2_vec, absM_vec, absM2_vec, num_susceptibility_vec, num_heatCapacity_vec); //temp er kalt temperature i Mortens program

        }

        //Print results:
        //output(n_spins, mcs, temp, average, E_vec, absM_vec); //temp er kalt temperature i Mortens program

        //free_matrix((void **) spin_matrix ); //free memory

        if (my_rank == 0){
            for (int i=0; i < mcs; i++){
                outfile << accepted_configs_vec[i] << endl; //Writing number of accepted configs to file "configs.txt"
            }
        }

    } //End Temp loop
    ofile.close(); //close output file (inneholder energier + magnetisering)
    outfile.close(); //contains number of accepted configs

    //End MPI
    MPI_Finalize();

    if (my_rank == 0){
        TotalTime = clock() - TimeStart;
        cout << "Runtime: " << TotalTime/1000000. << endl; //clock gives the time in microseconds
    }

    return 0;

}//End main

//All the needed functions:

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


//Calculates and prints numerical values:
void output(int n_spins, int mcs, double temp, double *average, double *total_average, vector <double> E_vec, vector <double> E2_vec, vector <double> absM_vec, vector <double> absM2_vec, vector <double> num_susceptibility_vec, vector <double> num_heatCapacity_vec)
{
  double norm = 1/((double) (mcs)*4);  // divided by total number of cycles
  double Eaverage = total_average[0]*norm;
  double E2average = total_average[1]*norm; //average of E^2
  double Maverage = total_average[2]*norm;
  double M2average = total_average[3]*norm;
  double Mabsaverage = total_average[4]*norm;
  double num_HeatCapacity = (E2average - Eaverage*Eaverage)/(temp*temp);
  double num_Susceptibility = (M2average - Mabsaverage*Mabsaverage)/((double) temp);
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average- Eaverage*Eaverage);
  double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  /*
  ofile << setw(15) << setprecision(8) << temp; //1 kolonne i filen output.txt
  ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins; //2 kolonne
  ofile << setw(15) << setprecision(8) << Evariance/temp/temp; //3 kolonne
  ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins; //4 kolonne
  ofile << setw(15) << setprecision(8) << Mvariance/temp; //5 kolonne
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl; //6 kolonne
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins << endl;
  */

  //Writing E- and M-values to file:
  //Plot i oppg c) (tror jeg):
  /*for (int i =1; i < mcs; i++){ //Skriver til fil temp, E, E2, absM, absM2, Cv, susceptibility
      ofile << temp << "  " << E_vec[i] << "  " << E2_vec[i]/(n_spins*n_spins) << "  " << absM_vec[i]<< "  " << absM2_vec[i]/(n_spins*n_spins) << "  " << num_heatCapacity_vec[i]/(n_spins*n_spins) << "  " << num_susceptibility_vec[i]/(n_spins*n_spins) << endl;  column: temp, 2 column: E, 3 column: abs(M), Cv, susceptibility
  //    //ofile << absM_vec[i] << endl;
  //}//End of writing E- and M-values to file */

  //Plot i oppg e):
  ofile << temp << "  " << Eaverage/(n_spins*n_spins) << "  " << E2average/(n_spins*n_spins) << "  " << Mabsaverage/(n_spins*n_spins) << "  " << M2average/(n_spins*n_spins) << " " << num_HeatCapacity/(n_spins*n_spins) << "  " << num_Susceptibility/(n_spins*n_spins) << endl;

  cout << "Numerical E average: " << Eaverage/n_spins/n_spins << endl;
  cout << "Numerical abs M average: " << Mabsaverage/n_spins/n_spins << endl;
  cout << "Numerical E^2 average: " << E2average/n_spins/n_spins << endl;
  cout << "Numerical M^2 average: " << M2average/n_spins/n_spins << endl;
  cout << "Numerical heat capacity: " << num_HeatCapacity/n_spins/n_spins << endl;
  cout << "Numerical susceptibility: " << num_Susceptibility/n_spins/n_spins << endl;
  cout << "Numerical E variance: " << Evariance << endl;
  cout << "Temp fra inni output: " << temp << endl;

} // end output function

//Function that performs the Metropolis algo:
void Metropolis(int n_spins, int& accepted_configs, int **spin_matrix, double & E, double & M, double *w, vector <int>& accepted_configs_vec, int& mc_counter){
    //Loop over all spins:
    for (int y = 0; y < n_spins; y++){
        for (int x = 0; x < n_spins; x++){  //For a 20x20 lattice we pick 400 random positions per mc cycle
            //find random position:
            int ix = floor(rrandom() * n_spins); //Random x and y which means a random spin is picked
            int iy = floor(rrandom() * n_spins);
            //Calculate energy difference: (likn 13.6 i forel notater) of total spin config compared to the last spin config(before the flip)

            int deltaE = 2*spin_matrix[iy][ix] *
            (spin_matrix[iy][periodic(ix, n_spins, -1)] +
            spin_matrix[periodic(iy, n_spins, -1)][ix] +
            spin_matrix[iy][periodic(ix, n_spins, 1)] +
            spin_matrix[periodic(iy, n_spins, 1)][ix] );



            //Performing the Metropolis test:
            if (rrandom() <= w[deltaE+8] ){
                spin_matrix[iy][ix] *= -1; //Flip one spin and accept new spin config
                //Update energy and magnetization:
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
                accepted_configs += 1; //Counts the accepted configs

            }
        }
    }
    accepted_configs_vec[mc_counter] = accepted_configs;
    mc_counter +=1;

    //return accepted_configs_vec;
}//End Metropolis function (Metropolis sampling over spins)

double analyticEavg(double Z, double beta){
    return (16*exp(-8*beta) - 16*exp(8*beta))/Z;
}

double analyticMavg(double Z, double beta){
    return (8*exp(8*beta) + 16)/Z;
}

double analyticE2avg(double Z, double beta){
    return 256*cosh(8*beta)/Z;
}

double analyticM2avg(double Z, double beta){
    return (32 + 32*exp(8*beta))/Z;
}

double analyticSusceptibility(double analytic_Mavg, double analytic_M2avg, double temp){
    return (analytic_M2avg - analytic_Mavg*analytic_Mavg)/temp;
}

double analyticHeatCapacity(double analytic_Eavg, double analytic_E2avg, double temp){
    return (analytic_E2avg - analytic_Eavg*analytic_Eavg)/(temp*temp);
}

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




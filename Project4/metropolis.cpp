//Function that performs the Metropolis algo:
void Metropolis(int n_spins, int& accepted_configs, int **spin_matrix, double & E, double & M, double *w, vector <int>& accepted_configs_vec, int& mc_counter){
    //Loop over all spins:
    for (int y = 0; y < n_spins; y++){
        for (int x = 0; x < n_spins; x++){ 
            int ix = floor(rrandom() * n_spins); 
            int iy = floor(rrandom() * n_spins);
            

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


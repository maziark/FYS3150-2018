int no_intervalls = mcs/numprocs;
int myloop_begin = my_rank*no_intervalls + 1;
int myloop_end = (my_rank+1)*no_intervalls;
if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;


for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
	// Find all the values!

}


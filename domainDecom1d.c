
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){

    int rank, size, i, provided;
    
    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2*3.1415; // Length of the domain
    

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells 
    int nxn_loc = nxc/size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);
    
    // define out function
    double *f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z

    for (i=1; i<(nxn_loc-1); i++)
      f[i] = sin(L_loc*rank + (i-1) * dx);
    
    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells
    int r_prev, r_next;
    if (rank == 0) {
        r_prev = size-1;
        r_next = 1;    
    }
    else if (rank == size-1) {
        r_prev = rank - 1;
        r_next = 0;
    }
    else {
        r_prev = rank - 1;
        r_next = rank + 1;
    }
    
    //printf("Rank = %d, r_prev = %d, r_next = %d \n", rank, r_prev, r_next);
    MPI_Send(&f[2], 1, MPI_DOUBLE, r_prev, 0, MPI_COMM_WORLD);
    MPI_Recv(&f[nxn_loc-1], 1, MPI_DOUBLE, r_next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&f[nxn_loc-3], 1, MPI_DOUBLE, r_next, 1, MPI_COMM_WORLD);
    MPI_Recv(&f[0], 1, MPI_DOUBLE, r_prev, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

    // here we finish the calculations

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i=1; i<(nxn_loc-1); i++)
      dfdx[i] = (f[i+1] - f[i-1])/(2*dx);

    
    // Print f values
    if (rank==3){ // print only rank 0 for convenience
        printf("My rank %d of %d\n", rank, size );
        printf("Here are my values for f including ghost cells\n");
        for (i=0; i<nxn_loc; i++)
	       printf("%f\n", f[i]);
        printf("Values for dfdx \n");
        for (i=1; i<nxn_loc-1; i++)
            printf("%f\n", dfdx[i]);
        printf("\n"); 
    }
        
    //Testing the correct values of the ghost cells
    double epsilon = 1e-2;
    if (fabs(f[0] - sin(L_loc*rank - dx)) > epsilon || fabs(f[nxn_loc-1] - sin(L_loc*(rank+1) + dx)) > epsilon) {
        printf("Incorrect ghost cell values, f[0] = %f, f[max] = %f, rank = %d \n correct 0 = %f, max = %f \n", f[0], f[nxn_loc-1], rank, sin(L_loc*rank - dx), sin(L_loc*(rank+1) + dx));
    }
    else printf("Correct ghost cell values rank %d \n", rank);
    

   
    //Testing derivatives on edges of the domain
    if (rank==0) {
        if (fabs(cos(0) - dfdx[1]) > epsilon) {
            printf("Left derivative incorrect %f != %f \n", cos(0), dfdx[1]);
        }
        else printf("Left derivative OK\n");
    }
    if (rank == size-1) {
        if (fabs(cos(L) - dfdx[nxn_loc-2]) > epsilon) {
            printf("Right derivative incorrect %f != %f \n", cos(L), dfdx[nxn_loc-2]);
        }
        else printf("Right derivative OK\n");
    }
    
    
    MPI_Finalize();
}







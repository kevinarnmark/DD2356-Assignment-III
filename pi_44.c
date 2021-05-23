
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    //int count = 0;
    double x, y, z, pi;
    int local_count = 0, flip = 1 << 24;
    int rank, num_ranks, i, iter, provided;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    
    double start_time, stop_time, elapsed_time;
    srand(SEED*rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    flip = flip/num_ranks;
    
    start_time = MPI_Wtime();
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < flip; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            local_count++;
        }
    }
    
    
    int global_count = 0;
    MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        // Estimate Pi and display the result
        pi = ((double)global_count / (double)(flip * num_ranks)) * 4.0;
    }
    
    stop_time = MPI_Wtime();
    elapsed_time = stop_time - start_time;
    
    if (rank == 0) {
        printf("pi = %f\n", pi);
        printf("Number of processes: %d\n", num_ranks);
        printf("Execution Time = %f\n", elapsed_time);
    }
    MPI_Finalize();
    //printf("The result is %f\n", pi);
    
    return 0;
}


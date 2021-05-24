
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
    
    int tree_depth = log2(num_ranks);
    int local_recv = 0;
    if (rank % 2 == 1) {
        MPI_Send(&local_count, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD); 
    }
    for (int i = 0; i < tree_depth; i++) {
        if (rank % (1 << i+1) == 0) {
            MPI_Recv(&local_recv, 1, MPI_INT, rank + pow(2, i), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_count += local_recv;
        }
        else if (rank % (1 << i) == 0 && rank % 2 != 1) {
            MPI_Send(&local_count, 1, MPI_INT, rank - pow(2, i), 0, MPI_COMM_WORLD); 
        }
    }
    if (rank == 0) {
        // Estimate Pi and display the result
        pi = ((double)local_count / (double)(flip * num_ranks)) * 4.0;
    
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


/*
Use MPI and OpenMP to approximate pi.
compile:
    mpicc hybrid_pi.c -o hybrid_pi -fopenmp -lm
execute:
    srun -N2 -n6 -c4 ./hybrid_pi 100000000
    
    -N2 means 2 nodes
    -n6 means 6 processes
    -c4 means 4 CPUs per processes
    You can use sbatch as well.
    Try different number of threads!
*/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv) {
    int rank, size;
    int numParts;
     
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    numParts = atoll(argv[1]);
    double myArea = 0.0;
    int i; 
    #pragma omp parallel for schedule(guided, 10)
    for(i = rank; i < numParts; i+=size){
        #pragma omp critical
        myArea += (sqrt(1 - pow((double)i/(double)numParts, 2))) / (double)numParts;
    }
    double area;
    MPI_Reduce(&myArea, &area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0){
        printf("%.10f\n", 4 * area);       	
    }

    MPI_Finalize();

    return 0;
}

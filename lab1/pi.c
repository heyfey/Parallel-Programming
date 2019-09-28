#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv) {
    int rank, size;
    double numParts;
     
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    numParts = atoll(argv[1]);
    // printf("rank = %d numParts = %f\n", rank, numParts);

    double myArea = 0.0;
    double i; 

    for(i = rank; i < numParts; i+=size){
        myArea += (sqrt(1 - pow(i/numParts, 2))) / numParts;
    }
    // printf("rank = %d my area = %.10f\n", rank, myArea);

    double area;
    MPI_Reduce(&myArea, &area, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    if(rank == 0){
        // printf("rank = %d total area = %.10f\n", rank, 4 * area);
        printf("%.10f\n", 4 * area);       	
    }
  
    // printf("size= %d  my rank= %d\n", size, rank);

    MPI_Finalize();

    return 0;
}

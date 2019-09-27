#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv) {
  int rank, size;
  long long numParts;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  numParts = atoll(argv[1]);

  MPI_Finalize();
  
  return 0;
}

#include <cstdio>
#include <mpi.h>
int main(int argc, char** argv) {
    //printf("%s\n", argv[1]);
	MPI_Init(&argc,&argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_File f;
    int rc;
    // filename: cases/01.in
	rc = MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
    if(rc){
        printf("faied to read file.\n");
    }
	float data[1];
	MPI_File_read_at(f, sizeof(float) * rank, &data, 1, MPI_FLOAT, MPI_STATUS_IGNORE);
	printf("rank %d got float: %f\n", rank, data[0]);
	// MPI_File_write_at(f, ...
	MPI_Finalize();
}

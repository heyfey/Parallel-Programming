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

	float data[2];
	MPI_File_read_at(f, sizeof(float) * rank * 2, &data, 2, MPI_FLOAT, MPI_STATUS_IGNORE); //JUST TESTING.
	printf("rank %d got float: %f %f\n", rank, data[0], data[1]);
	
    MPI_File out;
    rc = MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out);
    if(rc){
        printf("open output file failed.\n");
    }
    MPI_File_write_at(out, sizeof(float) * rank * 2, &data, 2, MPI_FLOAT, MPI_STATUS_IGNORE); // JUST TESTING.

	MPI_Finalize();
}

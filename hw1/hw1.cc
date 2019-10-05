#include <cstdio>
#include <cstdlib>
#include <mpi.h>
int main(int argc, char** argv) {
    
	MPI_Init(&argc,&argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int N = atoll(argv[1]);
    int buf_size = N / size;
    int buf_size_last;
    if (N % size == 0){
        buf_size_last = N / size;
    }else{
        buf_size_last = N % size;
    }

    //Read file.
    MPI_File f_in;
    int rc;
	rc = MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &f_in);
    if(rc){
        printf("faied to read file.\n");
    }

	float * data = new float[buf_size];
	MPI_File_read_at(f_in, sizeof(float) * rank * buf_size, data, buf_size, MPI_FLOAT, MPI_STATUS_IGNORE);
	
    // if only 1 process?
    // if processes > N ?

    //Main work here.
    int phase = 0;
    while(1){
        if(phase % 2 == 0){ // even phase.

        }else{ // odd phase.


        }
        //reduce & break....
        break;
        phase++;
    }

	//Write file.
    MPI_File f_out;
    rc = MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_out);
    if(rc){
        printf("open output file failed.\n");
    }
    MPI_File_write_at(f_out, sizeof(float) * rank * buf_size, data, buf_size, MPI_FLOAT, MPI_STATUS_IGNORE);

	MPI_Finalize();
}

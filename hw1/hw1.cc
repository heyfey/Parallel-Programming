#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <mpi.h>

int find_partner(int phase, int rank, int size);
void print_array(float *arr, int len);
void swap_array(float **a, float **b);
int merge_low(float *my_buf, float *recv_buf, float *tmp_buf, int size);
int merge_high(float *my_buf, float *recv_buf, float *tmp_buf, int size);
int merge(int phase, int rank, float *my_buf, float *recv_buf, float *tmp_buf, int size);

int main(int argc, char** argv) {
    
    MPI_Init(&argc,&argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;  
    MPI_Request req;
    int N = atoll(argv[1]); 
    
    MPI_File f_in, f_out;
    MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &f_in);
    MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_out);

    double begin, end;
    double IO_time, computing_time, communication_time = 0;
    
    // If only oen process, or small input. 
    if(size == 1 || size >= N){
        if(rank == 0){
            float *buf;
            buf = new float[N];
            begin = MPI_Wtime();
            MPI_File_read(f_in, buf, N, MPI_FLOAT, MPI_STATUS_IGNORE);
            end = MPI_Wtime();
            IO_time += (end - begin);

            begin = MPI_Wtime();
            std::sort(buf, buf + N);
            end = MPI_Wtime();
            computing_time += (end - begin);

            begin = MPI_Wtime();
            MPI_File_write(f_out, buf, N, MPI_FLOAT, MPI_STATUS_IGNORE);
            end = MPI_Wtime();
            IO_time += (end - begin);

            printf("[%2d] IO time: %f, Computing time: %f, Communication time: %f, Total time: %f\n", rank, IO_time, computing_time, communication_time, IO_time + computing_time + communication_time);
        }
        MPI_File_close(&f_in);
        MPI_File_close(&f_out);
        MPI_Finalize();
        return 0;
    }

    // Determine buffer size.
    int buf_size = N / size;
    int buf_size_last;
    if (N % size == 0){
        buf_size_last = N / size;
    }else{
        buf_size = buf_size + 1;
        buf_size_last = N % (buf_size);
    }
    // printf("rank = %d, size = %d, N = %d, buf_size = %d, buf_size_last = %d\n", rank, size, N, buf_size, buf_size_last);

    float *my_buf, *recv_buf, *tmp_buf;
    my_buf = new float[buf_size];
    recv_buf = new float[buf_size];
    tmp_buf = new float[buf_size];
    begin = MPI_Wtime();
    MPI_File_read_at(f_in, sizeof(float) * rank * buf_size, my_buf, buf_size, MPI_FLOAT, MPI_STATUS_IGNORE);	
    end = MPI_Wtime();
    IO_time += (end - begin);

    // Deal with the last process: fill the buffer to buf_size with INF.
    if(rank == size - 1 && buf_size != buf_size_last){
        for(int i = buf_size_last; i < buf_size; i++){
            my_buf[i] = INFINITY;
        }
    }
    
    // Local sort.
    begin = MPI_Wtime();
    std::sort(my_buf, my_buf + buf_size);
    end = MPI_Wtime();
    computing_time += (end - begin);

    // Main work here.
    int phase, done, alldone = 0;
    int partner;
    while(!alldone){
        partner = find_partner(phase, rank, size);
        
        begin = MPI_Wtime();
        MPI_Isend(my_buf, buf_size, MPI_FLOAT, partner, 0, MPI_COMM_WORLD, &req);
        MPI_Irecv(recv_buf, buf_size, MPI_FLOAT, partner, MPI_ANY_TAG, MPI_COMM_WORLD, &req); 
        MPI_Wait(&req, &status);    
        end = MPI_Wtime();
        communication_time += (end - begin);
        
        
        if(partner != MPI_PROC_NULL){
            begin = MPI_Wtime();
            done = merge(phase, rank, my_buf, recv_buf, tmp_buf, buf_size);
            end = MPI_Wtime();
            computing_time += (end - begin);

            swap_array(&my_buf, &tmp_buf);
            if(done == 2) swap_array(&my_buf, &tmp_buf);
            if(done == 3){ swap_array(&my_buf, &tmp_buf); swap_array(&my_buf, &recv_buf); done = 0;}
        }
        
        begin = MPI_Wtime();
        MPI_Allreduce(&done, &alldone, 1, MPI_FLOAT, MPI_LAND, MPI_COMM_WORLD);
        end = MPI_Wtime();
        communication_time += (end - begin);

        phase++;
    }
    
    begin = MPI_Wtime();  
    if(rank == size - 1 && buf_size != buf_size_last){
        MPI_File_write_at(f_out, sizeof(float) * rank * buf_size, my_buf, buf_size_last, MPI_FLOAT, MPI_STATUS_IGNORE);
    }else{
        MPI_File_write_at(f_out, sizeof(float) * rank * buf_size, my_buf, buf_size, MPI_FLOAT, MPI_STATUS_IGNORE);
    }
    end = MPI_Wtime();
    IO_time += (end - begin);

    printf("[%2d] IO time: %f, Computing time: %f, Communication time: %f, Total time: %f\n", rank, IO_time, computing_time, communication_time, IO_time + computing_time + communication_time);
    MPI_File_close(&f_in);
    MPI_File_close(&f_out);
    MPI_Finalize();
  
}


int find_partner(int phase, int rank, int size){
    int partner;
    if(phase % 2 == 0){ // even phase
        if(rank % 2 == 0){
            partner = rank + 1;
        }else{
            partner = rank - 1;
        }
    }else{ // odd phase
        if(rank % 2 == 0){
            partner = rank - 1;
        }else{
            partner = rank + 1;
        }
    }
    if(partner == -1 || partner == size){
        partner = MPI_PROC_NULL;
    }
    return partner;
}

int merge(int phase, int rank, float *my_buf, float *recv_buf, float *tmp_buf, int size){
    int done = 0;
    if(phase % 2 == 0){ // even phase
        if(rank % 2 == 0){
            done = merge_low(my_buf, recv_buf, tmp_buf, size);
        }else{
            done = merge_high(my_buf, recv_buf, tmp_buf, size);
        }
    }else{ // odd phase
        if(rank % 2 == 0){
            done = merge_high(my_buf, recv_buf, tmp_buf, size);
        }else{
            done = merge_low(my_buf, recv_buf, tmp_buf, size);
        }
    }
    return done;
}

int merge_low(float *my_buf, float *recv_buf, float *tmp_buf, int size){
    if(my_buf[size-1] <= recv_buf[0]) return 2;
    if(recv_buf[size-1] <= my_buf[0]) return 3;
    int done = 1;
    int i = 0, j = 0;
    for(int idx = 0; idx < size; idx++){
        if(my_buf[i] <= recv_buf[j]){
            tmp_buf[idx] = my_buf[i];
            i++;
        }else{
            tmp_buf[idx] = recv_buf[j];
            j++;
            done = 0;
        }
    }
    return done;
}

int merge_high(float *my_buf, float *recv_buf, float *tmp_buf, int size){
    if(my_buf[0] >= recv_buf[size-1]) return 2;
    if(recv_buf[0] >= my_buf[size-1]) return 3;
    int done = 1;
    int i = size - 1;
    int j = size - 1;
    for(int idx = size-1; idx >= 0; idx--){
        if(my_buf[i] >= recv_buf[j]){
            tmp_buf[idx] = my_buf[i];
            i--;
        }else{
            tmp_buf[idx] = recv_buf[j];
            j--;
            done = 0;
        }
    }
    return done;
}

void swap_array(float **a, float **b){
    float *tmp;
    tmp = *a;
    *a = *b;
    *b = tmp;
}

void print_array(float *arr, int len){
    for(int i = 0; i < len; i++){
        printf(" %f", arr[i]);
    }
    printf("\n");
}



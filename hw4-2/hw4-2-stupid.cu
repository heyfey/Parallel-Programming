#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <chrono>
#include <iostream>
#include <omp.h>

#define BLOCK_FACTOR 32
#define NUM_THREAD 1024
#define PART (BLOCK_FACTOR * BLOCK_FACTOR) / NUM_THREAD

const int INF = ((1 << 30) - 1);
void input(char* inFileName);
void output(char* outFileName);
void printArray(int* Dist);

void block_FW();
int ceil(int a, int b);

const int V = 40010;
int n, m;

int numDevs = 0;
int* Dist;
size_t pitch;

clock_t begin, end;
double IO_time = 0;
double kernel_time = 0;

int main(int argc, char* argv[]) {
    cudaGetDeviceCount(&numDevs);
    printf("devices count: %d\n", numDevs);

    cudaMallocManaged(&Dist, V * V * sizeof(int));

    begin = clock();
    input(argv[1]);
    end = clock();
    IO_time += (double) (end - begin) / CLOCKS_PER_SEC;
    //printArray(Dist);

    begin = clock();
    block_FW();
    end = clock();
    kernel_time += (double) (end - begin) / CLOCKS_PER_SEC;
    //printArray(Dist);

    begin = clock();
    output(argv[2]);
    end = clock();
    IO_time += (double) (end - begin) / CLOCKS_PER_SEC;

    cudaFree(&Dist);

    printf("I/O time: %f secs.\n", IO_time);
    printf("GPU kernel time: %f secs.\n", kernel_time);
    return 0;
}

void input(char* infile) {
    FILE* file = fopen(infile, "rb");
    fread(&n, sizeof(int), 1, file);
    fread(&m, sizeof(int), 1, file);
    printf("V = %d, E = %d\n", n, m);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                Dist[i * n + j] = 0;
            } else {
                Dist[i * n + j] = INF;
            }
        }
    }
    int pair[3];
    for (int i = 0; i < m; ++i) {
        fread(pair, sizeof(int), 3, file);
        Dist[pair[0] * n + pair[1]] = pair[2];
    }
    fclose(file);
}

void output(char* outFileName) {
    FILE* outfile = fopen(outFileName, "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (Dist[i * n + j] >= INF)
                Dist[i * n + j] = INF;
            fwrite(&Dist[i * n + j], sizeof(int), 1, outfile);
        }
    }
    fclose(outfile);
}

void printArray(int* Dist){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            printf("%6d ", Dist[i * n + j]);
        }
        printf("\n");
    } 
}

int ceil(int a, int b) { return (a + b - 1) / b; }

/*
__device__ void block_cal(int* C, int* A, int* B, int j, int i) {
    for (int k = 0; k < BLOCK_FACTOR; k++) {
        int sum = A[i*BLOCK_FACTOR + k] + B[k*BLOCK_FACTOR + j];
        if (C[i*BLOCK_FACTOR + j] > sum) {
            C[i*BLOCK_FACTOR + j] = sum;
        }
        __syncthreads();
    }
}
*/

__global__ void cal_phase1(int Round, int block_start_x, int block_start_y, int* Dist, int n) {
    int x = threadIdx.y;
    int y = threadIdx.x;
    int i = (block_start_x + blockIdx.x) * BLOCK_FACTOR + x;
    int j = (block_start_y + blockIdx.y) * BLOCK_FACTOR + y;

    __shared__ int C[BLOCK_FACTOR][BLOCK_FACTOR+1];
    int i_start = i * n;
    int rb = Round * BLOCK_FACTOR;
    C[x][y] = Dist[i_start + j];
    __syncthreads();

    if(i >= n || j >= n) return;

    for (int k = 0; k < BLOCK_FACTOR && rb + k < n; k++) {
        int sum = C[x][k] + C[k][y];
        if (sum < C[x][y])
            C[x][y] = sum;
        __syncthreads();
    }
    Dist[i_start + j] = C[x][y];
}

__global__ void cal_phase2A(int Round, int block_start_x, int block_start_y, int* Dist, int n) {
    int x = threadIdx.y;
    int y = threadIdx.x;
    int i = (block_start_x + blockIdx.x) * BLOCK_FACTOR + x;
    int j = (block_start_y + blockIdx.y) * BLOCK_FACTOR + y;

    __shared__ int C[BLOCK_FACTOR][BLOCK_FACTOR+1];
    __shared__ int A[BLOCK_FACTOR][BLOCK_FACTOR+1];
    int i_start = i * n;
    int rb = Round * BLOCK_FACTOR;
    C[x][y] = Dist[i_start + j];
    A[x][y] = Dist[i_start + rb + y];
    __syncthreads();

    if(i >= n || j >= n) return;
    
    for (int k = 0; k < BLOCK_FACTOR && rb + k < n; k++) {
        int sum = A[x][k] + C[k][y];
        if (sum < C[x][y])
            C[x][y] = sum;
        __syncthreads();
    }
    Dist[i_start + j] = C[x][y];
}


__global__ void cal_phase2B(int Round, int block_start_x, int block_start_y, int* Dist, int n) {
    int x = threadIdx.y;
    int y = threadIdx.x;
    int i = (block_start_x + blockIdx.x) * BLOCK_FACTOR + x;
    int j = (block_start_y + blockIdx.y) * BLOCK_FACTOR + y;

    __shared__ int C[BLOCK_FACTOR][BLOCK_FACTOR+1];
    __shared__ int B[BLOCK_FACTOR][BLOCK_FACTOR+1];
    int i_start = i * n;
    int rb = Round * BLOCK_FACTOR;
    C[x][y] = Dist[i_start + j];
    B[x][y] = Dist[((rb + x) * n) + j];
    __syncthreads();

    if(i >= n || j >= n) return;

    for (int k = 0; k < BLOCK_FACTOR && rb + k < n; k++) {
        int sum = C[x][k] + B[k][y];
        if (sum < C[x][y])
            C[x][y] = sum;
        __syncthreads();
    }
    Dist[i_start + j] = C[x][y];
}

__global__ void cal_phase3(int Round, int block_start_x, int block_start_y, int* Dist, int n) {
    int nn = n;    
    int x = threadIdx.y;
    int y = threadIdx.x;
    int i = (block_start_x + blockIdx.x) * BLOCK_FACTOR + x;
    int j = (block_start_y + blockIdx.y) * BLOCK_FACTOR + y;

    __shared__ int A[BLOCK_FACTOR][BLOCK_FACTOR+1];
    __shared__ int B[BLOCK_FACTOR][BLOCK_FACTOR+1];
    int i_start = i * nn;
    int rb = Round * BLOCK_FACTOR;
    int c = Dist[i_start + j];
    A[x][y] = Dist[i_start + rb + y];
    B[x][y] = Dist[((rb + x) * nn) + j];
    __syncthreads();

    if(i >= nn || j >= nn) return;

    for (int k = 0; k < BLOCK_FACTOR && rb + k < nn; k++) {
        int sum = A[x][k] + B[k][y];
        if (sum < c) 
            c = sum;
    }
    Dist[i_start + j] = c;
}

void block_FW() {
    int round = ceil(n, BLOCK_FACTOR);
    dim3 block_dim(BLOCK_FACTOR, BLOCK_FACTOR);
    dim3 grid;
    cudaThreadSynchronize();
    for (int r = 0; r < round; ++r) {
        //printf("%d %d\n", r, round);
        //fflush(stdout);
        int bs = round - r - 1;
        cudaThreadSynchronize();
        
        /* phase 1 */
        cudaSetDevice(0);
        grid.x = 1; 
        grid.y = 1;
        cal_phase1 <<<grid, block_dim >>> (r, r, r, Dist, n);
        cudaThreadSynchronize();
        /* phase 2 */
        if (r > 0) {
            grid.x = 1;
            grid.y = r;
            cal_phase2A <<<grid, block_dim>>> (r, r, 0, Dist, n);
        }
        if (bs > 0) {
            grid.x = 1;
            grid.y = bs;
            cal_phase2A <<<grid, block_dim>>> (r, r, r + 1, Dist, n);
        }
        if (r > 0) {
            grid.x = r;
            grid.y = 1;
            cal_phase2B <<<grid, block_dim>>> (r, 0, r, Dist, n);
        }
        if (bs > 0) {
            grid.x = bs;
            grid.y = 1;
            cal_phase2B <<<grid, block_dim>>> (r, r + 1, r, Dist, n);
        }
        cudaThreadSynchronize();
        /* phase 3 */
#pragma omp parallel num_threads(2)
{
        int tid = omp_get_thread_num();
        cudaSetDevice(tid);
        if (tid == 0) {
        if (r > 0) {
            grid.x = r;
            grid.y = r;
            cal_phase3 <<<grid, block_dim>>> (r, 0, 0, Dist, n);
            cudaDeviceSynchronize();
        }
        }
        if (tid == 1) {
        if (r > 0 && bs > 0) {
            grid.x = r;
            grid.y = bs;
            cal_phase3 <<<grid, block_dim>>> (r, 0, r + 1, Dist, n);
            cudaDeviceSynchronize();
            
            grid.x = bs;
            grid.y = r;
            cal_phase3 <<<grid, block_dim>>> (r, r + 1, 0, Dist, n);
            cudaDeviceSynchronize();
        }
        }
        if (tid == 0) {
        if (bs > 0) {
            grid.x = bs;
            grid.y = bs;
            cal_phase3 <<<grid, block_dim>>> (r, r + 1, r + 1, Dist, n);
            cudaDeviceSynchronize();
        }
        }
}
    }
    cudaDeviceSynchronize();
}























# Parallel programming

course at NTHU CS 2019 fall

## 1. [odd-even sort](https://en.wikipedia.org/wiki/Odd%E2%80%93even_sort)

+ **MPI**
+ optimized algorithm to minimize message size
+ asynchronous communication (non-blocking send/recv)

## 2. [Mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set)

`hw2/hw2_hybrid_dynamic_p_v.c`

+ **MPI + pthreads + OpenMP**
+ leader/follower architecture
+ load balance with dynamic scheduling
+ overlapped computing and file writing
+ vetorization with Intel SSE3 (SIMD)


## 3. [all-pairs shortest path](https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm) (cpu)

+ **OpenMP**
+ implemented blocked-Floyd-Warshall algorithm to utilize cache locality

## 4. all-pairs shortest path (CUDA)

### 4-1: single-GPU
+ utilized NVIDIA Pascal GPU memory hierarchy : shared memory, registers
+ fine-tuned block size and kernel size
+ resolved bank conflicts
### 4-2: multi-GPU
+ minimized peer-to-peer communication

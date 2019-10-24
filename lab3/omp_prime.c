/*
Use OpenMP to count prime numbers <= N.
compile:
    gcc -lm omp_prime.c -o omp_prime -fopenmp
execute:
    srun -c4 -n1 ./omp_prime 100000000

    -n1 means 1 process
    -c4 means 4 CPUs per processes
    100000000 = N
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
int isPrime(int i) {
    int sqrt_i = (int)sqrt((double)i);
    int j;
    for (j = 2; j <= sqrt_i; ++ j) {
        if (i % j == 0) return 0;
    }
    return 1;
}
int main(int argc, char** argv) {
    assert(argc == 2);
    int N = atoi(argv[1]);
    int count = 0;
    int i;
#pragma omp parallel for schedule(guided, 100)
    for (i = 2; i<= N; i++) {
        if (isPrime(i)) {
    #pragma omp critical
        count++;
        }
    }
    
    printf("There are %d prime numbers <= %d\n", count, N);
    return 0;
}

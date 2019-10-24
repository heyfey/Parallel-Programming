/*
Use pthread to approximate pi.
compile:
    gcc pi_pthread.c -o pi_pthread -pthread -lm
execute:
    srun -c4 -n1 ./pi_pthread 4 500000000

    number of threads = 4
    number of points = 500000000
*/
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double pi = 0.0;
pthread_mutex_t mutex;

struct thread_data {
    double num_threads;
    double num_parts;
    double t;
};

void* hello(void* threadid) {
    struct thread_data* tid = (struct thread_data*)threadid;
    double area = 0.0;
    double i;
    for (i = tid->t; i < tid->num_parts; i += tid->num_threads){
        area += (sqrt(1 - pow(i/tid->num_parts, 2))) / tid->num_parts;
    }
    pthread_mutex_lock(&mutex);
    pi += area;
    pthread_mutex_unlock(&mutex);
    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
    assert(argc == 3);
    int num_threads = atoi(argv[1]);
    int num_parts = atoi(argv[2]);
    pthread_t threads[num_threads];
    pthread_mutex_init(&mutex, NULL);
    int rc;
    struct thread_data ID[num_threads];
    int t;
    for (t = 0; t < num_threads; t++) {
        ID[t].t = t;
        ID[t].num_threads = num_threads;
        ID[t].num_parts = num_parts;
        rc = pthread_create(&threads[t], NULL, hello, (void*)&ID[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
    for (t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL);
    }
    printf("%.10f\n", 4 * pi);    
    pthread_mutex_destroy(&mutex);
    pthread_exit(NULL);
}

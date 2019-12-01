#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_V 6000 
#define MAX_DIST 1073741823

int graph[MAX_V][MAX_V];
//int d[MAX_V][MAX_V];
int total_v, total_e;

void init_graph() {
    for (int i = 0; i < total_v; i++) {
        for (int j = 0; j < total_v; j++) {
            graph[i][j] = MAX_DIST;
//            d[i][j] = MAX_DIST;
        }
    }
    for (int i = 0; i < total_v; i++) {
        graph[i][i] = 0;
//        d[i][i] = 0;
    }
}

void print_graph() {
    for (int i = 0; i < total_v; i++) {
        for (int j = 0; j < total_v; j++) {
            printf("%11d ", graph[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void output_graph(FILE *fp) {
    for (int i = 0; i < total_v; i++) {
        for (int j = 0; j < total_v; j++) {
            fwrite(&graph[i][j], sizeof(int), 1, fp);
            //fwrite(&d[i][j], sizeof(int), 1, fp); //dijskatra
        }
    }
}

void APSP() {
    //volatile int done = 1;
    //int done = 1;
    for (int k = 0; k < total_v; k++){
        //if (done && k != 0) continue;
        //done = 1;
#pragma omp parallel for schedule(static) //shared(done)
        for (int i = 0; i < total_v; i++) {
            for (int j = 0; j < total_v; j++) {
                if (graph[i][k] + graph[k][j] < graph[i][j])
                    graph[i][j] = graph[i][k] + graph[k][j];
                    //done = 0;
            }
        }
        //if (done) break;
        //done = 1;
    }
}
/*
void dijkstra(int v) {
    int visit[total_v];
    for (int i = 0; i < total_v; i++) visit[i] = 0;
    for (int k = 0; k < total_v; k++) {
        int a = -1, b = -1, min = MAX_DIST;
        for (int i = 0; i < total_v; i++) {
            if (!visit[i] && d[v][i] < min){
                a = i;
      volatile int done = 1;          min = d[v][i];
            }
        }
        if (a == -1) break;
        visit[a] = 1;
        for (int b = 0; b < total_v; b++) {
            if (!visit[b] && d[v][a] + graph[a][b] < d[v][b]) {
                d[v][b] = d[v][a] + graph[a][b];
            }
        }
    }
}
*/
int main(int argc, char** argv) {
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    printf("%d cpus available\n", CPU_COUNT(&cpu_set));

    assert(argc == 3);
    FILE *f_in;
    FILE *f_out;
    f_in = fopen(argv[1], "rb");
    f_out = fopen(argv[2], "wb+");

    fread(&total_v, sizeof(int), 1, f_in);
    fread(&total_e, sizeof(int), 1, f_in);
    printf("total v: %d, total e: %d\n", total_v, total_e);
    fflush(stdout);

    init_graph();
    for (int i = 0; i < total_e; i++) {
        int src, dst, w;
        fread(&src, sizeof(int), 1, f_in);
        fread(&dst, sizeof(int), 1, f_in);
        fread(&w, sizeof(int), 1, f_in);
        graph[src][dst] = w;
        //printf("%d %d %d\n", src, dst, w);
    }

    //print_graph();
    APSP();
    /*
#pragma omp parallel for schedule(static)
    for (int i = 0; i < total_v; i++) {
        dijkstra(i);
    }
    */
    //print_graph();
    output_graph(f_out);
    
    fclose(f_in);
    fclose(f_out);
}

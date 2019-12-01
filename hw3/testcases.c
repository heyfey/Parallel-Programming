#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_V 6000
#define MIN_DIST 0
#define MAX_DIST 1000

void main () {
    FILE *fp;
    fp = fopen("testcase.in", "wb+");
    srand(time(NULL));

    int v = MAX_V;
    int e = v * (v - 1);
    fwrite(&v, sizeof(int), 1, fp);
    fwrite(&e, sizeof(int), 1, fp);
    
    for (int i = 0; i < v; i++) {
        for (int j = 0; j < v; j++) {
            if (i != j) {
                fwrite(&i, sizeof(int), 1, fp);
                fwrite(&j, sizeof(int), 1, fp);
                int r = rand() % 1000;
                fwrite(&r, sizeof(int), 1, fp);
            }
        }
    }



}

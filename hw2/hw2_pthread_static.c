#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void write_png(const char* filename, int iters, int width, int height, const int* buffer) {
    FILE* fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
    png_write_info(png_ptr, info_ptr);
    png_set_compression_level(png_ptr, 1);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; ++x) {
            int p = buffer[(height - 1 - y) * width + x];
            png_bytep color = row + x * 3;
            if (p != iters) {
                if (p & 16) {
                    color[0] = 240;
                    color[1] = color[2] = p % 16 * 16;
                } else {
                    color[0] = p % 16 * 16;
                }
            }
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

typedef struct thread_data {
    int num_threads;
    int t;
    int iters;
    double left;
    double right;
    double lower;
    double upper;
    int width;
    int height;
    int* image;
} thread_data;

void* mandelbrot_set_t(void* thread_d) {
    thread_data* td = (thread_data*)thread_d;
    /* mandelbrot set */
    for (int j = td->t; j < td->height; j += td->num_threads) {
        double y0 = j * ((td->upper - td->lower) / td->height) + td->lower;
        for (int i = 0; i < td->width; ++i) {
            double x0 = i * ((td->right - td->left) / td->width) + td->left;

            int repeats = 0;
            double x = 0;
            double y = 0;
            double length_squared = 0;
            while (repeats < td->iters && length_squared < 4) {
                double temp = x * x - y * y + x0;
                y = 2 * x * y + y0;
                x = temp;
                length_squared = x * x + y * y;
                ++repeats;
            }
            td->image[j * td->width + i] = repeats;
        }
    }
}

int main(int argc, char** argv) {
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    printf("%d cpus available\n", CPU_COUNT(&cpu_set));

    /* argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    int iters = strtol(argv[2], 0, 10);
    double left = strtod(argv[3], 0);
    double right = strtod(argv[4], 0);
    double lower = strtod(argv[5], 0);
    double upper = strtod(argv[6], 0);
    int width = strtol(argv[7], 0, 10);
    int height = strtol(argv[8], 0, 10);

    /* allocate memory for image */
    int* image = (int*)malloc(width * height * sizeof(int));
    assert(image);

    pthread_t threads[CPU_COUNT(&cpu_set)];
    thread_data td[CPU_COUNT(&cpu_set)];
    int rc;
    int t;
    for (t = 0; t < CPU_COUNT(&cpu_set); t++) {
        td[t].t = t;
        td[t].num_threads = CPU_COUNT(&cpu_set);
        td[t].iters = iters; 
        td[t].left = left;
        td[t].right = right;
        td[t].lower = lower;
        td[t].upper = upper;
        td[t].width = width;
        td[t].height = height;
        td[t].image = image;
        rc = pthread_create(&threads[t], NULL, mandelbrot_set_t, (void*)&td[t]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
    for (t = 0; t < CPU_COUNT(&cpu_set); t++) {
        pthread_join(threads[t], NULL);
    }

    /* draw and cleanup */
    write_png(filename, iters, width, height, image);
    free(image);
}

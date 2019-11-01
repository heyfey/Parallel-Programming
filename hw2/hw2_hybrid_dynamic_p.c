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
#include <mpi.h>

typedef struct write_png_data {
    const char* filename;
    int iters;
    int width;
    int height;
    const int* buffer;
    volatile int* ready_to_write;
} write_png_data;

void* write_png_t(void* thread_d) {
    write_png_data* td = (write_png_data*)thread_d;

    FILE* fp = fopen(td->filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, td->width, td->height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
    png_write_info(png_ptr, info_ptr);
    png_set_compression_level(png_ptr, 1);
    size_t row_size = 3 * td->width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);

    for (int y = td->height-1; y >= 0; y--) {
        while(!td->ready_to_write[y]);
        memset(row, 0, row_size);
        for (int x = 0; x < td->width; ++x) {
            int p = td->buffer[y * td->width + x];
            png_bytep color = row + x * 3;
            if (p != td->iters) {
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

void mandelbrot_set(int rank, int size, int iters, double left, double right, double lower, double upper, int width, int height, int* image_row, int num_threads, int j) {
  double y0 = j * ((upper - lower) / height) + lower;
#pragma omp parallel for schedule(dynamic, 1) num_threads(num_threads)
  for (int i = 0; i < width; ++i) {
      double x0 = i * ((right - left) / width) + left;
  
      int repeats = 0;
      double x = 0;
      double y = 0;
      double length_squared = 0;
      while (repeats < iters && length_squared < 4) {
          double temp = x * x - y * y + x0;
          y = 2 * x * y + y0;
          x = temp;
          length_squared = x * x + y * y;
          ++repeats;
      }
      image_row[i] = repeats;
  }
}

typedef struct thread_data {
    int num_threads;
    int size;
    int width;
    int height;
    int* image;
    int waiting;
    double communication_time;
    double waiting_time;
    volatile int* ready_to_write;
} thread_data;

int next_row_global;
pthread_mutex_t mutex;

void* recv_t(void* thread_d) {
    thread_data* td = (thread_data*)thread_d;
    MPI_Status status;
    int source;
    int j;
    int next_row;
    double begin, end;
    while (td->waiting) {
        begin = MPI_Wtime();
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        end = MPI_Wtime();
        td->waiting_time += (end - begin);
        
        source = status.MPI_SOURCE;
        j = status.MPI_TAG;

        begin = MPI_Wtime();
        MPI_Recv(&td->image[j * td->width], td->width, MPI_INT, source, j, MPI_COMM_WORLD, &status);
        end = MPI_Wtime();
        td->communication_time += (end - begin);

        td->ready_to_write[j] = 1;

        pthread_mutex_lock(&mutex);
        next_row = next_row_global;
        next_row_global--;
        pthread_mutex_unlock(&mutex);

        begin = MPI_Wtime();
        MPI_Send(&next_row, 1, MPI_INT, source, 123, MPI_COMM_WORLD);
        end = MPI_Wtime();
        td->communication_time += (end - begin);
        
        if (next_row < 0) {
            td->waiting--;
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

    int rank, size, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double begin, end;
    double computing_time, communication_time, waiting_time, write_time = 0;

    next_row_global = (height-1) - (size-1);

    if (rank == size-1) {
        /* allocate memory for image */
        int* image = (int*)malloc(width * height * sizeof(int));
        assert(image);

        volatile int ready_to_write[height];
        for (int i = 0; i < height; i++) {
            ready_to_write[i] = 0;
        }

        /* use a thread to receive message */
        pthread_t thread;
        thread_data td;
        td.size = size;
        td.width = width;
        td.height = height;
        td.image = image;
        td.waiting = size - 1;
        td.communication_time = 0;
        td.waiting_time = 0;
        td.ready_to_write = ready_to_write;
        int rc;
        rc = pthread_create(&thread, NULL, recv_t, (void*)&td);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
       
        /* use a thread to keep writng png */
        pthread_t thread_wrt;
        write_png_data wd;
        wd.filename = filename;
        wd.iters = iters;
        wd.width = width;
        wd.height = height;
        wd.buffer = image;
        wd.ready_to_write = ready_to_write;
        rc = pthread_create(&thread_wrt, NULL, write_png_t, (void*)&wd);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
       
        /* the rest: do computing */ 
        pthread_mutex_lock(&mutex);
        int next_row = next_row_global;
        next_row_global--;
        pthread_mutex_unlock(&mutex);
        while (next_row >= 0) { 
            begin = MPI_Wtime();
            mandelbrot_set(rank, size, iters, left, right, lower, upper, width, height, &image[next_row * width], CPU_COUNT(&cpu_set)-2, next_row);
            end = MPI_Wtime();
            computing_time += (end - begin);
            
            ready_to_write[next_row] = 1;

            pthread_mutex_lock(&mutex);
            next_row = next_row_global;
            next_row_global--;
            pthread_mutex_unlock(&mutex);
        }

        pthread_join(thread, NULL);
        communication_time = td.communication_time;
        waiting_time = td.waiting_time;
        
        begin = MPI_Wtime();
        pthread_join(thread_wrt, NULL);
        end = MPI_Wtime();
        write_time += (end - begin);

        free(image);
        printf("[%2d] Communication time: %f, Computing time: %f, Extra write time: %f, Waiting time: %f\n", rank, communication_time, computing_time, write_time, waiting_time);
    } else {
        /* allocate memory for image */
        int* image_row = (int*)malloc(width * sizeof(int));
        assert(image_row);

        /* mandelbrot set */
        MPI_Status status;
        int next_row = (height - 1) - rank;
        while (next_row >= 0) {
            begin = MPI_Wtime();
            mandelbrot_set(rank, size, iters, left, right, lower, upper, width, height, image_row, CPU_COUNT(&cpu_set), next_row);
            end = MPI_Wtime();
            computing_time += (end - begin);

            begin = MPI_Wtime();
            MPI_Send(image_row, width, MPI_INT, size-1, next_row, MPI_COMM_WORLD);
            MPI_Recv(&next_row, 1, MPI_INT, size-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            end = MPI_Wtime();
            communication_time += (end - begin);
        }

        free(image_row);
        printf("[%2d] Communication time: %f, Computing time: %f\n", rank, communication_time, computing_time);
    }
    MPI_Finalize();
}

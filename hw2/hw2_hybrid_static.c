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

void mandelbrot_set(int rank, int size, int iters, double left, double right, double lower, double upper, int width, int height, int* image_row, int num_threads, double* computing, double* communication) {
  double begin, end;
  double begin2, end2;
  begin = MPI_Wtime(); 
  for (int j = rank; j < height; j += size) {
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
      if (rank != size-1) {
          begin2 = MPI_Wtime();
          MPI_Send(image_row, width, MPI_INT, size-1, j, MPI_COMM_WORLD);
          end2 = MPI_Wtime();
          *communication += (end2 - begin2);
      } else {
          image_row += (size * width); // last rank is master
      }
  }
  end = MPI_Wtime();
  *computing += ((end - begin) - *communication);
}

typedef struct thread_data {
    int num_threads;
    int size;
    int width;
    int height;
    int* image;
    double communication_time;
} thread_data;

void* recv_t(void* thread_d) {
    thread_data* td = (thread_data*)thread_d;
    MPI_Status status;
    int source;
    double begin, end;
    for (int j = 0; j < td->height; j++) {
        source = j % td->size;
        if (source != (td->size - 1)) { // recv message except for the last rank (master)
            begin = MPI_Wtime();
            MPI_Recv(&td->image[j * td->width], td->width, MPI_INT, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            end = MPI_Wtime();
            td->communication_time += (end - begin);
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
    double computing_time, communication_time, write_time = 0;

    if (rank == size-1) {
        /* allocate memory for image */
        int* image = (int*)malloc(width * height * sizeof(int));
        assert(image);

        /* use a thread to receive message */
        pthread_t thread;
        thread_data td;
        td.size = size;
        td.width = width;
        td.height = height;
        td.image = image;
        td.communication_time = 0;
        int rc;
        rc = pthread_create(&thread, NULL, recv_t, (void*)&td);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
        
        mandelbrot_set(rank, size, iters, left, right, lower, upper, width, height, &image[(size-1) * width], CPU_COUNT(&cpu_set)-1, &computing_time, &communication_time);
        
        pthread_join(thread, NULL);
        communication_time = td.communication_time;

        /* draw and cleanup */
        begin = MPI_Wtime();
        write_png(filename, iters, width, height, image);
        end = MPI_Wtime();
        write_time += (end - begin);
        free(image);
    } else {
        /* allocate memory for image */
        int* image_row = (int*)malloc(width * sizeof(int));
        assert(image_row);

        /* mandelbrot set */
        mandelbrot_set(rank, size, iters, left, right, lower, upper, width, height, image_row, CPU_COUNT(&cpu_set), &computing_time, &communication_time);

        free(image_row);
    }
    MPI_Finalize();
    printf("[%2d] Communication time: %f,Computing time: %f ,Write time: %f\n", rank, communication_time, computing_time, write_time);
}

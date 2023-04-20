#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    int num_evals = atoi(argv[1]);

    int in = 0;

    double start_time = omp_get_wtime();

    for (int i = 0; i < num_evals; i++) {
        double x = (double) rand() / RAND_MAX;
        double y = (double) rand() / RAND_MAX;

        if (sqrt(x*x + y*y) < 1.0) in++;
    }

    double end_time = omp_get_wtime();

    double pi = (in / (double) num_evals) * 4.0;

    printf("The value of pi is: %.10lf (diff: %.10lf)\n", pi, fabs(M_PI - pi));
    printf("The calculation took: %.10lf seconds\n", end_time - start_time);

}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// Modified from Numerical Recipes Page 340 http://numerical.recipes/book/book.html
typedef struct {
    unsigned long long u;
    unsigned long long v;
    unsigned long long w;
} random_generator;

random_generator init_random_generator(unsigned long long seed) {
    random_generator r;
    r.v = 4101842887655102017LL;
    r.w = 1;
    r.u = seed ^ r.v;
    r.v = r.u;
    r.w = r.v;
    return r;
}

double random_double(random_generator *r) {
    //Return 64-bit random integer.
    r->u = r->u * 2862933555777941757LL + 7046029254386353087LL;
    r->v ^= r->v >> 17;
    r->v ^= r->v << 31;
    r->v ^= r->v >> 8;
    r->w = 4294957665U*(r->w & 0xffffffff) + (r->w >> 32);
    unsigned long long x = r->u ^ (r->u << 21);
    x ^= x >> 35;
    x ^= x << 4;
    return 5.42101086242752217E-20 * ((x + r->v) ^ r->w);
}

int main(int argc, char *argv[]) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    random_generator my_random = init_random_generator(time(NULL) + rank);
    double x = random_double(&my_random);

    double tol = 1e-7;
    double pi = 0.0;
    long in = 0;
    long all_in = 0l;
    long total = 0l;
    long all_total = 0l;
    long round = 2454631l;
    long rounds = 0l;

    while (fabs(pi - M_PI) > tol) {
        rounds++;
        for (long i = 0; i < round; i++) {
            double x = random_double(&my_random);
            double y = random_double(&my_random);

            if (sqrt(x*x + y*y) < 1.0) in++;
        }
        total += round;

        /* We could reduce over the total number in and the total number thrown */
        //MPI_Allreduce(&in, &all_in, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        //MPI_Allreduce(&total, &all_total, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

        /* or we could take the average of all of the Pi values found */

        pi = (in / (double) total) * 4.0;
        double avg_pi = 0;;
        MPI_Allreduce(&pi, &avg_pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        pi = avg_pi / size;

        if (rank == 0) {
            printf("Round %ld, and the average calculated value of pi is: %.10lf (diff: %.10lf)\n", rounds, pi, fabs(M_PI - pi));
        }
    }

    MPI_Finalize();
}

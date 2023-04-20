#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

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
    int num_evals = atoi(argv[1]);


    int in = 0;

    double start_time = omp_get_wtime();

    #pragma omp parallel
    {
        random_generator my_random = init_random_generator(omp_get_thread_num() * time(NULL));

        #pragma omp for reduction(+:in)
        for (int i = 0; i < num_evals; i++) {
            double x = random_double(&my_random);
            double y = random_double(&my_random);

            if (sqrt(x*x + y*y) < 1.0) in++;
        }
    }

    double end_time = omp_get_wtime();

    double pi = (in / (double) num_evals) * 4.0;

    printf("The value of pi is: %.10lf (diff: %.10lf)\n", pi, fabs(M_PI - pi));
    printf("The calculation took: %.10lf seconds\n", end_time - start_time);

}

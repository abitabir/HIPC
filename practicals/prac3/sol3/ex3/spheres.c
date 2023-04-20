#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define MAX_RAND 10000

struct timespec timer;

struct sphere_t {
	double x;
	double y;
	double z;
	double r;
};

double random_number() {
	// returns a random floating point number between 0.0 and RAND_MAX
	return fmod(rand() * ((double) rand() / RAND_MAX), MAX_RAND);
}

double get_time() {
	clock_gettime(CLOCK_MONOTONIC, &timer); 
	return (double) (timer.tv_sec + timer.tv_nsec / 1000000000.0);
}

int main(int argc, char *argv[]) {
	double total_time = get_time();
	double fill_time;
	double calc_area_time;
	double calc_vol_time;
	double calc_intersect_time;


	// read N from the first command line argument
	int N = atoi(argv[1]);

	fill_time = get_time();
	struct sphere_t * sphere = malloc(sizeof(struct sphere_t) * N);
	// fill with random numbers
	for (int i = 0; i < N; i++) {
		sphere[i].x = random_number();
		sphere[i].y = random_number();
		sphere[i].z = random_number();
		sphere[i].r = random_number() / 4.0;
	}
	fill_time = get_time() - fill_time;

	// calculate areas
	calc_area_time = get_time();
	double * area = malloc(N * sizeof(double));
	double four_pi = 4.0 * M_PI;
	for (int i = 0; i < N; i++) {
		area[i] = four_pi * sphere[i].r * sphere[i].r;
	}
	calc_area_time = get_time() - calc_area_time;

	// calculate volume
	calc_vol_time = get_time();
	double * volume = malloc(N * sizeof(double));
	double four_over_three_pi = (4.0 / 3.0) * M_PI;
	for (int i = 0; i < N; i++) {
		volume[i] = four_over_three_pi * sphere[i].r * sphere[i].r * sphere[i].r;
	}
	calc_vol_time = get_time() - calc_vol_time;

	// calculate the number of spheres each sphere intersects
	calc_intersect_time = get_time();
	int * intersects = calloc(N, sizeof(int));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			//if (i == j) continue; // same circle

			// calculate distance between two spheres
			double d = sqrt((sphere[j].x - sphere[i].x) * (sphere[j].x - sphere[i].x) 
					      + (sphere[j].y - sphere[i].y) * (sphere[j].y - sphere[i].y) 
					      + (sphere[j].z - sphere[i].z) * (sphere[j].z - sphere[i].z));
			// if the distance is less than the sum of the radii, they intersect
			if (d < (sphere[j].r + sphere[i].r)) {
				intersects[i]++;
				intersects[j]++;
			}
		}
	}
	calc_intersect_time = get_time() - calc_intersect_time;

	// print out all information to the screen (consider piping this to a file)
	printf("x, y, z, r, area, volume, intersects\n");
	for (int i = 0; i < N; i++) {
		printf("%lf, %lf, %lf, %lf, %lf, %lf, %d\n", sphere[i].x, sphere[i].y, sphere[i].z, sphere[i].r, area[i], volume[i], intersects[i]);
	}

	total_time = get_time() - total_time;

	fprintf(stderr, "Timing Summary\n");
	fprintf(stderr, " Fill Time: %lf\n", fill_time);
	fprintf(stderr, " Calc Area Time: %lf\n", calc_area_time);
	fprintf(stderr, " Calc Vol Time: %lf\n", calc_vol_time);
	fprintf(stderr, " Calc Intersect Time: %lf\n", calc_intersect_time);
	fprintf(stderr, " Total Time: %lf\n", total_time);
}


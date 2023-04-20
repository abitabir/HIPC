#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <string.h>

#define MAX_RAND 10000

struct timespec timer;

struct sphere_t {
	double *x;
	double *y;
	double *z;
	double *r;
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
	double calc_small_time;


	// read N from the first command line argument
	int N = atoi(argv[1]);

	fill_time = get_time();
	struct sphere_t sphere;
	sphere.x = _mm_malloc(sizeof(double) * N, 64);
	sphere.y = _mm_malloc(sizeof(double) * N, 64);
	sphere.z = _mm_malloc(sizeof(double) * N, 64);
	sphere.r = _mm_malloc(sizeof(double) * N, 64);
	// fill with random numbers
	for (int i = 0; i < N; i++) {
		sphere.x[i] = random_number();
		sphere.y[i] = random_number();
		sphere.z[i] = random_number();
		sphere.r[i] = random_number() / 4.0;
	}
	fill_time = get_time() - fill_time;

	// calculate areas
	calc_area_time = get_time();
	double * area = _mm_malloc(N * sizeof(double), 64);
	double four_pi = 4.0 * M_PI;
	for (int i = 0; i < N; i++) {
		area[i] = four_pi * sphere.r[i] * sphere.r[i];
	}
	calc_area_time = get_time() - calc_area_time;

	// calculate volume
	calc_vol_time = get_time();
	double * volume = _mm_malloc(N * sizeof(double), 64);
	double four_over_three_pi = (4.0 / 3.0) * M_PI;
	for (int i = 0; i < N; i++) {
		volume[i] = four_over_three_pi * sphere.r[i] * sphere.r[i] * sphere.r[i];
	}
	calc_vol_time = get_time() - calc_vol_time;

	// calculate smallest sphere
	calc_small_time = get_time();
	double smallest = MAXFLOAT;
	#pragma parallel for reduce(min:smallest)
	for (int i = 0; i < N; i++) {
		if (volume[i] < smallest) {
			smallest = volume[i];
		}
	}
	calc_small_time = get_time() - calc_small_time;

	// calculate the number of spheres each sphere intersects
	calc_intersect_time = get_time();
	// make this a long just so its 64 bits like a double
	long * intersects = _mm_malloc(N * sizeof(long), 64);
	memset(intersects, 0, N); // zero the array (since it is no longer calloced)

	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		int j = 0;

		// load sphere[i] into 4 blocks
		__m256d i_x_v = _mm256_broadcast_sd(sphere.x+i);
		__m256d i_y_v = _mm256_broadcast_sd(sphere.y+i);
		__m256d i_z_v = _mm256_broadcast_sd(sphere.z+i);
		__m256d i_r_v = _mm256_broadcast_sd(sphere.r+i);

		//__m256i i_acc = {0, 0, 0, 0};

		for (j = 0; j+4-1 < i; j += 4) {
			// calculate distance between two spheres
			
			// load 4 values from sphere[j]
			__m256d j_x_v = _mm256_load_pd(sphere.x+j);
			__m256d j_y_v = _mm256_load_pd(sphere.y+j);
			__m256d j_z_v = _mm256_load_pd(sphere.z+j);
			__m256d j_r_v = _mm256_load_pd(sphere.r+j);

			// Calculate distance between each sphere[j] amd sphere[i] (subtraction phase)
			__m256d j_minus_i_x = _mm256_sub_pd(j_x_v, i_x_v);
			__m256d j_minus_i_y = _mm256_sub_pd(j_y_v, i_y_v);
			__m256d j_minus_i_z = _mm256_sub_pd(j_z_v, i_z_v);

			// Calculate distance between each sphere[j] amd sphere[i] (pow(2) phase)
			__m256d j_m_i_x_sq = _mm256_mul_pd(j_minus_i_x, j_minus_i_x);
			__m256d j_m_i_y_sq = _mm256_mul_pd(j_minus_i_y, j_minus_i_y);
			__m256d j_m_i_z_sq = _mm256_mul_pd(j_minus_i_z, j_minus_i_z);

			// Calculate distance between each sphere[j] amd sphere[i] (addition phase)
			__m256d x_plus_y = _mm256_add_pd(j_m_i_x_sq, j_m_i_y_sq);
			__m256d x_p_y_p_z = _mm256_add_pd(x_plus_y, j_m_i_z_sq); 

			// Calculate distance between each sphere[j] amd sphere[i] (sqrt phase)
			__m256d d_sqrt = _mm256_sqrt_pd(x_p_y_p_z);

			/* The next two lines are only used for the branching version that is commented out */
			//double d[4];
			//_mm256_store_pd(d, d_sqrt);

			// Add the radii of sphere i and j together
			__m256d j_plus_i_r = _mm256_add_pd(j_r_v, i_r_v);

			/* The next two lines are only used for the branching version that is commented out */
			//double r[4];
			//_mm256_store_pd(r, j_plus_i_r);

			// Check if the distance between two spheres is less than the sum of the radii
			__m256d d_lt_r = _mm256_cmp_pd(d_sqrt, j_plus_i_r, _CMP_LT_OQ);

			// The result above will be 0xFFFFFFFFFFFFFFFF if less than, 0 otherwise
			// Lets use this as a mask with a vector with {1, 1, 1, 1} in.
			// This will give us a long int vector of 1s to add to count intersections for j
			__m256i mask = _mm256_castpd_si256(d_lt_r);
			__m256i ones = {1, 1, 1, 1};
			__m256i j_acc = _mm256_and_si256(ones, mask);

			// Load the intersects for [j to j+3] into a long integer vector, add accumulator, store result 
			__m256i j_intersects = _mm256_load_si256((const __m256i *) &(intersects[j]));
			__m256i j_i_p_acc = _mm256_add_epi64(j_intersects, j_acc);
			_mm256_store_si256((__m256i *) &(intersects[j]), j_i_p_acc);
			
			// For the intersects for i, we have to add 0, 1, 2, 3 or 4
			// The easiest way to do this is just extract the 1s from the 
			// j_acc register, and sum them.
			intersects[i] += _mm256_extract_epi64(j_acc, 0)
			                 + _mm256_extract_epi64(j_acc, 1)
			                 + _mm256_extract_epi64(j_acc, 2)
							 + _mm256_extract_epi64(j_acc, 3);
			/* 
			// NOTE: This code below is actually slightly faster than the AVX2 version of above..!
			
			// if the distance is less than the sum of the radii, they intersect
			if (d[0] < r[0]) {
				intersects[i]++;
				intersects[j]++;
			}
			if (d[1] < r[1]) {
				intersects[i]++;
				intersects[j+1]++;
			}
			if (d[2] < r[2]) {
				intersects[i]++;
				intersects[j+2]++;
			}
			if (d[3] < r[3]) {
				intersects[i]++;
				intersects[j+3]++;
			}
			*/
		}
		for (; j < i; j++) { // loop to cleanup the last few iterations when N is not nicely divisible
			// calculate distance between two spheres
			double d = sqrt((sphere.x[j] - sphere.x[i]) * (sphere.x[j] - sphere.x[i]) 
					      + (sphere.y[j] - sphere.y[i]) * (sphere.y[j] - sphere.y[i]) 
					      + (sphere.z[j] - sphere.z[i]) * (sphere.z[j] - sphere.z[i]));
			// if the distance is less than the sum of the radii, they intersect
			if (d < (sphere.r[j] + sphere.r[i])) {
				intersects[i]++;
				intersects[j]++;
			}
		}
	}
	calc_intersect_time = get_time() - calc_intersect_time;

	// print out all information to the screen (consider piping this to a file)
	printf("x, y, z, r, area, volume, intersects\n");
	for (int i = 0; i < N; i++) {
		printf("%lf, %lf, %lf, %lf, %lf, %lf, %ld\n", sphere.x[i], sphere.y[i], sphere.z[i], sphere.r[i], area[i], volume[i], intersects[i]);
	}
	fprintf(stderr, " Smallest sphere volume: %lf\n", smallest);

	total_time = get_time() - total_time;

	fprintf(stderr, "Timing Summary\n");
	fprintf(stderr, " Fill Time: %lf\n", fill_time);
	fprintf(stderr, " Calc Area Time: %lf\n", calc_area_time);
	fprintf(stderr, " Calc Vol Time: %lf\n", calc_vol_time);
	fprintf(stderr, " Calc Intersect Time: %lf\n", calc_intersect_time);
	fprintf(stderr, " Total Time: %lf\n", total_time);
}


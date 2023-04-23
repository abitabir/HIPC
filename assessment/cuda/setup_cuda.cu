#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "data.h"

/**
 * @brief Allocate a managed 1d array.
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @return double** A 2D array
 */
double *alloc_1d_array(int m, int n) {
  	double *x;
    cudaMallocManaged(&x, m*n * sizeOf(double));
	return x;
}

/**
 * @brief Allocate a managed char array
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @return char** A array
 */
char *alloc_1d_char_array(int m, int n) {
  	char *x;
    cudaMallocManaged(&x, m*n * sizeOf(char));
	return x;
}

/**
 * @brief Free a managed array
 * 
 * @param array The array to free
 */
void free_1d_array(void * array) {
	cudaFree(array);
}

__global__ void problem_set_up_1(double * u, double * v, double * p, double ui, double vi, double nought, int iters, int ydim) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int initial_index = i * ydim + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        u[index] = ui;
		v[index] = vi;
		p[index] = nought;
    }    
}

__global__ void problem_set_up_2(double * flag, double mx, double my, double delx, double dely, double rad1, int jmax, int iters) {
	double x, y;
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index <= iters; index += stride) {
        x = (i - 0.5) * delx - mx;
		y = (j - 0.5) * dely - my;
		flag[index] = (x*x + y*y <= rad1*rad1) ? C_B : C_F;
    }    
}


__global__ void problem_set_up_3(char * flag, int jmax) {
	int i = threadIdx.x;
	flag[i * jmax]      = C_B;
	int index = i * jmax + jmax + 1;
	flag[index] = C_B;
}

	
__global__ void problem_set_up_4(char * flag, int imax, int jmax) {
	int j = threadIdx.x + 1;
	int index = (imax + 1) * jmax + j;
	flag[j]      = C_B;
	flag[index] = C_B;
}

__global__ void problem_set_up_5(char * flag, int jmax, int iters, int * fluid_cells) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    int index_left, index_right, index_up, index_down;
    for (int index = initial_index; index <= iters; index += stride) {
		if (!(flag[index] & C_F)) {
			atomicSub(fluid_cells, 1);;
			index_left = (i-1) * jmax + j;
			index_right = (i+1) * jmax + j;
			index_down = i * jmax + (j-1);
			index_up = i * jmax + (j+1);
			if (flag[index_left] & C_F) flag[index] |= B_W;
			if (flag[index_right] & C_F) flag[index] |= B_E;
			if (flag[index_down] & C_F) flag[index] |= B_S;
			if (flag[index_up] & C_F) flag[index] |= B_N;
		}
	}
}


#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "data.h"

__global__ void boundary_1(int imax, int jmax) {
	int j = threadIdx.x;
    
    /* Fluid freely flows in from the west */
    index_1 = j; // index_1 = 0 * jmax + j;
    index_2 = jmax + j; // index_2 = 1 * jmax + j;
    u[index_1] = u[index_2];
    v[index_1] = v[index_2];

    /* Fluid freely flows out to the east */
    
    index_1 = imax * jmax + j;
    index_2 = (imax-1) * jmax + j;
    u[index_1] = u[index_2];
    v[index_2+jmax+jmax] = v[index_1];  // index_2 = (imax+1) * jmax + j;

}

__global__ void boundary_2(int jmax) {
	int i = threadIdx.x;
    index_1 = i * jmax + jmax;
    v[index_1] = 0.0;
    u[index_1+1] = u[index_1];
    index_1 =  i * jmax;  // index_1 =  i * jmax + 0;
    v[index_1] = 0.0;
    u[index_1] = u[index_1+1];  // index_1+1 =  i * jmax + jmax+1;
}


__global__ void boundary_3(int imax, int jmax, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int initial_index = i * imax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index_1 = initial_index; index_1 < iters; index_1 += stride) {
        if (flag[index_1] & B_NSEW) {
            switch (flag[index_1]) {
                case B_N: 
                    v[index_1]   = 0.0;
                    u[index_1]   = -u[index_1+1];
                    u[index_1-jmax] = -u[index_1-jmax+1];
                    break;
                case B_E: 
                    u[index_1]   = 0.0;
                    v[index_1]   = -v[index_1+jmax];
                    v[index_1-1] = -v[index_1+jmax-1];
                    break;
                case B_S:
                    v[index_1-1] = 0.0;
                    u[index_1]   = -u[index_1-1];
                    u[index_1-jmax] = -u[index_1-jmax-1];
                    break;
                case B_W: 
                    u[index_1-jmax] = 0.0;
                    v[index_1]   = -v[index_1-jmax];
                    v[index_1-1] = -v[index_1-jmax-1];
                    break;
                case B_NE:
                    v[index_1]   = 0.0;
                    u[index_1]   = 0.0;
                    v[index_1-1] = -v[index_1+jmax-1];
                    u[index_1-jmax] = -u[index_1-jmax+1];
                    break;
                case B_SE:
                    v[index_1-1] = 0.0;
                    u[index_1]   = 0.0;
                    v[index_1]   = -v[index_1+jmax];
                    u[index_1-jmax] = -u[index_1-jmax-1];
                    break;
                case B_SW:
                    v[index_1-1] = 0.0;
                    u[index_1-jmax] = 0.0;
                    v[index_1]   = -v[index_1-jmax];
                    u[index_1]   = -u[index_1-1];
                    break;
                case B_NW:
                    v[index_1]   = 0.0;
                    u[index_1-jmax] = 0.0;
                    v[index_1-1] = -v[index_1-jmax-1];
                    u[index_1]   = -u[index_1+1];
                    break;
            }
        }
    }    
}

__global__ void boundary_4(int jmax, double ui, double vi) {
    int j = threadIdx.x;
    u[j] = ui;  // index_1 =  0 * jmax + j;
    v[j] = 2 * vi - v[jmax + j];  // index_1 = 1 * jmax + j;
}
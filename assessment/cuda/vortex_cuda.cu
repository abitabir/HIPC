#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "data.h"

__global__ void vortex_1(int jmax, double Re, double dely, double delx, double del_t, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        if ((flag[index] & C_F) && (flag[index+jmax] & C_F)) {
            double du2dx = ((u[index] + u[index+jmax]) * (u[index] + u[index+jmax]) +
                            y * fabs(u[index] + u[index+jmax]) * (u[index] - u[index+jmax]) -
                            (u[index-jmax] + u[index]) * (u[index-jmax] + u[index]) -
                            y * fabs(u[index-jmax] + u[index]) * (u[index-jmax]-u[index]))
                            / (4.0 * delx);
            double duvdy = ((v[index] + v[index+jmax]) * (u[index] + u[index+1]) +
                            y * fabs(v[index] + v[index+jmax]) * (u[index] - u[index+1]) -
                            (v[index-1] + v[index+jmax-1]) * (u[index-1] + u[index]) -
                            y * fabs(v[index-1] + v[index+jmax-1]) * (u[index-1] - u[index]))
                            / (4.0 * dely);
            double laplu = (u[index+jmax] - 2.0 * u[index] + u[index-jmax]) / delx / delx +
                            (u[index+1] - 2.0 * u[index] + u[index-1]) / dely / dely;

            f[index] = u[index] + del_t * (laplu / Re - du2dx - duvdy);
        } else {
            f[index] = u[index];
        }
    }
}

__global__ void vortex_2(int jmax, double Re, double dely, double delx, double del_t, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
    /* only if both adjacent cells are fluid cells */
        if ((flag[index] & C_F) && (flag[index+1] & C_F)) {
            double duvdx = ((u[index] + u[index+1]) * (v[index] + v[index+jmax]) +
                            y * fabs(u[index] + u[index+1]) * (v[index] - v[index+jmax]) -
                            (u[index-jmax] + u[index-jmax+1]) * (v[index-jmax] + v[index]) -
                            y * fabs(u[index-jmax] + u[index-jmax+1]) * (v[index-jmax]-v[index]))
                            / (4.0 * delx);
            double dv2dy = ((v[index] + v[index+1]) * (v[index] + v[index+1]) +
                            y * fabs(v[index] + v[index+1]) * (v[index] - v[index+1]) -
                            (v[index-1] + v[index]) * (v[index-1] + v[index]) -
                            y * fabs(v[index-1] + v[index]) * (v[index-1] - v[index]))
                            / (4.0 * dely);
            double laplv = (v[index+jmax] - 2.0 * v[index] + v[index-jmax]) / delx / delx +
                            (v[index+1] - 2.0 * v[index] + v[index-1]) / dely / dely;

            g[index] = v[index] + del_t * (laplv / Re - duvdx - dv2dy);
        } else {
            g[index] = v[index];
        }
    }
}

 __global__ void vortex_3<<<1, jmax>>>(imax, jmax) {
    /* f & g at external boundaries */
    int j = threadIdx.x + 1;
    index = imax * jmax + j;
    f[j]    = u[j];
    f[index] = u[index];
 }

 __global__ void vortex_4<<<1, imax>>>(jmax) {
    int i = threadIdx.x + 1;
    index = i * jmax + jmax;
    g[i*jmax]    = v[i*jmax];
    g[index] = v[index];
 }




__global__ void vortex_5(int jmax, double Re, double dely, double delx, double del_t, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        if (flag[index] & C_F) {
            /* only for fluid and non-surface cells */
            rhs[index] = ((f[index] - f[index-jmax]) / delx + 
                            (g[index] - g[index-1]) / dely)
                            / del_t;
        }
    }
}

__global__ void vortex_6(int jmax, double Re, double dely, double delx, double del_t, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        if (flag[index] & C_F) {
            double eps_E = ((flag[index+jmax] & C_F) ? 1.0 : 0.0);
            double eps_W = ((flag[index-jmax] & C_F) ? 1.0 : 0.0);
            double eps_N = ((flag[index+1] & C_F) ? 1.0 : 0.0);
            double eps_S = ((flag[index-1] & C_F) ? 1.0 : 0.0);

            /* only fluid cells */
            double add = (eps_E * (p[index+jmax] - p[index]) - 
                eps_W * (p[index] - p[index-jmax])) * rdx2  +
                (eps_N * (p[index+1] - p[index]) -
                eps_S * (p[index] - p[index-1])) * rdy2  -  rhs[index];
            res[index] = add * add;
        }
    }
}

__global__ void vortex_8(int jmax, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        if (flag[index] & C_F) { p0[index] = p[index] * p[index]; }
    }
}

__global__ void vortex_7(int rb, int jmax, double Re, double dely, double delx, double del_t, double omega, double beta_2, double rdx2, double rdy2, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    // int k = blockIdx.z * blockDim.z + threadIdx.z;  // is 1
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        if ((i + j) % 2 != rb) { continue; }
        if (flag[index] == (C_F | B_NSEW)) {
            /* five point star for interior fluid cells */
            p[index] = (1.0 - omega) * p[index] - 
                    beta_2 * ((p[index+jmax] + p[index-jmax] ) *rdx2
                        + (p[index+1] + p[index-1]) * rdy2
                        - rhs[index]);
        } else if (flag[index] & C_F) { 
            /* modified star near boundary */

            double eps_E = ((flag[index+jmax] & C_F) ? 1.0 : 0.0);
            double eps_W = ((flag[index-jmax] & C_F) ? 1.0 : 0.0);
            double eps_N = ((flag[index+1] & C_F) ? 1.0 : 0.0);
            double eps_S = ((flag[index-1] & C_F) ? 1.0 : 0.0);

            double beta_mod = -omega / ((eps_E + eps_W) * rdx2 + (eps_N + eps_S) * rdy2);
            p[index] = (1.0 - omega) * p[index] -
                beta_mod * ((eps_E * p[index+jmax] + eps_W * p[index-jmax]) * rdx2
                    + (eps_N * p[index+1] + eps_S * p[index-1]) * rdy2
                    - rhs[index]);
        }
    }
}

__global__ void vortex_9<<<grid, block>>>(int jmax, double delx, double del_t, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        /* only if both adjacent cells are fluid cells */
        if ((flag[index] & C_F) && (flag[index+jmax] & C_F)) {
            u[index] = f[index] - (p[index+jmax] - p[index]) * del_t / delx;
        }
    }
}

__global__ void vortex_10<<<grid, block>>>(int jmax, double dely, double del_t, int iters) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    
    int initial_index = i * jmax + j;  // column first index of one dimensional array, as GPUs are column major
    int stride = blockDim.x;

    for (int index = initial_index; index < iters; index += stride) {
        /* only if both adjacent cells are fluid cells */
        if ((flag[index] & C_F) && (flag[index+1] & C_F)) {
            v[index] = g[index] - (p[index+1] - p[index]) * del_t / dely;
        }
    }
}


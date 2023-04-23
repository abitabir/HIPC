#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "vortex.h"
#include "data.h"
#include "vtk.h"
#include "setup.h"
#include "boundary.h"
#include "args.h"
#include "vortex.cu"

/**
 * @brief Computation of tentative velocity field (f, g)
 * 
 */

void compute_tentative_velocity() {
    vortex_1<<<grid, block>>>(jmax, Re, dely, delx, del_t, (imax)*(jmax+1))
    // int index;
    // for (int i = 1; i < imax; i++) {  //! perfect loop
    //     for (int j = 1; j < jmax+1; j++) {
    //         /* only if both adjacent cells are fluid cells */
    //         index = i * jmax + j;
    //         if ((flag[index] & C_F) && (flag[index+jmax] & C_F)) {
    //             double du2dx = ((u[index] + u[index+jmax]) * (u[index] + u[index+jmax]) +
    //                             y * fabs(u[index] + u[index+jmax]) * (u[index] - u[index+jmax]) -
    //                             (u[index-jmax] + u[index]) * (u[index-jmax] + u[index]) -
    //                             y * fabs(u[index-jmax] + u[index]) * (u[index-jmax]-u[index]))
    //                             / (4.0 * delx);
    //             double duvdy = ((v[index] + v[index+jmax]) * (u[index] + u[index+1]) +
    //                             y * fabs(v[index] + v[index+jmax]) * (u[index] - u[index+1]) -
    //                             (v[index-1] + v[index+jmax-1]) * (u[index-1] + u[index]) -
    //                             y * fabs(v[index-1] + v[index+jmax-1]) * (u[index-1] - u[index]))
    //                             / (4.0 * dely);
    //             double laplu = (u[index+jmax] - 2.0 * u[index] + u[index-jmax]) / delx / delx +
    //                             (u[index+1] - 2.0 * u[index] + u[index-1]) / dely / dely;
   
    //             f[index] = u[index] + del_t * (laplu / Re - du2dx - duvdy);
    //         } else {
    //             f[index] = u[index];
    //         }
    //     }
    // }
    vortex_2<<<grid, block>>>(jmax, Re, dely, delx, del_t, (imax+1)*(jmax))
    // for (int i = 1; i < imax+1; i++) {
    //     for (int j = 1; j < jmax; j++) {  //! perfect loop
    //         /* only if both adjacent cells are fluid cells */
    //         index = i * jmax + j;
    //         if ((flag[index] & C_F) && (flag[index+1] & C_F)) {
    //             double duvdx = ((u[index] + u[index+1]) * (v[index] + v[index+jmax]) +
    //                             y * fabs(u[index] + u[index+1]) * (v[index] - v[index+jmax]) -
    //                             (u[index-jmax] + u[index-jmax+1]) * (v[index-jmax] + v[index]) -
    //                             y * fabs(u[index-jmax] + u[index-jmax+1]) * (v[index-jmax]-v[index]))
    //                             / (4.0 * delx);
    //             double dv2dy = ((v[index] + v[index+1]) * (v[index] + v[index+1]) +
    //                             y * fabs(v[index] + v[index+1]) * (v[index] - v[index+1]) -
    //                             (v[index-1] + v[index]) * (v[index-1] + v[index]) -
    //                             y * fabs(v[index-1] + v[index]) * (v[index-1] - v[index]))
    //                             / (4.0 * dely);
    //             double laplv = (v[index+jmax] - 2.0 * v[index] + v[index-jmax]) / delx / delx +
    //                             (v[index+1] - 2.0 * v[index] + v[index-1]) / dely / dely;

    //             g[index] = v[index] + del_t * (laplv / Re - duvdx - dv2dy);
    //         } else {
    //             g[index] = v[index];
    //         }
    //     }
    // }

    vortex_3<<<1, jmax>>>(imax, jmax)
    /* f & g at external boundaries */
    // for (int j = 1; j < jmax+1; j++) {
    //     index = imax * jmax + j;
    //     f[j]    = u[j];
    //     f[index] = u[index];
    // }
    vortex_4<<<1, jmax>>>(imax, jmax)
    // for (int i = 1; i < imax+1; i++) {
    //     index = i * jmax + jmax;
    //     g[i*jmax]    = v[i*jmax];
    //     g[index] = v[index];
    // }

    cudaDeviceSynchronize();
}


/**
 * @brief Calculate the right hand side of the pressure equation 
 * 
 */
void compute_rhs() {
    vortex_5<<<grid, block>>>(jmax, Re, dely, delx, del_t, (imax+1)*(jmax+1))
    // int index;
    // for (int i = 1; i < imax+1; i++) {
    //     for (int j = 1; j < jmax+1; j++) {  //! perfect loop
    //         index = (i * jmax) + j;
    //         if (flag[index] & C_F) {
    //             /* only for fluid and non-surface cells */
    //             rhs[index] = ((f[index] - f[index-jmax]) / delx + 
    //                          (g[index] - g[index-1]) / dely)
    //                          / del_t;
    //         }
    //     }
    // }

    cudaDeviceSynchronize();
}


/**
 * @brief Red/Black SOR to solve the poisson equation.
 * 
 * @return Calculated residual of the computation
 * 
 */
double poisson() {
    double rdx2 = 1.0 / (delx * delx);
    double rdy2 = 1.0 / (dely * dely);
    double beta_2 = -omega / (2.0 * (rdx2 + rdy2));

    double p0 = 0.0;
    /* Calculate sum of squares */
    int index;

    
    double * p0_array, * host_p0;
    cudaMalloc((double *) &p0_array, sizeof(double) * imax * jmax);


    vortex_8<<<grid, block>>>(jmax, (imax+1)*(jmax+1));
    // for (int i = 1; i < imax+1; i++) {  //! perfect loop
    //     for (int j = 1; j < jmax+1; j++) {
    //         index = i * jmax + j;
    //         if (flag[index] & C_F) { p0 += p[index] * p[index]; }
    //     }
    // }

    
    cudaDeviceSynchronize();

    cudaMemcpy(host_p0, res_p0, sizeof(double) * imax * jmax, cudaMemcpyDeviceToHost);

    for (int i = 0; i < imax*jmax; i++) {
        p0 += host_p0[i];
    }

    cudaFree(p0_array);
    free(host_p0);

    p0 = sqrt(p0 / fluid_cells); 
    if (p0 < 0.0001) { p0 = 1.0; }

    /* Red/Black SOR-iteration */
    int iter;
    double res = 0.0;

    dim3 grid_3d, block_3d;
    block_3d.x = 2;  // or 8
    block_3d.y = 128; // or 32  // better to be longer for memory coalesing
    block_3d.z = 1;  // or 8
    grid_3d.x = jmax / block.x;  // 128 / 2 = 64 or 128 / 8 = 16
    grid_3d.y = imax / block.y;
    grid_3d.z = iter / block.z;


    for (int rb = 0; rb < 2; rb++) {
        vortex_7<<<grid_3d, block_3d>>>(rb, imax, jmax, Re, dely, delx, del_t, omega, beta_2, rdx2, rdy2, iters*(imax+1)*(jmax+1));

    // for (iter = 0; iter < itermax; iter++) {
    //     for (int rb = 0; rb < 2; rb++) {
    //         for (int i = 1; i < imax+1; i++) {
    //             for (int j = 1; j < jmax+1; j++) {  //! perfect loop
    //                 index = i * jmax + j;
    //                 if ((i + j) % 2 != rb) { continue; }
    //                 if (flag[index] == (C_F | B_NSEW)) {
    //                     /* five point star for interior fluid cells */
    //                     p[index] = (1.0 - omega) * p[index] - 
    //                           beta_2 * ((p[index+jmax] + p[index-jmax] ) *rdx2
    //                               + (p[index+1] + p[index-1]) * rdy2
    //                               - rhs[index]);
    //                 } else if (flag[index] & C_F) { 
    //                     /* modified star near boundary */

    //                     double eps_E = ((flag[index+jmax] & C_F) ? 1.0 : 0.0);
    //                     double eps_W = ((flag[index-jmax] & C_F) ? 1.0 : 0.0);
    //                     double eps_N = ((flag[index+1] & C_F) ? 1.0 : 0.0);
    //                     double eps_S = ((flag[index-1] & C_F) ? 1.0 : 0.0);

    //                     double beta_mod = -omega / ((eps_E + eps_W) * rdx2 + (eps_N + eps_S) * rdy2);
    //                     p[index] = (1.0 - omega) * p[index] -
    //                         beta_mod * ((eps_E * p[index+jmax] + eps_W * p[index-jmax]) * rdx2
    //                             + (eps_N * p[index+1] + eps_S * p[index-1]) * rdy2
    //                             - rhs[index]);
    //                 }
    //             }
    //         }
    //     }

        double * res_array, * host_res;
        cudaMalloc((double *) &res_array, sizeof(double) * imax * jmax);

        vortex_6<<<grid, block>>>(jmax, Re, dely, delx, del_t, (imax+1)*(jmax+1));
        /* computation of residual */
        // for (int i = 1; i < imax+1; i++) {  //! perfect loop
        //     for (int j = 1; j < jmax+1; j++) {
        //         index = i * jmax + j;
        //         if (flag[index] & C_F) {
        //             double eps_E = ((flag[index+jmax] & C_F) ? 1.0 : 0.0);
        //             double eps_W = ((flag[index-jmax] & C_F) ? 1.0 : 0.0);
        //             double eps_N = ((flag[index+1] & C_F) ? 1.0 : 0.0);
        //             double eps_S = ((flag[index-1] & C_F) ? 1.0 : 0.0);

        //             /* only fluid cells */
        //             double add = (eps_E * (p[index+jmax] - p[index]) - 
        //                 eps_W * (p[index] - p[index-jmax])) * rdx2  +
        //                 (eps_N * (p[index+1] - p[index]) -
        //                 eps_S * (p[index] - p[index-1])) * rdy2  -  rhs[index];
        //             res[index] = add * add;
        //         }
        //     }
        // }
        
        cudaDeviceSynchronize();

        cudaMemcpy(host_res, res_array, sizeof(double) * imax * jmax, cudaMemcpyDeviceToHost);

        for (int i = 0; i < imax*jmax; i++) {
            res += host_res[i];
        }

        cudaFree(res_array);
        free(host_res);

        res = sqrt(res / fluid_cells) / p0;
        
        /* convergence? */
        if (res < eps) break;
    }

    
    cudaDeviceSynchronize();

    return res;
}


/**
 * @brief Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity() {   
    int index;
    
    vortex_9<<<grid, block>>>(jmax, delx, del_t, (imax-2)*(jmax-1));
    // for (int i = 1; i < imax-2; i++) {  //! perfect loop
    //     for (int j = 1; j < jmax-1; j++) {
    //         index = i * jmax + j;
    //         /* only if both adjacent cells are fluid cells */
    //         if ((flag[index] & C_F) && (flag[index+jmax] & C_F)) {
    //             u[index] = f[index] - (p[index+jmax] - p[index]) * del_t / delx;
    //         }
    //     }
    // }
    
    vortex_10<<<grid, block>>>(jmax, dely, del_t, (imax-1)*(jmax-2));
    // for (int i = 1; i < imax-1; i++) {
    //     for (int j = 1; j < jmax-2; j++) {  //! perfect loop
    //         /* only if both adjacent cells are fluid cells */
    //         index = i * jmax + j;
    //         if ((flag[index] & C_F) && (flag[index+1] & C_F)) {
    //             v[index] = g[index] - (p[index+1] - p[index]) * del_t / dely;
    //         }
    //     }
    // }

    
    cudaDeviceSynchronize();
}


/**
 * @brief Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions. Otherwise the simulation becomes unstable.
 */
void set_timestep_interval() {
    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        int index;
        double umax = 1.0e-10;
        double vmax = 1.0e-10; 
        // vortex_1<<<grid, block>>>(flag, jmax, imax*jmax); // better to not parallelise due to overhead - would have to use atomic
        for (int i = 0; i < imax+2; i++) {  // start from j, as when i=0 index=0*jmax+j
            for (int j = 1; j < jmax+2; j++) {  //! perfect loop
                index = i * jmax + j;
                umax = fmax(fabs(u[index]), umax);
            }
        }

        for (int i = 1; i < imax+2; i++) {  // start from i+jmax, as when j=0 index=i*jmax+j
            for (int j = 0; j < jmax+2; j++) {  //! perfect loop
                index = i * jmax + j;
                vmax = fmax(fabs(v[index]), vmax);
            }
        }

        double deltu = delx / umax;
        double deltv = dely / vmax; 
        double deltRe = 1.0 / (1.0 / (delx * delx) + 1 / (dely * dely)) * Re / 2.0;

        if (deltu < deltv) {
            del_t = fmin(deltu, deltRe);
        } else {
            del_t = fmin(deltv, deltRe);
        }
        del_t = tau * del_t; /* multiply by safety factor */
    }
}

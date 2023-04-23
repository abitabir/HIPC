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

/**
 * @brief Computation of tentative velocity field (f, g)
 * 
 */

void compute_tentative_velocity() {
    int index;
    for (int i = 1; i < imax; i++) {  //! perfect loop
        for (int j = 1; j < jmax+1; j++) {
            /* only if both adjacent cells are fluid cells */
            index = i * jmax + j;
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
    for (int i = 1; i < imax+1; i++) {
        for (int j = 1; j < jmax; j++) {  //! perfect loop
            /* only if both adjacent cells are fluid cells */
            index = i * jmax + j;
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

    /* f & g at external boundaries */
    for (int j = 1; j < jmax+1; j++) {
        index = imax * jmax + j;
        f[j]    = u[j];
        f[index] = u[index];
    }
    for (int i = 1; i < imax+1; i++) {
        index = i * jmax + jmax;
        g[i*jmax]    = v[i*jmax];
        g[index] = v[index];
    }
}


/**
 * @brief Calculate the right hand side of the pressure equation 
 * 
 */
void compute_rhs() {
    int index;
    for (int i = 1; i < imax+1; i++) {
        for (int j = 1; j < jmax+1; j++) {  //! perfect loop
            index = (i * jmax) + j;
            if (flag[index] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[index] = ((f[index] - f[index-jmax]) / delx + 
                             (g[index] - g[index-1]) / dely)
                             / del_t;
            }
        }
    }
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
    for (int i = 1; i < imax+1; i++) {  //! perfect loop
        for (int j = 1; j < jmax+1; j++) {
            index = i * jmax + j;
            if (flag[index] & C_F) { p0 += p[index] * p[index]; }
        }
    }
   
    p0 = sqrt(p0 / fluid_cells); 
    if (p0 < 0.0001) { p0 = 1.0; }

    /* Red/Black SOR-iteration */
    int iter;
    double res = 0.0;
    for (iter = 0; iter < itermax; iter++) {
        for (int rb = 0; rb < 2; rb++) {
            for (int i = 1; i < imax+1; i++) {
                for (int j = 1; j < jmax+1; j++) {  //! perfect loop
                    index = i * jmax + j;
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
        }
        
        /* computation of residual */
        for (int i = 1; i < imax+1; i++) {  //! perfect loop
            for (int j = 1; j < jmax+1; j++) {
                index = i * jmax + j;
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
                    res += add * add;
                }
            }
        }
        res = sqrt(res / fluid_cells) / p0;
        
        /* convergence? */
        if (res < eps) break;
    }

    return res;
}


/**
 * @brief Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity() {   
    int index;
    for (int i = 1; i < imax-2; i++) {  //! perfect loop
        for (int j = 1; j < jmax-1; j++) {
            index = i * jmax + j;
            /* only if both adjacent cells are fluid cells */
            if ((flag[index] & C_F) && (flag[index+jmax] & C_F)) {
                u[index] = f[index] - (p[index+jmax] - p[index]) * del_t / delx;
            }
        }
    }
    
    for (int i = 1; i < imax-1; i++) {
        for (int j = 1; j < jmax-2; j++) {  //! perfect loop
            /* only if both adjacent cells are fluid cells */
            index = i * jmax + j;
            if ((flag[index] & C_F) && (flag[index+1] & C_F)) {
                v[index] = g[index] - (p[index+1] - p[index]) * del_t / dely;
            }
        }
    }
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
        for (int i = 0; i < imax+2; i++) {
            for (int j = 1; j < jmax+2; j++) {  //! perfect loop
                index = i * jmax + j;
                umax = fmax(fabs(u[index]), umax);
            }
        }

        for (int i = 1; i < imax+2; i++) {
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

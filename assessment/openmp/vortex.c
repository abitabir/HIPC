#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>

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
    
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < imax; i++) {  //! perfect loop
        for (int j = 1; j < jmax+1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                double du2dx = ((u[i][j] + u[i+1][j]) * (u[i][j] + u[i+1][j]) +
                                y * fabs(u[i][j] + u[i+1][j]) * (u[i][j] - u[i+1][j]) -
                                (u[i-1][j] + u[i][j]) * (u[i-1][j] + u[i][j]) -
                                y * fabs(u[i-1][j] + u[i][j]) * (u[i-1][j]-u[i][j]))
                                / (4.0 * delx);
                double duvdy = ((v[i][j] + v[i+1][j]) * (u[i][j] + u[i][j+1]) +
                                y * fabs(v[i][j] + v[i+1][j]) * (u[i][j] - u[i][j+1]) -
                                (v[i][j-1] + v[i+1][j-1]) * (u[i][j-1] + u[i][j]) -
                                y * fabs(v[i][j-1] + v[i+1][j-1]) * (u[i][j-1] - u[i][j]))
                                / (4.0 * dely);
                double laplu = (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / delx / delx +
                                (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / dely / dely;
   
                f[i][j] = u[i][j] + del_t * (laplu / Re - du2dx - duvdy);
            } else {
                f[i][j] = u[i][j];
            }
        }
    }
    
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < imax+1; i++) {
        for (int j = 1; j < jmax; j++) {  //! perfect loop
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                double duvdx = ((u[i][j] + u[i][j+1]) * (v[i][j] + v[i+1][j]) +
                                y * fabs(u[i][j] + u[i][j+1]) * (v[i][j] - v[i+1][j]) -
                                (u[i-1][j] + u[i-1][j+1]) * (v[i-1][j] + v[i][j]) -
                                y * fabs(u[i-1][j] + u[i-1][j+1]) * (v[i-1][j]-v[i][j]))
                                / (4.0 * delx);
                double dv2dy = ((v[i][j] + v[i][j+1]) * (v[i][j] + v[i][j+1]) +
                                y * fabs(v[i][j] + v[i][j+1]) * (v[i][j] - v[i][j+1]) -
                                (v[i][j-1] + v[i][j]) * (v[i][j-1] + v[i][j]) -
                                y * fabs(v[i][j-1] + v[i][j]) * (v[i][j-1] - v[i][j]))
                                / (4.0 * dely);
                double laplv = (v[i+1][j] - 2.0 * v[i][j] + v[i-1][j]) / delx / delx +
                                (v[i][j+1] - 2.0 * v[i][j] + v[i][j-1]) / dely / dely;

                g[i][j] = v[i][j] + del_t * (laplv / Re - duvdx - dv2dy);
            } else {
                g[i][j] = v[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    
    #pragma omp parallel for
    for (int j = 1; j < jmax+1; j++) {
        f[0][j]    = u[0][j];
        f[imax][j] = u[imax][j];
    }
    
    #pragma omp parallel for
    for (int i = 1; i < imax+1; i++) {
        g[i][0]    = v[i][0];
        g[i][jmax] = v[i][jmax];
    }
}


/**
 * @brief Calculate the right hand side of the pressure equation 
 * 
 */
void compute_rhs() {
    
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < imax+1; i++) {
        for (int j = 1;j < jmax+1; j++) {  //! perfect loop
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = ((f[i][j] - f[i-1][j]) / delx + 
                             (g[i][j] - g[i][j-1]) / dely)
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

    #pragma omp parallel for collapse(2) reduction(+:p0)
    for (int i = 1; i < imax+1; i++) {  //! perfect loop
        for (int j = 1; j < jmax+1; j++) {
            if (flag[i][j] & C_F) { p0 = p[i][j] * p[i][j]; }
        }
    }
   
    p0 = sqrt(p0 / fluid_cells); 
    if (p0 < 0.0001) { p0 = 1.0; }

    /* Red/Black SOR-iteration */
    int iter;
    double res = 0.0;
    
    #pragma omp parallel for collapse(2) reduction(max:vmax)
    for (iter = 0; iter < itermax; iter++) {
        for (int rb = 0; rb < 2; rb++) {
            for (int i = 1; i < imax+1; i++) {
                for (int j = 1; j < jmax+1; j++) {  //! perfect loop
                    if ((i + j) % 2 != rb) { continue; }
                    if (flag[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p[i][j] = (1.0 - omega) * p[i][j] - 
                              beta_2 * ((p[i+1][j] + p[i-1][j] ) *rdx2
                                  + (p[i][j+1] + p[i][j-1]) * rdy2
                                  - rhs[i][j]);
                    } else if (flag[i][j] & C_F) { 
                        /* modified star near boundary */

                        double eps_E = ((flag[i+1][j] & C_F) ? 1.0 : 0.0);
                        double eps_W = ((flag[i-1][j] & C_F) ? 1.0 : 0.0);
                        double eps_N = ((flag[i][j+1] & C_F) ? 1.0 : 0.0);
                        double eps_S = ((flag[i][j-1] & C_F) ? 1.0 : 0.0);

                        double beta_mod = -omega / ((eps_E + eps_W) * rdx2 + (eps_N + eps_S) * rdy2);
                        p[i][j] = (1.0 - omega) * p[i][j] -
                            beta_mod * ((eps_E * p[i+1][j] + eps_W * p[i-1][j]) * rdx2
                                + (eps_N * p[i][j+1] + eps_S * p[i][j-1]) * rdy2
                                - rhs[i][j]);
                    }
                }
            }
        }
        
        /* computation of residual */

        #pragma omp parallel for collapse(2) reduction(+:resz)
        for (int i = 1; i < imax+1; i++) {  //! perfect loop
            for (int j = 1; j < jmax+1; j++) {
                if (flag[i][j] & C_F) {
                    double eps_E = ((flag[i+1][j] & C_F) ? 1.0 : 0.0);
                    double eps_W = ((flag[i-1][j] & C_F) ? 1.0 : 0.0);
                    double eps_N = ((flag[i][j+1] & C_F) ? 1.0 : 0.0);
                    double eps_S = ((flag[i][j-1] & C_F) ? 1.0 : 0.0);

                    /* only fluid cells */
                    double add = (eps_E * (p[i+1][j] - p[i][j]) - 
                        eps_W * (p[i][j] - p[i-1][j])) * rdx2  +
                        (eps_N * (p[i][j+1] - p[i][j]) -
                        eps_S * (p[i][j] - p[i][j-1])) * rdy2  -  rhs[i][j];
                    res = add * add;
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
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < imax-2; i++) {  //! perfect loop
        for (int j = 1; j < jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                u[i][j] = f[i][j] - (p[i+1][j] - p[i][j]) * del_t / delx;
            }
        }
    }

    #pragma omp parallel for collapse(2)
    for (int i = 1; i < imax-1; i++) {
        for (int j = 1; j < jmax-2; j++) {  //! perfect loop
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                v[i][j] = g[i][j] - (p[i][j+1] - p[i][j]) * del_t / dely;
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
        double umax = 1.0e-10;
        double vmax = 1.0e-10; 
            
        #pragma omp parallel for collapse(2) reduction(max:umax)
        for (int i = 0; i < imax+2; i++) {
            for (int j = 1; j < jmax+2; j++) {
                umax = fabs(u[i][j]);
            }
        }

        #pragma omp parallel for collapse(2) reduction(max:vmax)
        for (int i = 1; i < imax+2; i++) {
            for (int j = 0; j < jmax+2; j++) {
                vmax = fabs(v[i][j]);
            }
        }
        
        double deltu = delx / umax;
        double deltv = dely / vmax; 
        double deltRe = 1.0 / (1.0 / (delx * delx) + 1 / (dely * dely)) * Re / 2.0;

        // #pragma omp parallel reduction(min:del_t) {
        //     #pragma omp single
        //     del_t = deltu;
        //     #pragma omp single
        //     del_t = deltv;
        //     #pragma omp single
        //     del_t = deltRe;
        // }

        if (deltu < deltv) {
            del_t = fmin(deltu, deltRe);
        } else {
            del_t = fmin(deltv, deltRe);
        }

        del_t = tau * del_t; /* multiply by safety factor */
    }
}

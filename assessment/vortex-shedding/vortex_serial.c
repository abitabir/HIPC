#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>

#include "data_serial.h"
// #include "vtk.h"
#include "setup_serial.h"
#include "boundary_serial.h"

/**
 * @brief Computation of tentative velocity field (f, g)
 * 
 */
void compute_tentative_velocity_serial() {
    for (int i = 1; i < imax_serial; i++) {  //! perfect loop
        for (int j = 1; j < jmax_serial+1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag_serial[i][j] & C_F) && (flag_serial[i+1][j] & C_F)) {
                double du2dx = ((u_serial[i][j] + u_serial[i+1][j]) * (u_serial[i][j] + u_serial[i+1][j]) +
                                y_serial * fabs(u_serial[i][j] + u_serial[i+1][j]) * (u_serial[i][j] - u_serial[i+1][j]) -
                                (u_serial[i-1][j] + u_serial[i][j]) * (u_serial[i-1][j] + u_serial[i][j]) -
                                y_serial * fabs(u_serial[i-1][j] + u_serial[i][j]) * (u_serial[i-1][j]-u_serial[i][j]))
                                / (4.0 * delx_serial);
                double duvdy = ((v_serial[i][j] + v_serial[i+1][j]) * (u_serial[i][j] + u_serial[i][j+1]) +
                                y_serial * fabs(v_serial[i][j] + v_serial[i+1][j]) * (u_serial[i][j] - u_serial[i][j+1]) -
                                (v_serial[i][j-1] + v_serial[i+1][j-1]) * (u_serial[i][j-1] + u_serial[i][j]) -
                                y_serial * fabs(v_serial[i][j-1] + v_serial[i+1][j-1]) * (u_serial[i][j-1] - u_serial[i][j]))
                                / (4.0 * dely_serial);
                double laplu = (u_serial[i+1][j] - 2.0 * u_serial[i][j] + u_serial[i-1][j]) / delx_serial / delx_serial +
                                (u_serial[i][j+1] - 2.0 * u_serial[i][j] + u_serial[i][j-1]) / dely_serial / dely_serial;
   
                f_serial[i][j] = u_serial[i][j] + del_t_serial * (laplu / Re_serial - du2dx - duvdy);
            } else {
                f_serial[i][j] = u_serial[i][j];
            }
        }
    }
    for (int i = 1; i < imax_serial+1; i++) {
        for (int j = 1; j < jmax_serial; j++) {  //! perfect loop
            /* only if both adjacent cells are fluid cells */
            if ((flag_serial[i][j] & C_F) && (flag_serial[i][j+1] & C_F)) {
                double duvdx = ((u_serial[i][j] + u_serial[i][j+1]) * (v_serial[i][j] + v_serial[i+1][j]) +
                                y_serial * fabs(u_serial[i][j] + u_serial[i][j+1]) * (v_serial[i][j] - v_serial[i+1][j]) -
                                (u_serial[i-1][j] + u_serial[i-1][j+1]) * (v_serial[i-1][j] + v_serial[i][j]) -
                                y_serial * fabs(u_serial[i-1][j] + u_serial[i-1][j+1]) * (v_serial[i-1][j]-v_serial[i][j]))
                                / (4.0 * delx_serial);
                double dv2dy = ((v_serial[i][j] + v_serial[i][j+1]) * (v_serial[i][j] + v_serial[i][j+1]) +
                                y_serial * fabs(v_serial[i][j] + v_serial[i][j+1]) * (v_serial[i][j] - v_serial[i][j+1]) -
                                (v_serial[i][j-1] + v_serial[i][j]) * (v_serial[i][j-1] + v_serial[i][j]) -
                                y_serial * fabs(v_serial[i][j-1] + v_serial[i][j]) * (v_serial[i][j-1] - v_serial[i][j]))
                                / (4.0 * dely_serial);
                double laplv = (v_serial[i+1][j] - 2.0 * v_serial[i][j] + v_serial[i-1][j]) / delx_serial / delx_serial +
                                (v_serial[i][j+1] - 2.0 * v_serial[i][j] + v_serial[i][j-1]) / dely_serial / dely_serial;

                g_serial[i][j] = v_serial[i][j] + del_t_serial * (laplv / Re_serial - duvdx - dv2dy);
            } else {
                g_serial[i][j] = v_serial[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    for (int j = 1; j < jmax_serial+1; j++) {
        f_serial[0][j]    = u_serial[0][j];
        f_serial[imax_serial][j] = u_serial[imax_serial][j];
    }
    for (int i = 1; i < imax_serial+1; i++) {
        g_serial[i][0]    = v_serial[i][0];
        g_serial[i][jmax_serial] = v_serial[i][jmax_serial];
    }
}


/**
 * @brief Calculate the right hand side of the pressure equation 
 * 
 */
void compute_rhs_serial() {
    for (int i = 1; i < imax_serial+1; i++) {
        for (int j = 1;j < jmax_serial+1; j++) {  //! perfect loop
            if (flag_serial[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs_serial[i][j] = ((f_serial[i][j] - f_serial[i-1][j]) / delx_serial + 
                             (g_serial[i][j] - g_serial[i][j-1]) / dely_serial)
                             / del_t_serial;
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
double poisson_serial() {
    double rdx2 = 1.0 / (delx_serial * delx_serial);
    double rdy2 = 1.0 / (dely_serial * dely_serial);
    double beta_2 = -omega_serial / (2.0 * (rdx2 + rdy2));

    double p0 = 0.0;
    /* Calculate sum of squares */
    for (int i = 1; i < imax_serial+1; i++) {  //! perfect loop
        for (int j = 1; j < jmax_serial+1; j++) {
            if (flag_serial[i][j] & C_F) { p0 += p_serial[i][j] * p_serial[i][j]; }
        }
    }
   
    p0 = sqrt(p0 / fluid_cells_serial); 
    if (p0 < 0.0001) { p0 = 1.0; }

    /* Red/Black SOR-iteration */
    int iter;
    double res = 0.0;
    for (iter = 0; iter < itermax_serial; iter++) {
        for (int rb = 0; rb < 2; rb++) {
            for (int i = 1; i < imax_serial+1; i++) {
                for (int j = 1; j < jmax_serial+1; j++) {  //! perfect loop
                    if ((i + j) % 2 != rb) { continue; }
                    if (flag_serial[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p_serial[i][j] = (1.0 - omega_serial) * p_serial[i][j] - 
                              beta_2 * ((p_serial[i+1][j] + p_serial[i-1][j] ) *rdx2
                                  + (p_serial[i][j+1] + p_serial[i][j-1]) * rdy2
                                  - rhs_serial[i][j]);
                    } else if (flag_serial[i][j] & C_F) { 
                        /* modified star near boundary */

                        double eps_E = ((flag_serial[i+1][j] & C_F) ? 1.0 : 0.0);
                        double eps_W = ((flag_serial[i-1][j] & C_F) ? 1.0 : 0.0);
                        double eps_N = ((flag_serial[i][j+1] & C_F) ? 1.0 : 0.0);
                        double eps_S = ((flag_serial[i][j-1] & C_F) ? 1.0 : 0.0);

                        double beta_mod = -omega_serial / ((eps_E + eps_W) * rdx2 + (eps_N + eps_S) * rdy2);
                        p_serial[i][j] = (1.0 - omega_serial) * p_serial[i][j] -
                            beta_mod * ((eps_E * p_serial[i+1][j] + eps_W * p_serial[i-1][j]) * rdx2
                                + (eps_N * p_serial[i][j+1] + eps_S * p_serial[i][j-1]) * rdy2
                                - rhs_serial[i][j]);
                    }
                }
            }
        }
        
        /* computation of residual */
        for (int i = 1; i < imax_serial+1; i++) {  //! perfect loop
            for (int j = 1; j < jmax_serial+1; j++) {
                if (flag_serial[i][j] & C_F) {
                    double eps_E = ((flag_serial[i+1][j] & C_F) ? 1.0 : 0.0);
                    double eps_W = ((flag_serial[i-1][j] & C_F) ? 1.0 : 0.0);
                    double eps_N = ((flag_serial[i][j+1] & C_F) ? 1.0 : 0.0);
                    double eps_S = ((flag_serial[i][j-1] & C_F) ? 1.0 : 0.0);

                    /* only fluid cells */
                    double add = (eps_E * (p_serial[i+1][j] - p_serial[i][j]) - 
                        eps_W * (p_serial[i][j] - p_serial[i-1][j])) * rdx2  +
                        (eps_N * (p_serial[i][j+1] - p_serial[i][j]) -
                        eps_S * (p_serial[i][j] - p_serial[i][j-1])) * rdy2  -  rhs_serial[i][j];
                    res += add * add;
                }
            }
        }
        res = sqrt(res / fluid_cells_serial) / p0;
        
        /* convergence? */
        if (res < eps_serial) break;
    }

    return res;
}


/**
 * @brief Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity_serial() {   
    for (int i = 1; i < imax_serial-2; i++) {  //! perfect loop
        for (int j = 1; j < jmax_serial-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag_serial[i][j] & C_F) && (flag_serial[i+1][j] & C_F)) {
                u_serial[i][j] = f_serial[i][j] - (p_serial[i+1][j] - p_serial[i][j]) * del_t_serial / delx_serial;
            }
        }
    }
    
    for (int i = 1; i < imax_serial-1; i++) {
        for (int j = 1; j < jmax_serial-2; j++) {  //! perfect loop
            /* only if both adjacent cells are fluid cells */
            if ((flag_serial[i][j] & C_F) && (flag_serial[i][j+1] & C_F)) {
                v_serial[i][j] = g_serial[i][j] - (p_serial[i][j+1] - p_serial[i][j]) * del_t_serial / dely_serial;
            }
        }
    }
}


/**
 * @brief Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions. Otherwise the simulation becomes unstable.
 */
void set_timestep_interval_serial() {
    /* del_t_serial satisfying CFL conditions */
    if (tau_serial >= 1.0e-10) { /* else no time stepsize control */
        double umax = 1.0e-10;
        double vmax = 1.0e-10; 
        
        for (int i = 0; i < imax_serial+2; i++) {
            for (int j = 1; j < jmax_serial+2; j++) {  //! perfect loop
                umax = fmax(fabs(u_serial[i][j]), umax);
            }
        }

        for (int i = 1; i < imax_serial+2; i++) {
            for (int j = 0; j < jmax_serial+2; j++) {  //! perfect loop
                vmax = fmax(fabs(v_serial[i][j]), vmax);
            }
        }

        double deltu = delx_serial / umax;
        double deltv = dely_serial / vmax; 
        double deltRe = 1.0 / (1.0 / (delx_serial * delx_serial) + 1 / (dely_serial * dely_serial)) * Re_serial / 2.0;

        if (deltu < deltv) {
            del_t_serial = fmin(deltu, deltRe);
        } else {
            del_t_serial = fmin(deltv, deltRe);
        }
        del_t_serial = tau_serial * del_t_serial; /* multiply by safety factor */
    }
}

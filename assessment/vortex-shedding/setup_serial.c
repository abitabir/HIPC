#include <stdio.h>
#include <stdlib.h>

#include "data_serial.h"
#include "boundary_serial.h"

// /**
//  * @brief Set up some default values before arguments are parsed.
//  * 
//  */
// void set_defaults_serial() {
// 	set_default_base_serial();
// }

/**
 * @brief Set up some values after arguments have been parsed.
 * 
 */
void setup_serial() {
	delx_serial = xlength_serial/imax_serial;
    dely_serial = ylength_serial/jmax_serial;
}
 
/**
 * @brief Allocate all of the arrays used by the computation.
 * 
 */
void allocate_arrays_serial() {
	    /* Allocate arrays */
    u_size_x_serial = imax_serial+2; u_size_y_serial = jmax_serial+2;
    u_serial = alloc_2d_array_serial(u_size_x_serial, u_size_y_serial);
    v_size_x_serial = imax_serial+2; v_size_y_serial = jmax_serial+2;
    v_serial = alloc_2d_array_serial(v_size_x_serial, v_size_y_serial);
    f_size_x_serial = imax_serial+2; f_size_y_serial = jmax_serial+2;
    f_serial = alloc_2d_array_serial(f_size_x_serial, f_size_y_serial);
    g_size_x_serial = imax_serial+2; g_size_y_serial = jmax_serial+2;
    g_serial = alloc_2d_array_serial(g_size_x_serial, g_size_y_serial);
    p_size_x_serial = imax_serial+2; p_size_y_serial = jmax_serial+2;
    p_serial = alloc_2d_array_serial(p_size_x_serial, p_size_y_serial);
    rhs_size_x_serial = imax_serial+2; rhs_size_y_serial = jmax_serial+2;
    rhs_serial = alloc_2d_array_serial(rhs_size_x_serial, rhs_size_y_serial);
    flag_size_x_serial = imax_serial+2; flag_size_y_serial = jmax_serial+2;
    flag_serial = alloc_2d_char_array_serial(flag_size_x_serial, flag_size_y_serial);

    if (!u_serial || !v_serial || !f_serial || !g_serial || !p_serial || !rhs_serial || !flag_serial) {
        fprintf(stderr, "Couldn't allocate memory for matrices.\n");
		exit(1);
    }
}

/**
 * @brief Free all of the arrays used for the computation.
 * 
 */
void free_arrays_serial() {
	free_2d_array_serial((void**) u_serial);
    free_2d_array_serial((void**) v_serial);
    free_2d_array_serial((void**) f_serial);
    free_2d_array_serial((void**) g_serial);
    free_2d_array_serial((void**) p_serial);
    free_2d_array_serial((void**) rhs_serial);
    free_2d_array_serial((void**) flag_serial);
}

/**
 * @brief Initialise the velocity arrays and then initialize the flag array, 
 * marking any obstacle cells and the edge cells as boundaries. The cells 
 * adjacent to boundary cells have their relevant flags set too.
 */
void problem_set_up_serial() {
    for (int i = 0; i < imax_serial+2; i++) {
        for (int j = 0; j < jmax_serial+2; j++) {  //! perfect loop
            u_serial[i][j] = ui_serial;
            v_serial[i][j] = vi_serial;
            p_serial[i][j] = 0.0;
        }
    }

    /* Mark a circular obstacle as boundary cells, the rest as fluid */
    double mx = 20.0 / 41.0 * jmax_serial * dely_serial;
    double my = mx;
    double rad1 = 5.0 / 41.0 * jmax_serial * dely_serial;
    for (int i = 1; i <= imax_serial; i++) {
        for (int j = 1; j <= jmax_serial; j++) {
            double x = (i - 0.5) * delx_serial - mx;
            double y_serial = (j - 0.5) * dely_serial - my;
            flag_serial[i][j] = (x*x + y_serial*y_serial <= rad1*rad1) ? C_B : C_F;
        }
    }
    
    /* Mark the north & south boundary cells */
    for (int i = 0; i <= imax_serial + 1; i++) {
        flag_serial[i][0]      = C_B;
        flag_serial[i][jmax_serial+1] = C_B;
    }
    /* Mark the east and west boundary cells */
    for (int j = 1; j <= jmax_serial; j++) {
        flag_serial[0][j]      = C_B;
        flag_serial[imax_serial+1][j] = C_B;
    }

    fluid_cells_serial = imax_serial * jmax_serial;

    /* flags for boundary cells */
    for (int i = 1; i <= imax_serial; i++) {
        for (int j = 1; j <= jmax_serial; j++) {
            if (!(flag_serial[i][j] & C_F)) {
                fluid_cells_serial--;
                if (flag_serial[i-1][j] & C_F) flag_serial[i][j] |= B_W;
                if (flag_serial[i+1][j] & C_F) flag_serial[i][j] |= B_E;
                if (flag_serial[i][j-1] & C_F) flag_serial[i][j] |= B_S;
                if (flag_serial[i][j+1] & C_F) flag_serial[i][j] |= B_N;
            }
        }
    }

	apply_boundary_conditions_serial();
}
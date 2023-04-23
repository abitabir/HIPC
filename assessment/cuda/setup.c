#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "data.h"
#include "vtk.h"
#include "boundary.h"
#include "setup.cu"

/**
 * @brief Set up some default values before arguments are parsed.
 * 
 */
void set_defaults() {
	set_default_base();
}

/**
 * @brief Set up some values after arguments have been parsed.
 * 
 */
void setup() {
	delx = xlength/imax;
    dely = ylength/jmax;
}
 
/**
 * @brief Allocate all of the arrays used by the computation.
 * 
 */
void allocate_arrays() {
	    /* Allocate arrays */
    u_size_x = imax+2; u_size_y = jmax+2;
    u = alloc_1d_array(u_size_x, u_size_y);
    v_size_x = imax+2; v_size_y = jmax+2;
    v = alloc_1d_array(v_size_x, v_size_y);
    f_size_x = imax+2; f_size_y = jmax+2;
    f = alloc_1d_array(f_size_x, f_size_y);
    g_size_x = imax+2; g_size_y = jmax+2;
    g = alloc_1d_array(g_size_x, g_size_y);
    p_size_x = imax+2; p_size_y = jmax+2;
    p = alloc_1d_array(p_size_x, p_size_y);
    rhs_size_x = imax+2; rhs_size_y = jmax+2;
    rhs = alloc_1d_array(rhs_size_x, rhs_size_y);
    flag_size_x = imax+2; flag_size_y = jmax+2;
    flag = alloc_1d_char_array(flag_size_x, flag_size_y);

    if (!u || !v || !f || !g || !p || !rhs || !flag) {
        fprintf(stderr, "Couldn't allocate memory for matrices.\n");
		exit(1);
    }
}

/**
 * @brief Free all of the arrays used for the computation.
 * 
 */
void free_arrays() {
	free_1d_array((void*) u);
    free_1d_array((void*) v);
    free_1d_array((void*) f);
    free_1d_array((void*) g);
    free_1d_array((void*) p);
    free_1d_array((void*) rhs);
    free_1d_array((void*) flag);
}

/**
 * @brief Initialise the velocity arrays and then initialize the flag array, 
 * marking any obstacle cells and the edge cells as boundaries. The cells 
 * adjacent to boundary cells have their relevant flags set too.
 */
void problem_set_up() {
    problem_set_up_1<<<grid, block>>>(u, v, p, ui, vi, 0.0, (imax+2)*(jmax+2), jmax+2);
    // int index;
    // for (int i = 0; i < imax+2; i++) {
    //     for (int j = 0; j < jmax+2; j++) {  //! perfect loop
    //         index = i * jmax + j;
    //         u[index] = ui;
    //         v[index] = vi;
    //         p[index] = 0.0;
    //     }
    // }

    // don't need to sync here

    /* Mark a circular obstacle as boundary cells, the rest as fluid */
    double mx = 20.0 / 41.0 * jmax * dely;
    double my = mx;
    double rad1 = 5.0 / 41.0 * jmax * dely;

    problem_set_up_2<<<grid, block>>>(flag, mx, my, delx, dely, rad1, jmax, imax*jmax);
   
    // for (int i = 1; i <= imax; i++) {
    //     for (int j = 1; j <= jmax; j++) {
    //         double x = (i - 0.5) * delx - mx;
    //         double y = (j - 0.5) * dely - my;
    //         index = i * jmax + j;
    //         flag[index] = (x*x + y*y <= rad1*rad1) ? C_B : C_F;
    //     }
    // }

    problem_set_up_3<<<1, imax+2>>>(flag, jmax);
    problem_set_up_4<<<1, jmax>>>(flag, imax, jmax);

    // /* Mark the north & south boundary cells */
    // for (int i = 0; i <= imax + 1; i++) {
    //     flag[i * jmax]      = C_B;
    //     index = i * jmax + jmax + 1;
    //     flag[index] = C_B;
    // }
    // /* Mark the east and west boundary cells */
    // for (int j = 1; j <= jmax; j++) {
    //     index = (imax + 1) * jmax + j;
    //     flag[j]      = C_B;
    //     flag[index] = C_B;
    // }

    fluid_cells = imax * jmax;

    int * fluid_cell_ptr = fluid_cells;
    
    cudaMallocManaged(&fluid_cell_ptr, sizeOf(int));

    cudaDeviceSynchronize();
    
    problem_set_up_5<<<grid, block>>>(flag, jmax, imax*jmax, fluid_cell_ptr);

    // int index_left;
    // int index_right;
    // int index_up;
    // int index_down;
    // /* flags for boundary cells */
    // for (int i = 1; i <= imax; i++) {
    //     for (int j = 1; j <= jmax; j++) {
    //         index = i * jmax + j;
    //         if (!(flag[index] & C_F)) {
    //             fluid_cells--;
    //             index_left = (i-1) * jmax + j;
    //             index_right = (i+1) * jmax + j;
    //             index_down = i * jmax + (j-1);
    //             index_up = i * jmax + (j+1);
    //             if (flag[index_left] & C_F) flag[index] |= B_W;
    //             if (flag[index_right] & C_F) flag[index] |= B_E;
    //             if (flag[index_down] & C_F) flag[index] |= B_S;
    //             if (flag[index_up] & C_F) flag[index] |= B_N;
    //         }
    //     }
    // }

	apply_boundary_conditions();
    
    cudaDeviceSynchronize();

    fluid_cells = * fluid_cell_ptr;
}
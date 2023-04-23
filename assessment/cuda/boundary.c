#include <cuda.h>
#include <cuda_runtime.h>

#include "data.h"
#include "boundary.h"
#include "boundary.cu"

/**
 * @brief Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions() {
    boundary_1<<<1, jmax+2>>>(imax, jmax);
    // // index = i * jmax + j;
    // int index_1, index_2;
    // for (int j = 0; j < jmax+2; j++) {
    //     /* Fluid freely flows in from the west */
    //     index_1 = j; // index_1 = 0 * jmax + j;
    //     index_2 = jmax + j; // index_2 = 1 * jmax + j;
    //     u[index_1] = u[index_2];
    //     v[index_1] = v[index_2];

    //     /* Fluid freely flows out to the east */
        
    //     index_1 = imax * jmax + j;
    //     index_2 = (imax-1) * jmax + j;
    //     u[index_1] = u[index_2];
    //     v[index_2+jmax+jmax] = v[index_1];  // index_2 = (imax+1) * jmax + j;
    // }

    boundary_2<<<1, imax+2>>>(jmax);
    // for (int i = 0; i < imax+2; i++) {
    //     /* The vertical velocity approaches 0 at the north and south
    //      * boundaries, but fluid flows freely in the horizontal direction */
    //     index_1 = i * jmax + jmax;
    //     v[index_1] = 0.0;
    //     u[index_1+1] = u[index_1];
    //     index_1 =  i * jmax;  // index_1 =  i * jmax + 0;
    //     v[index_1] = 0.0;
    //     u[index_1] = u[index_1+1];  // index_1+1 =  i * jmax + jmax+1;
    // }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */
    boundary_3<<<grid, block>>>(imax, jmax, (imax+1)*(jmax+1));
    // for (int i = 1; i < imax+1; i++) {
    //     for (int j = 1; j < jmax+1; j++) {
    //         index_1 = i * jmax + j;
    //         if (flag[index_1] & B_NSEW) {
    //             switch (flag[index_1]) {
    //                 case B_N: 
    //                     v[index_1]   = 0.0;
    //                     u[index_1]   = -u[index_1+1];
    //                     u[index_1-jmax] = -u[index_1-jmax+1];
    //                     break;
    //                 case B_E: 
    //                     u[index_1]   = 0.0;
    //                     v[index_1]   = -v[index_1+jmax];
    //                     v[index_1-1] = -v[index_1+jmax-1];
    //                     break;
    //                 case B_S:
    //                     v[index_1-1] = 0.0;
    //                     u[index_1]   = -u[index_1-1];
    //                     u[index_1-jmax] = -u[index_1-jmax-1];
    //                     break;
    //                 case B_W: 
    //                     u[index_1-jmax] = 0.0;
    //                     v[index_1]   = -v[index_1-jmax];
    //                     v[index_1-1] = -v[index_1-jmax-1];
    //                     break;
    //                 case B_NE:
    //                     v[index_1]   = 0.0;
    //                     u[index_1]   = 0.0;
    //                     v[index_1-1] = -v[index_1+jmax-1];
    //                     u[index_1-jmax] = -u[index_1-jmax+1];
    //                     break;
    //                 case B_SE:
    //                     v[index_1-1] = 0.0;
    //                     u[index_1]   = 0.0;
    //                     v[index_1]   = -v[index_1+jmax];
    //                     u[index_1-jmax] = -u[index_1-jmax-1];
    //                     break;
    //                 case B_SW:
    //                     v[index_1-1] = 0.0;
    //                     u[index_1-jmax] = 0.0;
    //                     v[index_1]   = -v[index_1-jmax];
    //                     u[index_1]   = -u[index_1-1];
    //                     break;
    //                 case B_NW:
    //                     v[index_1]   = 0.0;
    //                     u[index_1-jmax] = 0.0;
    //                     v[index_1-1] = -v[index_1-jmax-1];
    //                     u[index_1]   = -u[index_1+1];
    //                     break;
    //             }
    //         }
    //     }
    // }

    /* Finally, fix the horizontal velocity at the  western edge to have
     * a continual flow of fluid into the simulation.
     */
    v[0] = 2 * vi-v[jmax]; // index_1 =  0 * jmax + 0; index_1 =  1 * jmax + 0;

    boundary_4<<<1, jmax>>>(jmax, ui, vi);
    // for (int j = 1; j < jmax+1; j++) {
    //     u[j] = ui;  // index_1 =  0 * jmax + j;
    //     v[j] = 2 * vi - v[jmax + j];  // index_1 = 1 * jmax + j;
    // }

    
    cudaDeviceSynchronize();
}
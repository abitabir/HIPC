#include "data_serial.h"
#include "boundary_serial.h"

/**
 * @brief Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions_serial() {
    for (int j = 0; j < jmax_serial+2; j++) {
        /* Fluid freely flows in from the west */
        u_serial[0][j] = u_serial[1][j];
        v_serial[0][j] = v_serial[1][j];

        /* Fluid freely flows out to the east */
        u_serial[imax_serial][j] = u_serial[imax_serial-1][j];
        v_serial[imax_serial+1][j] = v_serial[imax_serial][j];
    }

    for (int i = 0; i < imax_serial+2; i++) {
        /* The vertical velocity approaches 0 at the north and south
         * boundaries, but fluid flows freely in the horizontal direction */
        v_serial[i][jmax_serial] = 0.0;
        u_serial[i][jmax_serial+1] = u_serial[i][jmax_serial];

        v_serial[i][0] = 0.0;
        u_serial[i][0] = u_serial[i][1];
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */
    for (int i = 1; i < imax_serial+1; i++) {
        for (int j = 1; j < jmax_serial+1; j++) {
            if (flag_serial[i][j] & B_NSEW) {
                switch (flag_serial[i][j]) {
                    case B_N: 
                        v_serial[i][j]   = 0.0;
                        u_serial[i][j]   = -u_serial[i][j+1];
                        u_serial[i-1][j] = -u_serial[i-1][j+1];
                        break;
                    case B_E: 
                        u_serial[i][j]   = 0.0;
                        v_serial[i][j]   = -v_serial[i+1][j];
                        v_serial[i][j-1] = -v_serial[i+1][j-1];
                        break;
                    case B_S:
                        v_serial[i][j-1] = 0.0;
                        u_serial[i][j]   = -u_serial[i][j-1];
                        u_serial[i-1][j] = -u_serial[i-1][j-1];
                        break;
                    case B_W: 
                        u_serial[i-1][j] = 0.0;
                        v_serial[i][j]   = -v_serial[i-1][j];
                        v_serial[i][j-1] = -v_serial[i-1][j-1];
                        break;
                    case B_NE:
                        v_serial[i][j]   = 0.0;
                        u_serial[i][j]   = 0.0;
                        v_serial[i][j-1] = -v_serial[i+1][j-1];
                        u_serial[i-1][j] = -u_serial[i-1][j+1];
                        break;
                    case B_SE:
                        v_serial[i][j-1] = 0.0;
                        u_serial[i][j]   = 0.0;
                        v_serial[i][j]   = -v_serial[i+1][j];
                        u_serial[i-1][j] = -u_serial[i-1][j-1];
                        break;
                    case B_SW:
                        v_serial[i][j-1] = 0.0;
                        u_serial[i-1][j] = 0.0;
                        v_serial[i][j]   = -v_serial[i-1][j];
                        u_serial[i][j]   = -u_serial[i][j-1];
                        break;
                    case B_NW:
                        v_serial[i][j]   = 0.0;
                        u_serial[i-1][j] = 0.0;
                        v_serial[i][j-1] = -v_serial[i-1][j-1];
                        u_serial[i][j]   = -u_serial[i][j+1];
                        break;
                }
            }
        }
    }

    /* Finally, fix the horizontal velocity at the  western edge to have
     * a continual flow of fluid into the simulation.
     */
    v_serial[0][0] = 2 * vi_serial-v_serial[1][0];
    for (int j = 1; j < jmax_serial+1; j++) {
        u_serial[0][j] = ui_serial;
        v_serial[0][j] = 2 * vi_serial - v_serial[1][j];
    }
}
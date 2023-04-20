#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>

#include "data_serial.h"
#include "vtk_serial.h"
#include "setup_serial.h"
#include "boundary_serial.h"
#include "args_serial.h"
#include "vortex_serial.h"
//! ensure avoidance of data corruption with parallel workflow

void serial_looping(double ** u_return, double ** v_return) {
    allocate_arrays();
    problem_set_up();

    double res;

    /* Main loop */
    int iters = 0;
    double t;
    for (t = 0.0; t < t_end; t += del_t, iters++) {
        if (!fixed_dt)
            set_timestep_interval();

        compute_tentative_velocity();

        compute_rhs();

        res = poisson();

        update_velocity();

        apply_boundary_conditions();

        /*  // commenting out printings and writings
        if ((iters % output_freq == 0)) {
            printf("Step %8d, Time: %14.8e (del_t: %14.8e), Residual: %14.8e\n", iters, t+del_t, del_t, res);
 
            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+del_t);
        }
        */
    } /* End of main loop */

    u_return = u;
    v_return = v;
}

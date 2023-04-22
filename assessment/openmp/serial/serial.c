#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>

#include "serial.h"
#include "data_serial.h"
#include "setup_serial.h"
#include "boundary_serial.h"
#include "vortex_serial.h"
#include "args.h" // for fixed_dt
//! ensure avoidance of data corruption with parallel workflow

void serial_looping() {
    setup_serial();
    allocate_arrays_serial();
    problem_set_up_serial();
    double res_serial;

    /* Main loop */
    int iters_serial = 0;
    double t_serial;
    for (t_serial = 0.0; t_serial < t_end_serial; t_serial += del_t_serial, iters_serial++) {
        if (!fixed_dt)
            set_timestep_interval_serial();

        compute_tentative_velocity_serial();

        compute_rhs_serial();

        res_serial = poisson_serial();

        update_velocity_serial();

        apply_boundary_conditions_serial();

        /*  // commenting out printings and writings
        if ((iters % output_freq == 0)) {
            printf("Step %8d, Time: %14.8e (del_t_serial: %14.8e), Residual: %14.8e\n", iters, t+del_t_serial, del_t_serial, res);
 
            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+del_t_serial);
        }
        */
    } /* End of main loop */

    // * u_size_y_return = u_size_y_serial;
    // * u_size_x_return = u_size_x_serial;
    // * v_size_y_return = v_size_y_serial;
    // * v_size_x_return = v_size_x_serial;
    // u_return = u_serial;
    // v_return = v_serial;

}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <papi.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "data.h"
#include "data_serial.h"
#include "vtk.h"
#include "setup.h"
#include "boundary.h"
#include "args.h"
#include "vortex.h"
#include "serial.h"

struct timespec timer;

void parallel_looping(int * i_return, double * r_return, double * t_return) {
    allocate_arrays();
    problem_set_up();  //! perfect loop (all others in vortex.c)

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

        if ((iters % output_freq == 0)) {
            printf("Step %8d, Time: %14.8e (del_t: %14.8e), Residual: %14.8e\n", iters, t+del_t, del_t, res);

            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+del_t);
        }
    } /* End of main loop */

    *i_return = iters;
    *r_return = res;
    *t_return = t;
}

int validated() {
    #define MAX_ERR 1e-6
    // double ** u_serial, ** v_serial;
    // int u_size_y_serial, u_size_x_serial, v_size_y_serial, v_size_x_serial;
    int index;
    double temp_u, temp_u_serial, temp_v, temp_v_serial;
    serial_looping();
    printf("Validating...\n");
    for (int j = 0; j < u_size_y_serial; j++) {
        for (int i = 0; i < u_size_x_serial; i++) {
            index = i * jmax + j;
            temp_u = u[index];
            temp_u_serial = u_serial[i][j];  // do not replace serial 2d array access
            if (fabs(temp_u - temp_u_serial) > MAX_ERR) {
                printf("Invalidated at:");
                printf("u[i][j] = %lf\n u_serial[i][j] = %lf\n fabs(u[i][j] - u_serial[i][j]) = %lf\n", temp_u, temp_u_serial, fabs(temp_u - temp_u_serial));
                return 0;
            }
        }
    }

    for (int j = 0; j < v_size_y; j++) {
        for (int i = 0; i < v_size_x; i++) {
            index = i * jmax + j;
            temp_v = v[index];
            temp_v_serial = v_serial[i][j];  // do not replace serial 2d array access
            if (fabs(temp_v - temp_v_serial) > MAX_ERR) {
                printf("Invalidated at:");
                printf("v[i][j] = %lf\n v_serial[i][j] = %lf\n fabs(v[i][j] - v_serial[i][j]) = %lf\n", temp_v, temp_v_serial, fabs(temp_v - temp_v_serial));
                return 0;
            }
        }
    }

    printf("Validated!\n");

    return 1;
}

/**
 * @brief Return the current time in seconds.
 */
double get_time() {
	clock_gettime(CLOCK_MONOTONIC, &timer); 
	return (double) (timer.tv_sec + timer.tv_nsec / 1000000000.0);
}

void count_time_parallel_looping(int * i_return, double * r_return, double * t_return) {
    
	double total_problem_time = get_time();
	double allocate_arrays_time = 0.0;
	double problem_set_up_time = 0.0;
	double compute_tentative_velocity_time = 0.0;
	double total_compute_tentative_velocity_time = 0.0;
	double poisson_time = 0.0;
	double total_poisson_time = 0.0;
	double compute_rhs_time = 0.0;
	double total_compute_rhs_time = 0.0;
	double update_velocity_time = 0.0;
	double total_update_velocity_time = 0.0;
    double apply_boundary_conditions_time = 0.0;
    double total_apply_boundary_conditions_time = 0.0;
    double write_time = 0.0;
    double total_write_time = 0.0;

    allocate_arrays_time = get_time();
    allocate_arrays();
    allocate_arrays_time = get_time() - allocate_arrays_time;
    problem_set_up_time = get_time();
    problem_set_up();
    problem_set_up_time = get_time() - problem_set_up_time;

    double res;

    /* Main loop */
    int iters = 0;
    double t;
    for (t = 0.0; t < t_end; t += del_t, iters++) {
        if (!fixed_dt)
            set_timestep_interval();

        compute_tentative_velocity_time = get_time();
        compute_tentative_velocity();
        total_compute_tentative_velocity_time += get_time() - compute_tentative_velocity_time;

        compute_rhs_time = get_time();
        compute_rhs();
        total_compute_rhs_time += get_time() - compute_rhs_time;

        poisson_time = get_time();
        res = poisson();
        total_poisson_time += get_time() - poisson_time;

        update_velocity_time = get_time();
        update_velocity();
        total_update_velocity_time += get_time() - update_velocity_time;

        apply_boundary_conditions_time = get_time();
        apply_boundary_conditions();
        total_apply_boundary_conditions_time += get_time() - apply_boundary_conditions_time;

        if ((iters % output_freq == 0)) {
            printf("Step %8d, Time: %14.8e (del_t: %14.8e), Residual: %14.8e\n", iters, t+del_t, del_t, res);
 
            if ((!no_output) && (enable_checkpoints)) {
                write_time = get_time();
                write_checkpoint(iters, t+del_t);
                total_write_time += get_time() - write_time;
            }
        }
    } /* End of main loop */
    
	total_problem_time = get_time() - total_problem_time;

    printf("=======================================\n");
	printf("\n\nTiming Performance Summary\n");
    printf("=======================================\n");
	printf(" Allocate Arrays Time: %lf\n", allocate_arrays_time);
	printf(" Problem Set Up Time: %lf\n", problem_set_up_time);
	printf(" Average Compute Tentative Velocity Time: %lf\n", total_compute_tentative_velocity_time / iters);
    printf(" Average Poisson Time: %lf\n", total_poisson_time / iters);
    printf(" Average Compute RHS Time: %lf\n", total_compute_rhs_time / iters);
	printf(" Average Update Velocity Time: %lf\n", total_update_velocity_time / iters);
	printf(" Average Apply Boundary Conditions Time: %lf\n", total_apply_boundary_conditions_time / iters);
	printf(" Average Checkpoint Write Time: %lf\n", total_write_time / iters);
	printf("\n Total Problem Time: %lf\n", total_problem_time);
    printf("=======================================\n\n\n");

    *i_return = iters;
    *r_return = res;
    *t_return = t;
}

void count_misses_parallel_looping(int * i_return, double * r_return, double * t_return) {
    #define NUM_PAPI_EVENTS 6

    int events[NUM_PAPI_EVENTS] = { PAPI_L1_DCM, PAPI_L1_ICM, PAPI_L1_TCM, PAPI_L2_DCM, PAPI_L2_ICM, PAPI_L2_TCM };
    long long int remaining_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int allocate_arrays_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int problem_set_up_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int total_compute_tentative_velocity_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int total_poisson_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int total_compute_rhs_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int total_update_velocity_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int total_apply_boundary_conditions_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };
    long long int total_write_events[NUM_PAPI_EVENTS] = { 0, 0, 0, 0, 0, 0 };

    PAPI_start_counters(events, NUM_PAPI_EVENTS);

    PAPI_read_counters(remaining_events, NUM_PAPI_EVENTS);
    allocate_arrays();
    PAPI_read_counters(allocate_arrays_events, NUM_PAPI_EVENTS);

    // PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);
    problem_set_up();
    PAPI_read_counters(problem_set_up_events, NUM_PAPI_EVENTS);

    double res;

    /* Main loop */
    int iters = 0;
    double t;
    for (t = 0.0; t < t_end; t += del_t, iters++) {
        if (!fixed_dt)
            set_timestep_interval();

        PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);
        compute_tentative_velocity();
        PAPI_accum_counters(total_compute_tentative_velocity_events, NUM_PAPI_EVENTS);
        
        // PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);
        compute_rhs();
        PAPI_accum_counters(total_compute_rhs_events, NUM_PAPI_EVENTS);

        // PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);
        res = poisson();
        PAPI_accum_counters(total_poisson_events, NUM_PAPI_EVENTS);
        
        // PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);
        update_velocity();
        PAPI_accum_counters(total_update_velocity_events, NUM_PAPI_EVENTS);

        // PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);
        apply_boundary_conditions();
        PAPI_accum_counters(total_apply_boundary_conditions_events, NUM_PAPI_EVENTS);

        if ((iters % output_freq == 0)) {
            printf("Step %8d, Time: %14.8e (del_t: %14.8e), Residual: %14.8e\n", iters, t+del_t, del_t, res);
 
            if ((!no_output) && (enable_checkpoints)) {
                PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);
                write_checkpoint(iters, t+del_t);
                PAPI_accum_counters(total_write_events, NUM_PAPI_EVENTS);
            }
        }
    } /* End of main loop */

    PAPI_accum_counters(remaining_events, NUM_PAPI_EVENTS);

    printf("=======================================\n");
    printf("Cache Miss Summary");
    printf("=======================================\n");
    printf("\n\nAllocate Array Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", allocate_arrays_events[0], allocate_arrays_events[1], allocate_arrays_events[2]);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", allocate_arrays_events[3], allocate_arrays_events[4], allocate_arrays_events[5]);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", allocate_arrays_events[6], allocate_arrays_events[7], allocate_arrays_events[8]);

    printf("\n\nProblem Set Up Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", problem_set_up_events[0], problem_set_up_events[1], problem_set_up_events[2]);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", problem_set_up_events[3], problem_set_up_events[4], problem_set_up_events[5]);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", problem_set_up_events[6], problem_set_up_events[7], problem_set_up_events[8]);
    
    printf("\n\nAverage Compute Tentative Velocity Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_compute_tentative_velocity_events[0] / iters, total_compute_tentative_velocity_events[1] / iters, total_compute_tentative_velocity_events[2] / iters);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_compute_tentative_velocity_events[3] / iters, total_compute_tentative_velocity_events[4] / iters, total_compute_tentative_velocity_events[5] / iters);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_compute_tentative_velocity_events[6] / iters, total_compute_tentative_velocity_events[7] / iters, total_compute_tentative_velocity_events[8] / iters);
        
    printf("\n\nAverage Compute RHS Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_compute_rhs_events[0] / iters, total_compute_rhs_events[1] / iters, total_compute_rhs_events[2] / iters);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_compute_rhs_events[3] / iters, total_compute_rhs_events[4] / iters, total_compute_rhs_events[5] / iters);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_compute_rhs_events[6] / iters, total_compute_rhs_events[7] / iters, total_compute_rhs_events[8] / iters);
    
    printf("\n\nAverage Update Velocity Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_update_velocity_events[0] / iters, total_update_velocity_events[1] / iters, total_update_velocity_events[2] / iters);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_update_velocity_events[3] / iters, total_update_velocity_events[4] / iters, total_update_velocity_events[5] / iters);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_update_velocity_events[6] / iters, total_update_velocity_events[7] / iters, total_update_velocity_events[8] / iters);

    printf("\n\nAverage Apply Boundary Conditions Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_apply_boundary_conditions_events[0] / iters, total_apply_boundary_conditions_events[1] / iters, total_apply_boundary_conditions_events[2] / iters);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_apply_boundary_conditions_events[3] / iters, total_apply_boundary_conditions_events[4] / iters, total_apply_boundary_conditions_events[5] / iters);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_apply_boundary_conditions_events[6] / iters, total_apply_boundary_conditions_events[7] / iters, total_apply_boundary_conditions_events[8] / iters);
    
    printf("\n\nAverage Checkpoint Write Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_write_events[0], total_write_events[1], total_write_events[2]);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_write_events[3], total_write_events[4], total_write_events[5]);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", total_write_events[6], total_write_events[7], total_write_events[8]);

    printf("\n\nRemaining Looping Events:\n");
    printf("Level 1 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", remaining_events[0], remaining_events[1], remaining_events[2]);
    printf("Level 2 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", remaining_events[3], remaining_events[4], remaining_events[5]);
    // printf("Level 3 Cache Misses:\tData - %lld, Instruction  - %lld, Total %lld\n", remaining_events[6], remaining_events[7], remaining_events[8]);
    printf("=======================================\n\n\n");

    *i_return = iters;
    *r_return = res;
    *t_return = t;
}

/**
 * @brief The main routine that sets up the problem and executes the solving routines routines
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The return value of the application
 */
int main(int argc, char *argv[]) {
    set_defaults();
    parse_args(argc, argv);
    setup();
    // different loopings can be shared up to here

    int iters;
    double res, t;
    if (count_time) {
        count_time_parallel_looping(&iters, &res, &t);
    } else if (count_cache_misses) {
        count_misses_parallel_looping(&iters, &res, &t);
    } else {
        parallel_looping(&iters, &res, &t);
    }

    if (validate && !validated()) return 0;

    if (verbose) print_opts();

    printf("Simulation complete at:\n");
    printf("Step %8d, Time: %14.8e, Residual: %14.8e\n", iters, t, res);
    
    if (!no_output)
        write_result(iters, t);

    free_arrays();

    return 0;
}


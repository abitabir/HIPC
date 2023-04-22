#ifndef VORTEX_SERIAL_H
#define VORTEX_SERIAL_H

void compute_tentative_velocity_serial();
void compute_rhs_serial();
double poisson_serial();
void update_velocity_serial();
void set_timestep_interval_serial();

#endif
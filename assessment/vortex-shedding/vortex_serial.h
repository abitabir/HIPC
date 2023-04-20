#ifndef VORTEX_H
#define VORTEX_H

void compute_tentative_velocity();
void compute_rhs();
double poisson();
void update_velocity();
void set_timestep_interval();

#endif
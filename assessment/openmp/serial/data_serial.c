#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "data_serial.h"
#include "args.h"

double xlength_serial = 4.0;     /* Width of simulated domain */
double ylength_serial = 1.0;     /* Height of simulated domain */
int imax_serial = 512;           /* Number of cells horizontally */
int jmax_serial = 128;           /* Number of cells vertically */

double t_end_serial = 5.0;        /* Simulation runtime */
double del_t_serial = 0.003;      /* Duration of each timestep */
double tau_serial = 0.5;          /* Safety factor for timestep control */

int itermax_serial = 100;         /* Maximum number of iterations in SOR */
double eps_serial = 0.001;        /* Stopping error threshold for SOR */
double omega_serial = 1.7;        /* Relaxation parameter for SOR */
double y_serial = 0.9;            /* Gamma, Upwind differencing factor in PDE discretisation */

double Re_serial = 500.0;         /* Reynolds number */
double ui_serial = 1.0;           /* Initial X velocity */
double vi_serial = 0.0;           /* Initial Y velocity */

double delx_serial, dely_serial;

int fluid_cells_serial = 0;

// Grids used for veclocities, pressure, rhs, flag and temporary f and g arrays
int u_size_x_serial, u_size_y_serial;
double ** u_serial;
int v_size_x_serial, v_size_y_serial;
double ** v_serial;
int p_size_x_serial, p_size_y_serial;
double ** p_serial; 
int rhs_size_x_serial, rhs_size_y_serial;
double ** rhs_serial; 
int f_size_x_serial, f_size_y_serial;
double ** f_serial; 
int g_size_x_serial, g_size_y_serial;
double ** g_serial;
int flag_size_x_serial, flag_size_y_serial;
char ** flag_serial;

/**
 * @brief Allocate a 2D array that is addressable using square brackets
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @return double** A 2D array
 */
double **alloc_2d_array_serial(int m, int n) {
  	double **x;
  	int i;

  	x = (double **)malloc(m*sizeof(double *));
  	x[0] = (double *)calloc(m*n,sizeof(double));
  	for ( i = 1; i < m; i++ )
    	x[i] = &x[0][i*n];
	return x;
}


/**
 * @brief Allocate a 2D char array that is addressable using square brackets
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @return char** A 2D array
 */
char **alloc_2d_char_array_serial(int m, int n) {
  	char **x;
  	int i;

  	x = (char **)malloc(m*sizeof(char *));
  	x[0] = (char *)calloc(m*n,sizeof(char));
  	for ( i = 1; i < m; i++ )
    	x[i] = &x[0][i*n];
	return x;
}

/**
 * @brief Free a 2D array
 * 
 * @param array The 2D array to free
 */
void free_2d_array_serial(void ** array) {
	free(array[0]);
	free(array);
}

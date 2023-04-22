#ifndef DATA_SERIAL_H
#define DATA_SERIAL_H

// macros also defined in parallel
#define C_B      0x0000   // This cell is an obstacle/boundary cell
#define B_N      0x0001   // This obstacle cell has a fluid cell to the north
#define B_S      0x0002   // This obstacle cell has a fluid cell to the south
#define B_W      0x0004   // This obstacle cell has a fluid cell to the west
#define B_E      0x0008   // This obstacle cell has a fluid cell to the east
#define B_NW     (B_N | B_W)
#define B_SW     (B_S | B_W)
#define B_NE     (B_N | B_E)
#define B_SE     (B_S | B_E)
#define B_NSEW   (B_N | B_S | B_E | B_W)

#define C_F      0x0010    // This cell is a fluid cell


extern double xlength_serial;     /* Width of simulated domain */
extern double ylength_serial;     /* Height of simulated domain */
extern int imax_serial;           /* Number of cells horizontally */
extern int jmax_serial;           /* Number of cells vertically */

extern double t_end_serial;       /* Simulation runtime */
extern double del_t_serial;       /* Duration of each timestep */
extern double tau_serial;         /* Safety factor for timestep control */

extern int itermax_serial;        /* Maximum number of iterations in SOR */
extern double eps_serial;         /* Stopping error threshold for SOR */
extern double omega_serial;       /* Relaxation parameter for SOR */
extern double y_serial;           /* Gamma, Upwind differencing factor in PDE */

extern double Re_serial;          /* Reynolds number */
extern double ui_serial;          /* Initial X velocity */
extern double vi_serial;          /* Initial Y velocity */

extern int fluid_cells_serial;

extern double delx_serial, dely_serial;

// Grids used for veclocities, pressure, rhs, flag and temporary f and g arrays
extern int u_size_x_serial, u_size_y_serial;
extern double ** u_serial;
extern int v_size_x_serial, v_size_y_serial;
extern double ** v_serial;
extern int p_size_x_serial, p_size_y_serial;
extern double ** p_serial; 
extern int rhs_size_x_serial, rhs_size_y_serial;
extern double ** rhs_serial; 
extern int f_size_x_serial, f_size_y_serial;
extern double ** f_serial; 
extern int g_size_x_serial, g_size_y_serial;
extern double ** g_serial;
extern int flag_size_x_serial, flag_size_y_serial;
extern char ** flag_serial;

double **alloc_2d_array_serial(int m, int n);
char **alloc_2d_char_array_serial(int m, int n);
void free_2d_array_serial(void ** array);

#endif
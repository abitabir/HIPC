#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>

#include "args.h"
#include "data.h"
#include "data_serial.h"
#include "vtk.h"

int verbose = 0;
int no_output = 0;
int output_freq = 100;
int enable_checkpoints = 0;
int fixed_dt = 0;
int validate = 0;
int count_time = 0;
int count_cache_misses = 0;


static struct option long_options[] = {
	{"del-t",         required_argument, 0, 'd'},
	{"cellx",         required_argument, 0, 'x'},
	{"celly",         required_argument, 0, 'y'},
	{"freq",          required_argument, 0, 'f'},
	{"endtime",       required_argument, 0, 't'},
	{"noio",          no_argument,       0, 'd'},
	{"output",        required_argument, 0, 'o'},
	{"checkpoint",    no_argument,       0, 'c'},	
    {"verbose",       no_argument,       0, 'v'},	
    {"validate",      no_argument,       0, 'i'},	
    {"countseconds",  no_argument,       0, 's'},	
    {"countmisses",   no_argument,       0, 'm'},
    {"help",          no_argument,       0, 'h'},
	{0, 0, 0, 0}
};
#define GETOPTS "d:x:y:f:t:no:cvismh"

/**
 * @brief Print a help message
 * 
 * @param progname The name of the current application
 */
void print_help(char *progname) {
	fprintf(stderr, "A simple computational fluid dynamics solver based around a Karman vortex street.\n\n");
	fprintf(stderr, "Usage: %s [options]\n", progname);
	fprintf(stderr, "Options and arguments:\n");
	fprintf(stderr, "  -x N, --cellx=N         Cells in X-dimension\n");
	fprintf(stderr, "  -y N, --celly=N         Cells in Y-dimension\n");
	fprintf(stderr, "  -t N, --endtime=N       Set the end time (see -n)\n");
    fprintf(stderr, "  -d, --del-t=DELT        Set the simulation timestep size\n");
	fprintf(stderr, "  -f N, --freq=N          Output frequency (i.e. steps between output)\n");
	fprintf(stderr, "  -n, --noio              Disable file I/O\n");
	fprintf(stderr, "  -o FILE, --output=FILE  Set base filename for output (final output will be in BASENAME.vtk\n");
	fprintf(stderr, "  -c, --checkpoint        Enable checkpointing, checkpoints will be in BASENAME-ITERATION.vtk\n");
	fprintf(stderr, "  -v, --verbose           Set verbose output\n");
	fprintf(stderr, "  -i, --validate          Enable verification of the default code against the serial original.\n");
	fprintf(stderr, "  -s, --countseconds      Print timing in seconds of major functions. Will default to only s if used with m\n");
	fprintf(stderr, "  -m, --countmisses       Print data and instruction miss counts during major functions. Will default to only s if used with s\n");
	fprintf(stderr, "  -h, --help              Print this message and exit\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Report bugs to <steven.wright@york.ac.uk>\n");
}

/**
 * @brief Parse the argv arguments passed to the application
 * 
 * @param argc The number of arguments present
 * @param argv An array of the arguments presented
 */
void parse_args(int argc, char *argv[]) {
    int option_index = 0;
    char c;

    while ((c = getopt_long(argc, argv, GETOPTS, long_options, &option_index)) != -1) {
        switch (c) {
			case 'x':
				imax = atoi(optarg);
				imax_serial = imax;
				break;
			case 'y':
				jmax = atoi(optarg);
				imax_serial = jmax;
				break;
			case 't':
                t_end = atof(optarg);
				t_end_serial = t_end;
				break;
            case 'd':
                fixed_dt = 1;
                del_t = atof(optarg);
                del_t_serial = del_t;
                break;
			case 'f':
				output_freq = atoi(optarg);
				break;
			case 'n':
				no_output = 1;
				break;
			case 'o':
				set_basename(optarg);
				break;
			case 'c':
				enable_checkpoints = 1;
				break;
			case 'v':
				verbose = 1;
				break;
			case 'i':
				validate = 1;
				break;
			case 's':
				count_time = 1;
				break;
			case 'm':
				count_cache_misses = 1;
				break;
			case '?':
            case 'h':
				print_help(argv[0]);
				exit(1);
        }
    }
}

/**
 * @brief Print out the current parameters
 * 
 */
void print_opts() {
    printf("=======================================\n");
    printf("Started with the following options\n");
    printf("=======================================\n");
    printf("  del-t            = %14lf\n", del_t);
	printf("  cellx            = %14d\n", imax);
	printf("  celly            = %14d\n", jmax);
	printf("  freq             = %14d\n", output_freq);
	printf("  endtime          = %14.12lf\n", t_end);
	printf("  noio             = %14d\n", no_output);
	printf("  output           = %s\n", get_basename());
	printf("  checkpoint       = %14d\n", enable_checkpoints);	
	printf("  validate         = %14d\n", validate);	
	printf("  counttime        = %14d\n", count_time);	
	printf("  countmisses      = %14d\n", count_cache_misses);	
    printf("=======================================\n\n\n");
}
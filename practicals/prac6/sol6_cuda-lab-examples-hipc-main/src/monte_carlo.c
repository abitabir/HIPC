#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char* argv[])
{
    long niter = 1000000L;
    double x,y;
    int i;
    int count = 0;
    double z;
    double pi;
	struct timespec tstart={0,0}, tend={0,0};
    srand(time(NULL));

	clock_gettime(CLOCK_MONOTONIC, &tstart);

    //main loop
    for (i = 0; i < niter; ++i)
    {
        //get random points
        x = (double)random() / RAND_MAX;
        y = (double)random() / RAND_MAX;
        z = sqrt((x*x)+(y*y));
        //check to see if point is in unit circle
        if (z <= 1)
        {
            ++count;
        }
    }
    //p = 4(m/n)
    pi = ((double)count/(double)niter)*4.0;

	clock_gettime(CLOCK_MONOTONIC, &tend);

    printf("Pi: %f\n", pi);
    printf("Time takes %.6f ms\n",
           ((double)tend.tv_sec * 0.001 + 1.0e-6*tend.tv_nsec) - 
           ((double)tstart.tv_sec * 0.001 + 1.0e-6*tstart.tv_nsec));

	return 0;
}

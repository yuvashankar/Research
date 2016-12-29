#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary constants.
    double *data, *result, *period, *scales, *frequency;

    int n = DATA_SIZE;
    const int J = MAX_I - MIN_I;

    //Memory Allocations
    data =           (double*) malloc(n *     sizeof(double));
    result =         (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL);

    //Get Scales and Frequencies
    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY);
    frequency = IdentifyFrequencies(scales, J);

    //Populate the data array
    TestCases(data, 5);

    //Compute the ERSP
    ERSP (data, scales, FS, n, J, 77, 
    result);
    
    // Write to file
    WriteFile(result, frequency, J, n, "DATA.log");
    // Plot(result, period,  n, J);

    //Free up Memory
    free(data);  free(result);
    free(scales); free(frequency);
    t = omp_get_wtime() - t;
    printf("Execution Time: %f\n", t);
    return 0;
}

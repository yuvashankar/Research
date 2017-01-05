#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary arrays.
    double *data, *data_2D, *result, *scales, *frequency;

    int n = DATA_SIZE;
    const int J = (int) MAX_I - MIN_I;

    const int trials = 77;

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    data_2D =  (double*) malloc(n * J * sizeof(double));
    result  =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);

    TestCases(data, 2);

    // //Populate the data array
    // for (int i = 0; i < trials; ++i)
    // {
    //     TestCases(data, 2);
    //     for (int j = 0; j < n; ++j)
    //     {
    //         data_2D[i * n + j] = data[j];
    //     }
    // }

    
    Wavelet(data, scales, 
            FS, n, J,
            result);
    
    // Compute the ERSP
    // ERSP (data_2D, scales, FS, n, J, trials, PAD_FLAG, 
    // result);

    // Write to file
    WriteFile(result, frequency, J, n, "DATA.log");
    // Plot(result, frequency,  n, J);

    //Free up Memory
    free(data_2D);
    free(data);  free(result);
    free(scales); free(frequency);
    
    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}
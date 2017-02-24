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

    double *con_result, *wavelet_result; 

    int n = DATA_SIZE;
    const int J = (int) MAX_I - MIN_I;

    const int trials = 77;

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    data_2D =  (double*) malloc(n * trials * sizeof(double));

    con_result =     (double*) malloc(n * J * sizeof(double));
    wavelet_result = (double*) malloc(n * J * sizeof(double));

    result  =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);

    //Populate the data array
    for (int i = 0; i < trials; ++i)
    {
        TestCases(data, 5);
        for (int j = 0; j < n; ++j)
        {
            data_2D[i * n + j] = data[j];
        }
    }

    // Compute the ERSP
    ERSP (data_2D, scales, FS, n, J, trials, PAD_FLAG, 
    result);


    WriteFile(result, frequency, J, n, "DATA.log");

    //Free up Memory
    free(data_2D);
    free(data);  free(result);
    free(scales); free(frequency);

    free(con_result); free(wavelet_result);
    
    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}
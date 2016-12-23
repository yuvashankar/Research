#include "wavelet.h"
#include "math.h"
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    double *data, *result, *wavelet_result, *baseline_out, *period, *scales;

    // double t = omp_get_wtime();
    //Initialize the necessary constants.
    int n = DATA_SIZE;
    const int J = MAX_I - MIN_I;

    //Memory Allocations
    data =           (double*) malloc(n *     sizeof(double));
    result =         (double*) malloc(n * J * sizeof(double));
    wavelet_result = (double*) malloc(n * J * sizeof(double));
    baseline_out =   (double*) malloc(n * J * sizeof(double));
    period =         (double*) malloc(    J * sizeof(double));
    
    assert(data != NULL); assert(result != NULL); assert(period != NULL);
    assert(wavelet_result != NULL); assert(baseline_out != NULL);

    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY);

    //populate the data array
    TestCases(data, 1);
    // n = ReadFile(data, "sst_nino3.dat");
    


    CleanData(data, DATA_SIZE);

    
    Wavelet(data, period, scales, 
        FS, n, J,
        wavelet_result);
    
    
    // Write to file
    WriteFile(wavelet_result, period, J, n, "DATA.log");

    free(data); free(result); free(period); free(wavelet_result); free(baseline_out);
    free(scales);
    // t = omp_get_wtime() - t;
    // printf("Execution Time: %f\n", t);
    return 0;
}

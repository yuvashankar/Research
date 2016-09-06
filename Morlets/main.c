#include "wavelet.h"
#include "math.h"
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    double *data, *result, *wavelet_result, *baseline_out, *period;

    //Size of Data
    int n = DATA_SIZE;
    int J; 
    double dj, dt, s0;

    dj = 0.125;

    dt = 1.0/FS;
    s0 = 2 * dt;

    J = (int) ceil(log2 ( 1.0/(s0 * MIN_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) )/dj);
    printf("dt = %f, dj = %f, s0 = %f, J = %d\n", dt, dj, s0, J);

    //Memory Allocations
    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    wavelet_result = malloc( J * n * sizeof(double));
    baseline_out = malloc (J * n * sizeof(double));
    period = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(period != NULL);
    assert(wavelet_result != NULL); assert(baseline_out != NULL);

    //populate the data array
    TestCases(data, 2);

    // printf("Data Size: %d\n", DATA_SIZE);
    CleanData(data, DATA_SIZE);
    int start = (int) floor( log2( 1.0/(s0 * MAX_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) ) /dj);
    int bic = 77;

    for (int i = 0; i < bic; ++i)
    {
        // if (i == 0)
        // {
            //Compute wavelet analysis
            Wavelet(data, period,
            FS, n, dj, s0, J, MAX_FREQUENCY,
            wavelet_result);
        // }
        
        

        RemoveBaseline(wavelet_result, DATA_SIZE, J, 1, FS, baseline_out);

        for (int j = start; j < J; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                result[j * n + k] += baseline_out[j * n + k];
            }
        }
    }

    for (int i = start; i < J; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i * n + j] = result[i * n + j] / bic;
        }
    }

    //Write to file
    WriteFile(wavelet_result, period, J, n, "DATA.log");

    
    //sanitation engineering
    free(data); free(result); free(period); free(wavelet_result); free(baseline_out);
    return 0;
}

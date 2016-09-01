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

    //populate the data array
    TestCases(data, 1);

    // printf("Data Size: %d\n", DATA_SIZE);
    CleanData(data, DATA_SIZE);

    for (int i = 0; i < 77; ++i)
    {
        //Compute wavelet analysis
        Wavelet(data, period,
            FS, n, dj, s0, J, MAX_FREQUENCY,
            wavelet_result);

        RemoveBaseline(wavelet_result, DATA_SIZE, J, 1, FS, baseline_out);

        for (int j = 0; j < J * n ; ++j)
        {
            result[j] += baseline_out[j];
        }
    }

    for (int i = 0; i < J * n; ++i)
    {
        result[i] = result[i]/77;
    }

    //Write to file
    WriteFile(result, period, J, n, "DATA.log");

    
    //sanitation engineering
    free(data); free(result); free(period); free(wavelet_result); free(baseline_out);
    return 0;
}

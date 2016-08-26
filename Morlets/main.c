#include "wavelet.h"
#include "math.h"
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    double *data, *result, *period;

    //Size of Data
    int n = DATA_SIZE;
    int J; 
    double dj, dt, s0;

    dj = 0.0625;

    dt = 1.0/FS;
    s0 = 2 * dt;

    J = (int) ceil(log2 ( 1.0/(s0 * MIN_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) )/dj);
    printf("dt = %f, dj = %f, s0 = %f, J = %d\n", dt, dj, s0, J);

    //Memory Allocations
    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    period = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(period != NULL);

    //populate the data array
    // int DATA_SIZE = ReadFile(data, "Sample_EEG_data.dat");
    TestCases(data, 3);
    // printf("Data Size = %d\n", DATA_SIZE);

    double mean = gsl_stats_mean(data, 1, DATA_SIZE);
    double sDeviation = gsl_stats_sd_m(data, 1, DATA_SIZE, mean);
    // printf("Mean = %f, Standard Deviation = %f\n", mean, sDeviation);

    //Compute the Z-Score or Standard Score
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        data[i] = (data[i] - mean)/sDeviation;
    }

    //Compute wavelet analysis
    double execution_time = omp_get_wtime();
    Wavelet(data, period ,
        FS, n, dj, s0, J, MAX_FREQUENCY,
        result);
    execution_time = omp_get_wtime() - execution_time;
    printf("Execution Time: %f\n", execution_time);

    //Write to file
    WriteFile(result, period, J, n, "DATA.log");

    
    //sanitation engineering
    free(data); free(result); free(period);
    return 0;
}

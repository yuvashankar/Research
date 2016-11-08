#include "wavelet.h"
#include "math.h"
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    double *data, *result, *wavelet_result, *baseline_out, *period;

    //Initialize the necessary constants.
    const int n = DATA_SIZE;
    const int J = (log2(DATA_SIZE) + 1)/D_J;
    const double dt = 1.0/FS;
    const double s0 = 2 * dt;

    //Memory Allocations
    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    wavelet_result = malloc( J * n * sizeof(double));
    baseline_out = malloc (J * n * sizeof(double));
    period = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(period != NULL);
    assert(wavelet_result != NULL); assert(baseline_out != NULL);

    //populate the data array
    TestCases(data, 6);

    Wavelet(data, period, FS, n, s0, J, MAX_FREQUENCY, wavelet_result);

    //Write to file
    WriteFile(wavelet_result, period, J, n, "DATA.log");

    //sanitation engineering
    free(data); free(result); free(period); free(wavelet_result); free(baseline_out);
    return 0;
}
#include "wavelet.h"
#include "math.h"
#include <omp.h>

int main(void)
{
    double *data, *result, *frequency;

    //Size of Data
    int n = DATA_SIZE;
    double dj, dt, s0, J;
    dj = 0.25;

    dt = 1.0/FS;
    s0 = 2 * dt;
    J = ceil( log2( (W_0 * MAX_FREQUENCY)/(8 * M_PI * s0) )/dj);

    printf("dt = %f, dj = %f, s0 = %f, J = %f\n", dt, dj, s0, J);

    //Memory Allocations
    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    frequency = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(frequency != NULL);

    //populate the data array
    TestCases(data, 1);
    
    //Compute wavelet analysis
    double execution_time = omp_get_wtime();
    int out  = Wavelet(data, frequency ,
        FS, n, dj, s0, J, MIN_FREQUENCY,
        result);
    execution_time = omp_get_wtime() - execution_time;
    printf("Execution Time: %f\n", execution_time);

    //Write to file
    int writeFlag = WriteFile(result, frequency, J, n, "DATA.log");

    
    //sanitation engineering
    free(data); free(result); free(frequency);
    return 0;
}
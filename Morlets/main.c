#include "Morlet.h"
#include "math.h"
#include <omp.h>

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double *data, *result, *frequency;

    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    double dj, dt, s0, J;
    dt = 1.0/FS;
    dj = 0.25;
    s0 = 2 * dt;
    // J = (log2(n * dt)/s0)/dj;
    J = ceil( log2( (W_0 * MAX_FREQUENCY)/(8 * M_PI * s0) )/dj);

    // printf("dt = %f, dj = %f, s0 = %f, J = %f\n", dt, dj, s0, J);

    //Memory Allocation
    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    frequency = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(frequency != NULL);

    // populate the data array
    TestCases(data, 2);
    
    //Compute the CWT
    double execution_time = omp_get_wtime();
    int out  = Wavelet(data, FS, n, dj, s0, J, result, frequency);
    execution_time = omp_get_wtime() - execution_time;
    printf("Execution Time: %f\n", execution_time);

    //Write to File
    int writeFlag = WriteFile(result, frequency, J, n, "DATA.log");

    
    //Sanitation Engineering
    free(data); free(result); free(frequency);
    fclose(out_file);
    return 0;
}
#include "Morlet.h"
#include "math.h"
#include <omp.h>

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double *data, *result, *frequency;

    double dj, dt, s0, maxScale;
    double parallel_time, execution_time;
    int J;

    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);
    
    execution_time = omp_get_wtime();

    dt = 1.0/FS;
    dj = 0.25;
    s0 = 2 * dt;

    maxScale = MAX_FREQUENCY / FOURIER_WAVELENGTH_FACTOR;
    J = ceil(log2 (maxScale/s0)/dj);

    printf("dt = %f, dj = %f, s0 = %f, J = %f, Max Scale = %f\n", dt, dj, s0, J, maxScale);

    data = malloc(n * sizeof(double));
    result = malloc( (int)J * n * sizeof(double));
    frequency = malloc( (int)J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(frequency != NULL);

    // populate the data array
    FillData(data);
    
    parallel_time = omp_get_wtime();
    int out  = Wavelet(data, dt, n, dj, s0, J, result, frequency);
    parallel_time = omp_get_wtime() - parallel_time;

    int writeFlag = WriteFile(result, frequency, J, n, "DATA.log");
    printf("Finished Writing to file\n");
    execution_time = omp_get_wtime() - execution_time;
    printf("Execution Time: %f, Parallel Time: %f, Serial Time: %f\n", execution_time, parallel_time, (execution_time - parallel_time));

    free(data); free(result); free(frequency);
    fclose(out_file);
    return 0;
}
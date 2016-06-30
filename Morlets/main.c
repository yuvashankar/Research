#include "Morlet.h"
#include "math.h"
#include <omp.h>

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double *data, *result, *frequency;

<<<<<<< HEAD
=======
    double dj, dt, s0, maxScale;
    double parallel_time, execution_time;
    int J;

>>>>>>> parallel
    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);
    
    execution_time = omp_get_wtime();

    dt = 1.0/FS;
    dj = 0.25;
    s0 = 2 * dt;
<<<<<<< HEAD
    // J = (log2(n * dt)/s0)/dj;
    J = ceil( log2( (W_0 * MAX_FREQUENCY)/(8 * M_PI * s0) )/dj);
=======

    maxScale = MAX_FREQUENCY / FOURIER_WAVELENGTH_FACTOR;
    J = ceil(log2 (maxScale/s0)/dj);
>>>>>>> parallel

    printf("dt = %f, dj = %f, s0 = %f, J = %f, Max Scale = %f\n", dt, dj, s0, J, maxScale);

    data = malloc(n * sizeof(double));
    result = malloc( (int)J * n * sizeof(double));
    frequency = malloc( (int)J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(frequency != NULL);

    // populate the data array
<<<<<<< HEAD
    TestCases(data, 2);
=======
    // FillData(data);
    TestCases(data, 1);
>>>>>>> parent of f037623... "all of the tests came out positive i'm pretty confident that I have the right graphs and the right math now going to work on parallelizing it"
    
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
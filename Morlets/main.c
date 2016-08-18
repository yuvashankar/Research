#include "wavelet.h"
#include <math.h>
#include <omp.h>
#include <mpi.h>
//#include <mpi/mpi.h>
#include <gsl/gsl_statistics.h>
<<<<<<< HEAD


int main(void)
{
    FILE* out_file = fopen("data.log", "w");
    double *data; 


    // //Size of Data
    int n = DATA_SIZE;
    // int J; 
    double dj, dt, s0;

    // dj = 0.0625;

    dt = 1.0/FS;
    // s0 = 2 * dt;
=======
#include <cuda.h>

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
>>>>>>> parent of 295781a... I'm making some figures for my thesis and I'm going to be messing around with Morlets

    J = (int) ceil(log2 ( 1.0/(s0 * MIN_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) )/dj);
    printf("dt = %f, dj = %f, s0 = %f, J = %d\n", dt, dj, s0, J);

<<<<<<< HEAD
    // //Memory Allocations
    data = malloc(n * sizeof(double));

    data_fft = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
    // result = malloc(J * n * sizeof(double));
    // period = malloc(J * sizeof(double));
    // assert(data != NULL); assert(result != NULL); assert(period != NULL);

    // //populate the data array
    TestCases(data, 4);
    // int dataSize = ReadFile(data, "Main_debug.dat");
    // printf("Data Size = %d\n", dataSize);
=======
    //Memory Allocations
    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    period = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(period != NULL);

    //populate the data array
    // TestCases(data, 2);
    int dataSize = ReadFile(data, "Main_debug.dat");
    printf("Data Size = %d\n", dataSize);
>>>>>>> parent of 295781a... I'm making some figures for my thesis and I'm going to be messing around with Morlets

    double mean = gsl_stats_mean(data, 1, dataSize);
    double sDeviation = gsl_stats_sd_m(data, 1, dataSize, mean);
    // printf("Mean = %f, Standard Deviation = %f\n", mean, sDeviation);

    //Compute the Z-Score or Standard Score
    for (int i = 0; i < dataSize; ++i)
    {
        data[i] = (data[i] - mean)/sDeviation;
    }

    for (int i = 0; i < DATA_SIZE; ++i)
    {
        fprintf(out_file, "%d\t%f\n", i, data[i]);
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

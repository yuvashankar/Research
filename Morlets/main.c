#include "wavelet.h"
#include <math.h>
#include <omp.h>
#include <gsl/gsl_statistics.h>


int main(void)
{
    FILE* out_file = fopen("data.log", "w");
    double *data, 
    *period,
    *wavelet_result; 

    double dj, dt, s0;

    int J;



    fftw_complex *data_in, *data_fft;

    fftw_plan plan_forward;

    // //Size of Data
    int n = DATA_SIZE;
    
    //Begin Wavelet Analysis
    dj = 0.0625;
    dt = 1.0/FS;
    s0 = 2 * dt;
    J = (int) ceil(log2 ( 1.0/(s0 * MIN_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) )/dj);
    // //Memory Allocations
    data = malloc( n * sizeof(double) );
    wavelet_result = malloc( J * n * sizeof(double) );
    period = malloc(J * sizeof(double));
    
    
    //populate the data array
    TestCases(data, 5);
    // CleanData(data, DATA_SIZE);

    Wavelet(data, period,
            FS, DATA_SIZE, dj, s0, J, MAX_FREQUENCY,
            wavelet_result);

    WriteFile(wavelet_result, period, J, FS, "DATA.log");


    fclose(out_file);

    fftw_free(data_fft);
    fftw_free(data_in);
    fftw_destroy_plan(plan_forward);

    free(data);
    free(wavelet_result);


    return 0;
}

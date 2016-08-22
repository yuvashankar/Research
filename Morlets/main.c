 #include "wavelet.h"
#include <math.h>
#include <omp.h>
#include <gsl/gsl_statistics.h>


int main(void)
{
    FILE* out_file = fopen("data.log", "w");
    double *data; 

    fftw_complex *data_in, *data_fft;

    fftw_plan plan_forward;

    // //Size of Data
    int n = DATA_SIZE;
    // int J; 

    // //Memory Allocations
    data = malloc(n * sizeof(double));
    data_in =  (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * n );
    data_fft = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * n );

    plan_forward = fftw_plan_dft_1d(n, data_in, data_fft, 
        FFTW_FORWARD, FFTW_ESTIMATE);
    
    //populate the data array
    TestCases(data, 5);


    for (int i = 0; i < n; ++i)
    {
        data_in[i][0] = data[i];
    }

    fftw_execute(plan_forward);

    for (int i = 0; i < n; ++i)
    {
        fprintf(out_file, "%f\t%f\t%f\n", i/FS, data_in[i][0], data_in[i][1]);
    }

    fclose(out_file);

    fftw_free(data_fft);
    fftw_free(data_in);
    fftw_destroy_plan(plan_forward);


    return 0;
}

#include "Morlet.h"

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double * data, *result;

    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    double dj, dt, s0, J;
    dt = 1.0/FS;
    dj = 0.25;
    s0 = 2 * dt;
    J = 7/dj;
    // J = 10;

    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    assert(data != NULL); assert(result != NULL);

    //populate the data array
    fillData(data);
    

    // fftw_complex* fourier_morlet = fftw_alloc_complex(n);
    // fftw_complex* output = fftw_alloc_complex(n);
    
    // const double df= 1./n/dt;
    // double sign=1.;
    // for (int i = 0; i < n/2; ++i)
    // {
    //     fourier_morlet[i][0] = sign * NewFourierMorlet(i*df, 5.0, 0.5, n);
    //     fourier_morlet[i][1] = 0.0;
        
    //     fourier_morlet[n - i - 1][0] = 0.0;
    //     fourier_morlet[n - i - 1][1] = 0.0;
    //     sign *= -1.0;
    // }

    // fftw_plan plan_backward = fftw_plan_dft_1d(n, fourier_morlet, output, FFTW_FORWARD, FFTW_ESTIMATE);
    // fftw_execute(plan_backward);

    // for (int i = 0; i < n; ++i)
    // {
    //     fprintf(out_file, "%d\t%f\t%f\n", i, fourier_morlet[i][0], output[i][0]);
    // }
    

    int out  = Wavelet(data, dt, n, dj, s0, J, result);

    int writeFlag = WriteFile(result, J, n, "DATA.log");

    free(data); free(result);
    // fftw_free(fourier_morlet); fftw_free(output);
    // fftw_destroy_plan(plan_backward);

    fclose(out_file);

    return 0;
}
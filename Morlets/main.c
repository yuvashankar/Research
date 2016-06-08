#include "Morlet.h"

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double *data, *result;

    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    double dj, dt, s0, J;
    dt = 1.0/FS;
    dj = 0.25;
    s0 = 2 * dt;
    J = 7/dj; //Where did this 7 come from? dunno m8.
    // printf("dt = %f, dj = %f, s0 = %f, J = %f\n", dt, dj, s0, J);

    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    assert(data != NULL); assert(result != NULL);

    fftw_complex *data_in, *fft_data, *fftw_result;
    fftw_plan plan_backward;

    data_in = fftw_alloc_complex(n);
    fft_data = fftw_alloc_complex(n);
    fftw_result = fftw_alloc_complex(n);

    // populate the data array
    fillData(data);

    double df = 1.0/n/dt;
    printf("df = %f, fs = %f, n = %d\n", df, FS, n);

    for (int i = 0; i < n; ++i)
    {
        fft_data[i][0] = 0.0;
        fft_data[i][1] = 0.0;
    }

    // double sign = 1.0;
    // for (int i = 0; i < FS/2; ++i)
    // {
    //     fft_data[i][0] = sign * NewFourierMorlet(i*df, 5.0, 0.5, n);
    //     fft_data[i][1] = 0.0;

    //     sign *= -1.0;
    // }
    
    // plan_backward = fftw_plan_dft_1d(n, fft_data, fftw_result, FFTW_FORWARD, FFTW_ESTIMATE);
    // fftw_execute(plan_backward);

    // for (int i = 0; i < n; ++i)
    // {
    //     double tmp = Magnitude(fftw_result[i][0], fftw_result[i][1]);
    //     fprintf(out_file, "%d\t%f\t%f\t%f\t%f\n", i, fftw_result[i][0], fftw_result[i][1], tmp, fft_data[i][0]);
    // }

    int out  = Wavelet(data, dt, n, dj, s0, J, result);
    int writeFlag = WriteFile(result, J, n, "DATA.log");

    free(data); free(result);
    fclose(out_file);
    
    fftw_free(data_in); fftw_free(fft_data); fftw_free(fftw_result);
    fftw_destroy_plan(plan_backward);

    return 0;
}
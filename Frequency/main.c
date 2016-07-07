#include "wavelet.h"
#include "math.h"
#include <omp.h>

int main(void)
{
    fftw_complex *data_in, *fftw_result;
    fftw_plan plan_backward;

    FILE* out_file=fopen("DATA.log","w");
    if (out_file == NULL) return -1;

    //Size of Data
    int n = DATA_SIZE;
    
    const int pad = floor(log(n)/log(2.0) + 0.499);
    const double PADDED_SIZE = pow(2, pad + 1);

    const double dt = 1.0/FS;
    const double df = FS/PADDED_SIZE;

    const double scale = 1.0;
    double frequency = scale * (2 * M_PI)/(W_0);
    double period = 1.0/frequency;

    printf("Scale is: %f\n", scale);
    printf("Corrosponding Frequency is: %f\n", frequency);
    printf("Period is: %f\n", period);
    printf("Discrete Period is: %f\n", period * FS);


    //Memory Allocations
    data_in =            fftw_alloc_complex(PADDED_SIZE);
    double *data = malloc(PADDED_SIZE * sizeof(double)) ;
    
    fftw_result =        fftw_alloc_complex(PADDED_SIZE);
    assert(data_in != NULL); assert(fftw_result != NULL);

    const double k = exp(-0.5 * W_0_2);
    const double cSigma = pow(1.0 + exp(-W_0_2) - 2*exp(-0.75*W_0_2), -0.5);
    double normal = sqrt(2 * M_PI * scale * FS);

    for (int i = 0; i < PADDED_SIZE; ++i)
    {
        data_in[i][0] = FourierMorlet(i*df, scale, k, cSigma, normal);
        data_in[i][1] = 0.0;
    }

    for (int i = 0; i < PADDED_SIZE; ++i)
    {
        double value = Magnitude(fftw_result[i][0], fftw_result[i][1]);

        fprintf(out_file, "%d\t%f\t%f\t%f\n", i, fftw_result[i][0], fftw_result[i][1], 
            value);
    }

    plan_backward = fftw_plan_dft_1d(PADDED_SIZE, fftw_result, data_in, 
            FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan_backward);

    

    //sanitation engineering
    free(data);

    fftw_free(data_in); fftw_free(fftw_result);
    fftw_destroy_plan(plan_backward);

    fclose(out_file);

    return 0;
}
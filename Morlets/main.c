

#include "Morlet.h"

int main(void)
{
    //Allocate Memory for the necessary arrays.
    double *data = malloc(DATA_SIZE * sizeof(double));
    assert(data != NULL);
    
    double *result = malloc(DATA_SIZE * sizeof(double));
    // double *complexResult = malloc (DATA_SIZE*MAX_SCALES * sizeof(double));
    // assert(complexResult != NULL);
    assert(result != NULL); 


    double *conWindow = malloc(DATA_SIZE * sizeof(double));
    // double *complexWindow = malloc(MAX_CONV_SIZE * DATA_SIZE * sizeof(double));
    // assert(complexWindow != NULL); 
    assert(conWindow != NULL); 
    
    //FFTW allocations
    fftw_complex *data_in, *fft_result;
    fftw_plan plan_forward;

    data_in = fftw_alloc_complex(DATA_SIZE); fft_result = fftw_alloc_complex(DATA_SIZE);

    //Open up the Output File
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    //Populate the filter data and input data. 
    FillDataComplex(data_in);
    CreateComplexFilter(conWindow, FREQ);

    //Plan and calculate the Fourier Transform of the data. 
    plan_forward = fftw_plan_dft_1d(DATA_SIZE, data_in, fft_result, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    double value;
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        value = fft_result[i][0] * conWindow[i];
        result[i] = value;
    }

    for (int i = 0; i < DATA_SIZE; ++i)
    {
        fprintf(out_file, "%d\t%f\t%f\t%f\t%f\n", i, data_in[i][0], fft_result[i][0], fft_result[i][1], result[i]);
    }

    fclose(out_file);

    //Sanitation Engineering
    free(data);
    free(result);
    free(complexResult);
    free(conWindow);
    free(complexWindow);

    //FFTW sanitation. 
    fftw_destroy_plan(plan_forward);
    fftw_free(data_in); fftw_free(fft_result);

    return 0;
}
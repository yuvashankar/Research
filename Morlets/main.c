

#include "Morlet.h"

int main(void)
{
    //Allocate Memory for the necessary arrays.
    double *conWindow = malloc(DATA_SIZE * sizeof(double));
    assert(conWindow != NULL); 
    
    //FFTW allocations
    fftw_complex *data_in, *fft_data, *result;
    fftw_plan plan_forward, plan_backwards;

    data_in = fftw_alloc_complex(DATA_SIZE); 
    fft_data = fftw_alloc_complex(DATA_SIZE);
    result = fftw_alloc_complex(DATA_SIZE);
    

    //Open up the Output File
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    //Populate the filter data and input data. 
    FillDataComplex(data_in);
    CreateComplexFilter(conWindow, FREQ);

    //Plan and calculate the Fourier Transform of the data, the fft is stored in fft_data
    plan_forward = fftw_plan_dft_1d(DATA_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    double value;
    double value2;
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        value = fft_data[i][0] * conWindow[i];
        // value2 = fft_data[i][1] * conWindow[i];

        fft_data[i][0] = value;
        // fft_data[i][1] = value2;
    }
    
    //Calculate the backwards FFT of the result and daughter wavelets.
    plan_backwards = fftw_plan_dft_1d(DATA_SIZE, fft_data, result, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_backwards);

    //Print to file
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        fprintf(out_file, "%d\t%f\t%f\t%f\t%f\n", i, data_in[i][0], fft_data[i][0], fft_data[i][1], result[i][0]);
    }

    fclose(out_file);

    //Sanitation Engineering
    free(conWindow);

    //FFTW sanitation. 
    fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backwards);
    fftw_free(data_in); fftw_free(fft_data); fftw_free(result);

    return 0;
}
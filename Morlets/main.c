

#include "Morlet.h"

int main(void)
{
    //Find the closest 2 power of DATA_SIZE, and pad the arrays to the left with that. 
    int pad = floor(log(DATA_SIZE)/log(2.0) + 0.4999);
    double PADDED_SIZE = pow(2, pad + 1);

    //Allocate Memory for the necessary arrays.
    double *conWindow = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
    double *output = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
    double *complexOutput = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
    assert(conWindow != NULL); assert (output != NULL); assert(complexOutput != NULL);
    
    // //FFTW allocations
    // fftw_complex *data_in, *fft_data, *result;
    // fftw_plan plan_forward, plan_backwards;

    // data_in = fftw_alloc_complex(PADDED_SIZE); 
    // fft_data = fftw_alloc_complex(PADDED_SIZE);
    // result = fftw_alloc_complex(PADDED_SIZE);
    

    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    //Populate the filter data and input data. 
    // FillDataComplex(data_in);
    CreateComplexFilter(conWindow, FREQ);

    for (int i = 0; i < DATA_SIZE; ++i)
    {
        fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
            conWindow[0*DATA_SIZE + i], conWindow[1*DATA_SIZE + i], conWindow[2*DATA_SIZE + i],
            conWindow[3*DATA_SIZE + i], conWindow[4*DATA_SIZE + i], conWindow[5*DATA_SIZE + i],
            conWindow[6*DATA_SIZE + i], conWindow[7*DATA_SIZE + i], conWindow[8*DATA_SIZE + i],
            conWindow[9*DATA_SIZE + i]);
    }

    // //Plan and calculate the Fourier Transform of the data, the fft is stored in fft_data
    // plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
    // fftw_execute(plan_forward);

    // double value;
    // for (int i = 0; i < DATA_SIZE; ++i)
    // {
    //     value = fft_data[i][0] * conWindow[i];
    //     // value2 = fft_data[i][1] * conWindow[i];

    //     fft_data[i][0] = value;
    //     // fft_data[i][1] = value2;
    // }
    
    // //Calculate the backwards FFT of the result and daughter wavelets.
    // plan_backwards = fftw_plan_dft_1d(PADDED_SIZE, fft_data, result, FFTW_BACKWARD, FFTW_ESTIMATE);
    // fftw_execute(plan_backwards);

    // //Print to file
    // for (int i = 0; i < DATA_SIZE; ++i)
    // {
    //     value = Magnitude(result[i][0], result[i][1]);
    //     fprintf(out_file, "%d\t%f\t%f\t%f\t%f\n", i, data_in[i][0], result[i][0]/1500, result[i][1]/1500, value/1500);
    // }

    fclose(out_file);

    //Sanitation Engineering
    free(conWindow);
    free(output);
    free(complexOutput);

    // //FFTW sanitation. 
    // fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backwards);
    // fftw_free(data_in); fftw_free(fft_data); fftw_free(result);

    return 0;
}
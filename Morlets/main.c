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

    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    assert(data != NULL); assert(result != NULL);

    //populate the data array
    fillData(data);
    

    int out  = Wavelet(data, dt, n, dj, s0, J, result);
    printf("Wavelet Flag = %d\n", out);
    for (int i = 0; i < J; ++i)
    {
        for (int j = 0; j < DATA_SIZE; ++j)
        {
            // value = Magnitude(result[i*n + j], result[i*n + j]);
            fprintf(out_file, "%f\t", result[i*n + j]);
        }

        fprintf(out_file, "\n");
    }

    free(data); free(result);
    fclose(out_file);

    // //Find the closest 2 power of DATA_SIZE, and pad the arrays to the left with that. 
    // int pad = floor(log(DATA_SIZE)/log(2.0) + 0.4999);
    // double PADDED_SIZE = pow(2, pad + 1);

    // //Allocate Memory for the necessary arrays.
    // double *conWindow = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
    // double *output = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
    // double *complexOutput = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
    // assert(conWindow != NULL); assert (output != NULL); assert(complexOutput != NULL);
    
    // //FFTW allocations
    // fftw_complex *data_in, *fft_data, *result;
    // fftw_plan plan_forward, plan_backwards;

    // data_in = fftw_alloc_complex(PADDED_SIZE); 
    // fft_data = fftw_alloc_complex(PADDED_SIZE);
    // result = fftw_alloc_complex(PADDED_SIZE);
    

    // //Open the Output file
    // FILE* out_file=fopen("DATA.log","w");
    // assert(out_file != NULL);

    // //Populate the filter data and input data. 
    // int n = FillDataComplex(data_in);

    // int scales = CreateComplexFilter(conWindow);




    // //Plan and calculate the Fourier Transform of the data, the fft is stored in fft_data
    // plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
    // fftw_execute(plan_forward);

    // double value;
    // for (int i = 0; i < scales; ++i)
    // {
    //     for (int j = 0; j < DATA_SIZE; ++j)
    //     {
    //         value = fft_data[j][0] * conWindow[i * DATA_SIZE + j];
    //         // value2 = fft_data[i][1] * conWindow[i];

    //         fft_data[j][0] = value;
    //         // fft_data[i][1] = value2;
    //     }
        
    //     //Calculate the backwards FFT of the result and daughter wavelets.
    //     plan_backwards = fftw_plan_dft_1d(PADDED_SIZE, fft_data, result, FFTW_BACKWARD, FFTW_ESTIMATE);
    //     fftw_execute(plan_backwards);

    //     for (int j = 0; j < DATA_SIZE; ++j)
    //     {
    //         output[i * DATA_SIZE + j] = result[j][0];
    //         complexOutput[i * DATA_SIZE + j] = result[j][1];
    //     }
    // }



    // // for (int i = 0; i < scales; ++i)
    // // {
    // //     for (int j = 0; j < DATA_SIZE; ++j)
    // //     {
    // //         value = Magnitude(output[i*DATA_SIZE + j], complexOutput[i*DATA_SIZE + j]);

    // //         fprintf(out_file, "%f\t", value);
    // //     }

    // //     fprintf(out_file, "\n");
    // // }


    // for (int i = 0; i < DATA_SIZE; ++i)
    // {
    //     fprintf(out_file, "%f\n", data_in[i][0]);
    // }

    // // for (int i = 0; i < DATA_SIZE; ++i)
    // // {
    // //     fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
    // //         output[0*DATA_SIZE + i] + 0., output[1*DATA_SIZE + i] + 5., output[2*DATA_SIZE + i] + 10,
    // //         output[3*DATA_SIZE + i] + 15, output[4*DATA_SIZE + i] + 20, output[5*DATA_SIZE + i] + 25,
    // //         output[6*DATA_SIZE + i] + 30, output[7*DATA_SIZE + i] + 35, output[8*DATA_SIZE + i] + 40,
    // //         output[9*DATA_SIZE + i] + 45);
    // // }

    // // for (int i = 0; i < DATA_SIZE; ++i)
    // // {
    // //     fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
    // //         conWindow[0*DATA_SIZE + i] + 0., conWindow[1*DATA_SIZE + i] + 5., conWindow[2*DATA_SIZE + i] + 10,
    // //         conWindow[3*DATA_SIZE + i] + 15, conWindow[4*DATA_SIZE + i] + 20, conWindow[5*DATA_SIZE + i] + 25,
    // //         conWindow[6*DATA_SIZE + i] + 30, conWindow[7*DATA_SIZE + i] + 35, conWindow[8*DATA_SIZE + i] + 40,
    // //         conWindow[9*DATA_SIZE + i] + 45);
    // // }
    

    // fclose(out_file);

    // //Sanitation Engineering
    // free(conWindow);
    // free(output);
    // free(complexOutput);

    // //FFTW sanitation. 
    // fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backwards);
    // fftw_free(data_in); fftw_free(fft_data); fftw_free(result);

    return 0;
}
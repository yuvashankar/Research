#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary arrays.
    double *data, *data_2D, *result, *scales, *frequency;

    double *con_result, *wavelet_result; 

    int n = DATA_SIZE;
    const int J = (int) MAX_I - MIN_I;

    // const int trials = 77;

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    data_2D =  (double*) malloc(n * J * sizeof(double));

    con_result =     (double*) malloc(n * J * sizeof(double));
    wavelet_result = (double*) malloc(n * J * sizeof(double));

    result  =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);

    TestCases(data, 8);

    // //Populate the data array
    // for (int i = 0; i < trials; ++i)
    // {
    //     TestCases(data, 2);
    //     for (int j = 0; j < n; ++j)
    //     {
    //         data_2D[i * n + j] = data[j];
    //     }
    // }

    CWT_Convolution(data, scales, n, J, 
                    con_result);
    
    Wavelet(data, scales, 
            FS, n, J,
            wavelet_result);

    // printf("Wavelet_result = %f, con_result = %f\n", wavelet_result[466920], con_result[466920]);

    for (int i = 0; i < n * J; ++i)
    {
        // wavelet_result[i] = log(wavelet_result[i]);

        result[i] = abs(wavelet_result[i] - con_result[i]);
        // result[i] = abs(wavelet_result[i] - con_result[i])/con_result[i];
    }

    double max = result[0];
    int array_index = 0;
    int freq_index = 0;

    for (int i = 0; i < J; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (result[i* n + j] > max)
            {
                max = result[i* n + j];
                
                array_index = j;
                freq_index = i;
            }
        }
        
    }

    double error_time = (double) array_index/FS;
    printf("Worst error at %f Hz at %f s\n", frequency[freq_index], error_time);



    // Compute the ERSP
    // ERSP (data_2D, scales, FS, n, J, trials, PAD_FLAG, 
    // result);

    // Write to file
    // char filename[] = "DATA.log";



    WriteFile(result, frequency, J, n, "DATA.log");
    // WriteGnuplotScript("dee Herro" , "DATA.log");
    // Plot(result, frequency,  n, J);

    //Free up Memory
    free(data_2D);
    free(data);  free(result);
    free(scales); free(frequency);

    free(con_result); free(wavelet_result);
    
    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}
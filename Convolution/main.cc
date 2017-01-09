#include "wavelet.h"
#include <math.h>
#include <omp.h>
#include <stdlib.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    // Initialize the necessary arrays.
    double *data, *result;
    double *scales, *frequency;

    int n = DATA_SIZE;
    const int J = (int) MAX_I - MIN_I;

    //Memory Allocations
    data          =  (double*) malloc(n *     sizeof(double));
    result        =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);

    // int data_size = 8;
    // int conSize = 4;

    // int data[data_size] = {0, 1, 2, 3, 4, 5, 6, 7};
    // int realResult[data_size] = {0};
    // int complexResult[data_size] = {0};

    // int conWindow[conSize] = {0, 1, 2, 3};
    // int complexWindow[conSize] = {0, 1, 2, 3};

    // for (int i = 0; i < data_size; ++i) //For every element in the data file
    // {
    //     realResult[i]    = 0.0;
    //     complexResult[i] = 0.0;
    //     for (int j = -conSize + 1; j < conSize; ++j) 
    //     {
    //         if ( (i - j) >= 0 && (i - j) < data_size)
    //         {
    //             // printf("data[%d - %d] = %d, conWindow[%d] = %d\n", i,j,  data[i - j], j, conWindow[-j]);   
    //             if (j >= 0)
    //             {
    //                 realResult[i]    += data[i - j] * conWindow[j];
    //                 complexResult[i] += data[i - j] * complexWindow[j];
    //             }
                    
    //             if (j < 0)
    //             {
                    
    //                 realResult[i]     += data[i - j] * conWindow[-j];
    //                 complexResult[i] -= data[i - j] * complexWindow[-j];
    //             }
    //         }
    //         else
    //         {
    //             // int count = (i - j);
    //             int count = ( data_size + (i - j) )%data_size;
    //             if (j >= 0)
    //             {
    //                 realResult[i]    += data[count] * conWindow[j];
    //                 complexResult[i] += data[count] * complexWindow[j];
    //             }

    //             if (j < 0)
    //             {
    //                 realResult[i]  += data[count] * conWindow[-j];
    //                 complexResult[i] -= data[count] * conWindow[-j];
    //             }
    //         }
    //     }
    // }

    // for (int i = 0; i < data_size; ++i)
    // {
    //     printf("realResult = %d, complexResult = %d\n", realResult[i], complexResult[i]);
    // }


    TestCases(data, 8);
    CWT_Convolution(data, scales, n, J, result);

    // Write to file
    WriteFile(result, frequency, J, n, "DATA.log");

    // Free up Memory
    free(data);  free(result);

    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}
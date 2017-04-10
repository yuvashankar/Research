#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(int argc, char *argv[])
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary arrays.
    double *data, *result;
    double *scales, *frequency;

    //Get File Size
    int sampling_frequency = atoi( argv[1] );
    // char file_name = argv[1];
    int n = GetFileSize( argv[2] );

    printf("sampling_frequency = %d, argv = %d, file_name = %s\n", sampling_frequency, atoi( argv[1] ), argv[2]);

    const int J = (int) MAX_I - MIN_I;

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    result  =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales    = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);
    assert(scales != NULL); assert(frequency!= NULL);

    int readNumber = ReadFile(data, argv[2] );
    assert (readNumber == n);
    
    printf("Computing Wavelet\n");
    Wavelet(data, scales, 
            FS, n, J,
            result);

    // //Populate the data array
    // for (int i = 0; i < trials; ++i)
    // {
    //     TestCases(data, 5);
    //     for (int j = 0; j < n; ++j)
    //     {
    //         data_2D[i * n + j] = data[j];
    //     }
    // }

    // // Compute the ERSP
    // ERSP (data_2D, scales, FS, n, J, trials, PAD_FLAG, 
    // result);

    printf("Plotting Result\n");
    WriteFile(result, frequency, J, n, "DATA.log");

    // Plot_PNG(result, frequency, J, n, "Continuous Wavelet Transform of a song", 
    // "making_water.png");

    printf("Done Wavelet Analysis\n");
    //Free up Memory
    free(data);  free(result);
    free(scales); free(frequency);
    
    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}
#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(int argc, char *argv[])
{
    //Start the timer!
    double t = omp_get_wtime();

    //Test Input Arguments
    if (argc < 2)
    {
        printf("Not enough input arguments\n");
        return -1;
    }

    //Initialize the necessary arrays.
    double *data, *result;
    double *scales, *frequency;
    double* stft_result;



    //Get the file name
    char file_name[255];
    char * dot_pointer = strchr(argv[2], '.');
    int dot_location = dot_pointer - argv[2] + 1;
    memcpy(file_name, argv[2], dot_location - 1);

    //Get File Size
    int sampling_frequency = atoi( argv[1] );
    
    int n = GetFileSize( argv[2] );
    // int n = 3 * sampling_frequency;
    

    const int J = (int) MAX_I - MIN_I;

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    result  =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales    = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);
    assert(scales != NULL); assert(frequency!= NULL);

    int readNumber = ReadFile( data, argv[2] );
    assert (readNumber == n);
    // TestCases( data, 8, 128.0 , sampling_frequency, n);

    stft_result = ShortTimeFourierTransform(data, sampling_frequency, n, WINDOW_SIZE);
    WriteSTFTFile(stft_result, WINDOW_SIZE/2, ceil( (double) n / WINDOW_SIZE), sampling_frequency, "STFT_Result.log");

    

    printf("Computing Wavelet\n");
    Wavelet(data, scales, 
            sampling_frequency, n, J,
            result);

    Find_Peaks(result, frequency, n, J);

    ARRAY_DATA global_max = Max(result, n * J);
    printf("n = %d, global_max = %.16f\n", n, global_max.value);

    printf("Plotting Result\n");

    Plot(result, frequency, J, n, 0, sampling_frequency,
        file_name,
        file_name);

    printf("Done Wavelet Analysis\n");
    //Free up Memory
    free(stft_result);

    free(data);  free(result);
    free(scales); free(frequency);
    
    
    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}
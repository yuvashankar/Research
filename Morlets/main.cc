 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

#include "wavelet.h"
#include <math.h>
#include <omp.h>
#include <float.h>
#include <gsl/gsl_statistics.h>

int main(int argc, char *argv[])
{
 //Start the timer!
    double t = omp_get_wtime();

    // Test Input Arguments
    if (argc < 2)
    {
        printf("Not enough input arguments\n");
        return -1;
    }

    // Initialize the necessary arrays.
    double *data, *result;
    double *scales, *frequency;

    //Get the file name
    char file_name[255];
    char * dot_pointer = strchr(argv[2], '.');
    int dot_location = dot_pointer - argv[2] + 1;
    memcpy(file_name, argv[2], dot_location - 1);

    // //Get File Size
    // int sampling_frequency = atoi( argv[1] );
    int sampling_frequency = 2048;
    double d_t = 1.0/ (double) sampling_frequency;
    
    // int n = GetFileSize( argv[2] );
    int n = 3 * sampling_frequency;
    
    const int J = (int) freq_to_scale(MIN_FREQUENCY, d_t) - freq_to_scale(MAX_FREQUENCY, d_t);
    assert (J > 0);

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    result  =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales    = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, d_t);
    frequency = IdentifyFrequencies(scales, J);

    assert(scales != NULL); assert(frequency!= NULL);

    // int readNumber = ReadFile( data, argv[2] );
    // assert (readNumber == n);
    TestCases( data, 6, 16.0 , sampling_frequency, n);

    printf("Computing Wavelet\n");
    
    Wavelet(data, scales, 
        sampling_frequency, n, J, 2,
        result);


    // Find_Peaks(result, frequency, sampling_frequency, n, J);

    // ARRAY_DATA global_max = Max(result, n * J);
    // printf("n = %d, global_max = %.16e\n", n, global_max.value);

    // printf("Plotting Result\n");
    Plot(result, frequency, J, n, 1, sampling_frequency,
        "CWT of x(t)",
        file_name);

    printf("Done Wavelet Analysis\n");

    free(data);  free(result);
    free(scales); free(frequency);
    
    // //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("Wavelet Execution Time: %e\n", t/1000);
    return 0;
}

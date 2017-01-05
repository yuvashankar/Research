#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary arrays.
    double *data, *result;
    double *scales, *frequency;
    double *conWindow, *complexWindow;
    double *realResult, *complexResult;
    
    int n = DATA_SIZE;
    const int J = (int) MAX_I - MIN_I;

    //Memory Allocations
    data          =  (double*) malloc(n *     sizeof(double));
    result        =  (double*) malloc(n * J * sizeof(double));
    conWindow     =  (double*) malloc(n *     sizeof(double));
    complexWindow =  (double*) malloc(n *     sizeof(double));
    realResult    =  (double*) malloc(n *     sizeof(double));
    complexResult =  (double*) malloc(n *     sizeof(double));
    
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);

    TestCases(data, 2);

    // double temp = 0.0;
    // for (int i = 0; i < n; ++i)
    // {
    //     complexWindow[i] = - CompleteComplexMorlet(temp, 1.0);
    //     temp+= DT;
    // }
    // WriteDebug(complexWindow, n, FS, "debug.log");

    int conSize = 0;
    for (int i = 0; i < J; ++i)
    {
        memset(conWindow, 0.0, n *     sizeof(double));
        memset(complexWindow, 0.0, n *     sizeof(double));
        memset(realResult, 0.0, n *     sizeof(double));
        memset(complexResult, 0.0, n *     sizeof(double));
        
        conSize = (int) W_0/(2 * M_PI * scales[i]);
        conSize *= 4;
        printf("Scale[%d] = %f, conSize = %d\n", i, scales[i], conSize);
        
        double temp = 0.0;
        //Populate Convolution windows.
        for (int j = 0; j < conSize; ++j)
        {
            conWindow[j] = CompleteRealMorlet(temp, scales[i]);
            complexWindow[j] = - CompleteComplexMorlet(temp, scales[i]);
            temp += DT;
        }

        for (int j = conSize; j < n; ++j)
        {
            conWindow[j] = 0.0;
            complexWindow[j] = 0.0;
        }

        Convolute(data, conWindow, complexWindow, n, conSize,
            realResult, complexResult);

        for (int j = 0; j < n; ++j)
        {
            result[i * n + j] = MAGNITUDE(realResult[j], complexResult[j]);
        }
    }

    // Write to file
    WriteFile(result, frequency, J, n, "DATA.log");

    //Free up Memory
    free(data);  free(result);
    free(scales); free(frequency);
    free(conWindow); free(complexWindow);
    free(realResult); free(complexResult);
    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}
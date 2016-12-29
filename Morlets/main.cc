#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary constants.
    double *data, *result, *scales, *frequency;

    int n = DATA_SIZE;
    const int J = MAX_I - MIN_I;

    //Memory Allocations
    data =           (double*) malloc(n *     sizeof(double));
    result =         (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL);

    //Get Scales and Frequencies
    scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY);
    frequency = IdentifyFrequencies(scales, J);

    //Populate the data array
    // TestCases(data, 5);
    for (int i = 0; i < n/2; ++i)
    {
        data[i] = 0.5;
        data[i + n/2] = 1.0;
    }
    

    int pad = CalculatePaddingSize(n, 2);
    // printf("n = %d, Padding Size = %d, 2 * n = %d\n", n, pad, 2 * n);

    fftw_complex * data_in  =           (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * pad );
    PopulateDataArray(data, n, pad, 2, data_in);
    // WriteDebug(data, pad, "debug.log");
        FILE* debug_log = fopen("debug.log", "w");
        for (int i = 0; i < pad; ++i)
        {
            double value = (double) i/pad;
            fprintf(debug_log, "%f\t%f\n", value, data_in[i][0]);
        }
        fclose(debug_log);



    

    // //Compute the ERSP
    // ERSP (data, scales, FS, n, J, 77, 
    // result);
    
    // Write to file
    // WriteFile(result, frequency, J, n, "DATA.log");
    // Plot(result, period,  n, J);

    //Free up Memory
    free(data);  free(result);
    free(scales); free(frequency);
    t = omp_get_wtime() - t;
    printf("Execution Time: %f\n", t);
    return 0;
}

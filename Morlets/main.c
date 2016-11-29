#include "wavelet.h"
#include "math.h"
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    double *data, *result, *wavelet_result, *baseline_out, *period;

    double t = omp_get_wtime();
    //Initialize the necessary constants.
    int n = DATA_SIZE;
    const int J = MAX_I - MIN_I;

    //Memory Allocations
    data =           malloc(n *     sizeof(double));
    result =         malloc(n * J * sizeof(double));
    wavelet_result = malloc(n * J * sizeof(double));
    baseline_out =   malloc(n * J * sizeof(double));
    period =         malloc(    J * sizeof(double));
    
    assert(data != NULL); assert(result != NULL); assert(period != NULL);
    assert(wavelet_result != NULL); assert(baseline_out != NULL);

    double * scales = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY);

    //populate the data array
    TestCases(data, 4);
    // n = ReadFile(data, "sst_nino3.dat");
    FILE* debug_out = fopen("debug.log", "w");
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        fprintf(debug_out, "%d\t%f\n", i, data[i]);
    }

    CleanData(data, DATA_SIZE);

    
    Wavelet(data, period, scales, 
        FS, n, J,
        wavelet_result);
    
    

    // RemoveBaseline(wavelet_result, n, J,
    //     1, FS,
    //     baseline_out);



    //Write to file
    WriteFile(wavelet_result, period, J, n, "DATA.log");

    //sanitation engineering
    free(data); free(result); free(period); free(wavelet_result); free(baseline_out);
    free(scales);
    t = omp_get_wtime() - t;
    printf("Execution Time: %f\n", t/10);
    return 0;
}

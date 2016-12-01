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
    TestCases(data, 5);
    // n = ReadFile(data, "sst_nino3.dat");
    


    CleanData(data, DATA_SIZE);

    
    Wavelet(data, period, scales, 
        FS, n, J,
        wavelet_result);

    double signalFrequency = 6.0/FS;
    double dw = 2 * M_PI * signalFrequency;

    // int conSize = (int) 1.0/signalFrequency;
    double TEST_SCALE = 2.0;
    double continuousPriode =  ( 2 * M_PI * TEST_SCALE) / ( W_0 ) ;
    int discretePeriod = (int) continuousPriode*FS;

    printf("Period = %d\n", discretePeriod);
    
    double * conWindow = malloc(discretePeriod * sizeof(double));
    double * complexWindow = malloc(discretePeriod * sizeof(double));
    
    // printf("CompleteComplexMorlet = %f\n", CompleteComplexMorlet(0.0, 1.0));
    // FILE* debug_out = fopen("debug.log", "w");
    double mTime = 0.0;
    
    for (int i = 0; i < discretePeriod; ++i)
    {
        conWindow[i] = CompleteRealMorlet(mTime, TEST_SCALE);
        complexWindow[i] = CompleteComplexMorlet(mTime, TEST_SCALE);
        // fprintf(debug_out, "%f\t%f\t%f\n", mTime, conWindow[i], complexWindow[i]);
        mTime += DT;
    }

    // fclose(debug_out);

    // RemoveBaseline(wavelet_result, n, J,
    //     1, FS,
    //     baseline_out);

    //Write to file
    WriteFile(wavelet_result, period, J, n, "DATA.log");

    //sanitation engineering
    free(data); free(result); free(period); free(wavelet_result); free(baseline_out);
    free(scales);
    t = omp_get_wtime() - t;
    printf("Execution Time: %f\n", t);
    return 0;
}

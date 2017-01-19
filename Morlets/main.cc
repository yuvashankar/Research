#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>


double CosineRealMorlet(double time, double scale)
{
    double norm = sqrt(scale);
    double w_o = FREQ * 2 * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + K_SIGMA * exp( -0.5 * w_o * w_o );
    double cosine = cos(w_o * time);

    double out = 0.5 * C_SIGMA * QUAD_ROOT_PI * norm * mor * cosine;
    return(out);
}

double CosineComplexMorlet(double time, double scale)
{
    double norm = sqrt(scale);
    double w_o = FREQ * 2 * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + K_SIGMA * exp( -0.5 * w_o * w_o );
    double sine = sin(w_o * time);

    double out = 0.5 * C_SIGMA * QUAD_ROOT_PI * norm * mor * sine;
    return(out);
}

double DiracRealMorlet(double time, double scale)
{
    double impulse = 2.0;
    time = (impulse -time )/scale;
    double norm = 1.0/sqrt(scale);
    double out = exp( - 0.5 * time * time ) * ( cos( W_0 * time ) - K_SIGMA );
    out = norm * C_SIGMA * QUAD_ROOT_PI * out;
    return(out);
}

double DiracComplexMorlet(double time, double scale)
{
    double impulse = 2.0;
    time = (impulse -time )/scale;
    double norm = 1.0/sqrt(scale);
    double out = exp( - 0.5 * time * time ) * ( sin( W_0 * time ) - K_SIGMA );
    out = norm * C_SIGMA * QUAD_ROOT_PI * out;
    return(out);
}


int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary arrays.
    double *data, *data_2D, *result, *scales, *frequency;

    double *con_result, *wavelet_result; 

    int n = DATA_SIZE;
    const int J = (int) MAX_I - MIN_I;

    const int trials = 1;

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    data_2D =  (double*) malloc(n * trials * sizeof(double));

    con_result =     (double*) malloc(n * J * sizeof(double));
    wavelet_result = (double*) malloc(n * J * sizeof(double));

    result  =  (double*) malloc(n * J * sizeof(double));
    assert(data != NULL); assert(result != NULL); 

    //Get Scales and Frequencies
    scales    = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, S0);
    frequency = IdentifyFrequencies(scales, J);

    TestCases(data, 8);
    // data[0] = 1.0;

    for (int i = 0; i < J; ++i)
    {
        double tau = 0.0;
        for (int j = 0; j < n; ++j)
        {
            // double realVal = DiracRealMorlet(tau, scales[i]);
            // double compVal = DiracComplexMorlet(tau, scales[i]); //Complex Conjucate

            double realVal = CosineRealMorlet(tau, scales[i]);
            double compVal = CosineComplexMorlet(tau, scales[i]);

            con_result[i * n + j] = MAGNITUDE(realVal, compVal);
            
            tau += DT;
        }
    }

    // //Compute the Convolution
    // CWT_Convolution(data, scales, n, J, 
    //                 con_result);

    // Compute the Wavelet Transform
    Wavelet(data, scales, 
    FS, n, J,
    wavelet_result);

    // // Compute the ERSP
    // ERSP (data_2D, scales, FS, n, J, trials, PAD_FLAG, 
    // wavelet_result);
    

    for (int i = 0; i < n * J; ++i)
    {
        // wavelet_result[i] = log(wavelet_result[i]);
        // con_result[i] = log(con_result[i]);

        result[i] = abs(wavelet_result[i] - con_result[i]);
        // result[i] = abs(result[i] - con_result[i])/con_result[i];
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

    array_index = freq_index * n + array_index;
    printf("Wavelet_result = %f, con_result = %f\n", wavelet_result[array_index], con_result[array_index]);
    printf("Worst Error = %.16f\n", result[array_index]);

    WriteFile(result, frequency, J, n, "DATA.log");

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
/** 
    \file CleanData.cc
    \brief This file contains the code to remove any noise in the input data. 
*/


// #include "processEEG.h"
#include <gsl/gsl_statistics.h>
#include <wavelet.h>

void CleanData(double * data, double n)
{
	double mean = gsl_stats_mean(data, 1, n);
    double sDeviation = gsl_stats_sd_m(data, 1, n, mean);
    // printf("Mean: %f, SD: %f\n", mean, sDeviation);

    //Compute the Z-Score or Standard Score
    for (int i = 0; i < n; ++i)
    {
        data[i] = (data[i] - mean)/sDeviation;
    }
}

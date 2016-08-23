/*CleanData.c

Takes a 1xn array and preforms the Z-Score Calculation

This was made into a function to clean up main, too many lines of code in main.

Input:
	data: this is a pointer array, the buffer WILL be modified
	n: the size of data
Output:
	data will be rewritten
*/
#include "processEEG.h"
#include <gsl/gsl_statistics.h>

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
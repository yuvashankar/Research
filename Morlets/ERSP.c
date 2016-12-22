//Functions to calculate ERSP
//Using the results from Grandchamp et al. to calculate the ESRP
//Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials
//Romain Grandchamp1,2* and Arnaud Delorme
//Implemented by Vinay Yuvashanakar
//McMaster University

#include "processEEG.h"
#include <gsl/gsl_statistics.h>

#define ZERO_TEST 0.00001

int RemoveBaseline(double* data, int num_of_samples, int J,
	int trials, double sampling_frequency,
	double* output)
{
	int m;
	
	//Tells us if the z-score was going to be divided by zero and we had to fix that. 
	int div_by_zero = 0;
	double * pre_stimulus;

	//Basically the cardinal of t = 0 to the stimulus
	m = PRE_EVENT_TIME * sampling_frequency;
	

	pre_stimulus = (double*) malloc( m * sizeof(double) );

	// for (int i = 0; i < m; ++i)
	// {
	// 	pre_stimulus[i] = data[i];
	// }

	// double mean = gsl_stats_mean(pre_stimulus, 1, m);
	// double sDeviation = gsl_stats_sd_m(pre_stimulus, 1, m, mean);
	// // printf("Mean = %f, Standard Deviation = %f\n", mean, sDeviation);

	// for (int i = 0; i < num_of_samples; ++i)
	// {
	// 	output[i] = (abs((data[i] * data[i])) - mean) / sDeviation;
	// }
	double value;
	//For every frequency of the datablock. 
	for (int i = 0; i < J; ++i)
	{
		//Copy the pre trial results from each frequency block into pre_stimulus.
		// memcpy(pre_stimulus, &data[ i * num_of_samples ], m); 
		for (int j = 0; j < m; ++j)
		{
			pre_stimulus[j] = data[i * num_of_samples + j];
		}
		
		//Calculate the mean and the standard deviation
		double mean = gsl_stats_mean(pre_stimulus, 1, m);
    	double sDeviation = gsl_stats_sd_m(pre_stimulus, 1, m, mean);
    	
    	//This ensures that there is no division by zero. 
    	if (abs(sDeviation) > ZERO_TEST)
    	{
		    for (int j = 0; j < num_of_samples; ++j)
		    {
		    	value = data[i * num_of_samples + j] * data[i * num_of_samples + j];
		        output[ i * num_of_samples + j] = (abs(value) - mean) / sDeviation;
		    }
    	}
	}

	free(pre_stimulus);
	return(div_by_zero);
}

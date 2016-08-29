//Functions to calculate ERSP
//Using the results from Grandchamp et al. to calculate the ESRP
//Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials
//Romain Grandchamp1,2* and Arnaud Delorme
// Implemented by Vinay Yuvashanakar
// McMaster University

#include "processEEG.h"
#include <gsl/gsl_statistics.h>

#define ZERO_TEST 0.00001

int RemoveBaseline(double* data, int num_of_samples, int J, 
	int trials, double sampling_frequency)
{
	int m;
	
	//Tells us if the z-score was going to be divided by zero and we had to fix that. 
	int div_by_zero = 0;
	double * pre_stimulus;


	//Basically the cardinal of t = 0 to the stimulus.
	m = PRE_EVENT_TIME * sampling_frequency;
	

	pre_stimulus = malloc( m * sizeof(double) );


	//For every frequency of the datablock. 
	for (int i = 0; i < J; ++i)
	{
		//Copy the pre trial results from each frequency block into pre_stimulus.
		// memcpy(pre_stimulus, &data[ i * num_of_samples ], m); 
		//I don't trust memcpy... this is slower.. but atleast I understand it. 
		for (int j = 0; j < m; ++j)
		{
			pre_stimulus[j] = data[i * num_of_samples + j];
		}
		
		//Calculate the mean and the standard deviation
		double mean = gsl_stats_mean(pre_stimulus, 1, m);
    	double sDeviation = gsl_stats_sd_m(pre_stimulus, 1, m, mean);
    	
    	// //This must be a problem for the lower scales because we deal with higher frequencies than that.
    	if (abs(sDeviation) < ZERO_TEST)
    	{
    		//preventing a division by zero here. 
    		sDeviation = 1;
    		div_by_zero = 1;
    	}

    	// printf("mean: %f, SD = %f i = %d\n ", mean, sDeviation, i);

    	//Remove the baseline from the calculation.
	    for (int j = 0; j < num_of_samples; ++j)
	    {
	        data[i*num_of_samples + j] = (data[i*num_of_samples + j] - mean) / sDeviation;
	    }

	}

	free(pre_stimulus);

	return(div_by_zero);
}

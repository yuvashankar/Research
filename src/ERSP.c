//Functions to calculate ERSP
//Using the results from Grandchamp et al. to calculate the ESRP
//Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials
//Romain Grandchamp1,2* and Arnaud Delorme
// Implemented by Vinay Yuvashanakar
// McMaster University

#include "processEEG.h"
#include <gsl/gsl_statistics.h>

void RemoveBaseline(double* data, double num_of_samples, int J, 
	int trials, double sampling_frequency)
{
	int m, sampleNumber;

	double * pre_stimulus; 

	//Basically the cardinal of t = 0 to the stimulus.
	m = PRE_EVENT_TIME * sampling_frequency;
	sampleNumber = (int) num_of_samples;

	pre_stimulus = malloc( m * sizeof(double) );


	//For every frequency of the datablock. 
	for (int i = 0; i < J; ++i)
	{
		
		//Copy the pre trial results from each frequency block into pre_stimulus.
		memcpy(pre_stimulus, &data[ i* sampleNumber ], m);
		
		//Calculate the mean and the standard deviation
		double mean = gsl_stats_mean(pre_stimulus, 1, m);
    	double sDeviation = gsl_stats_sd_m(pre_stimulus, 1, m, mean);
    	// printf("mean: %f, SD = %f\n", mean, sDeviation);
    	//Remove the baseline from the calculation.
	    for (int j = 0; j < num_of_samples; ++j)
	    {
	        data[i*sampleNumber + j] = (data[i*sampleNumber + j] - mean)/sDeviation;
	    }

	}

	free(pre_stimulus);

}
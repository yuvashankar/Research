//Functions to calculate ERSP
//Using the results from Grandchamp et al. to calculate the ESRP
//Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials
//Romain Grandchamp1,2* and Arnaud Delorme
//Implemented by Vinay Yuvashanakar
//McMaster University

#include "processEEG.h"
// #include "wavelet.h"
#include <gsl/gsl_statistics.h>


// #define ZERO_TEST 0.00001

int RemoveBaseline(double* data, int num_of_samples, int J,
	int trials, double sampling_frequency,
	double* output)
{
	//Initializations
	int div_by_zero = 0;
	size_t stride = 1;
	double * pre_stimulus;
	int i, j, k, m;
	double value, mean, sDeviation = 0.0;

	//Cardinal from t = 0 to the stimulus.
	m = PRE_EVENT_TIME * sampling_frequency;
	
	//Allocate Memory
	pre_stimulus = (double*) malloc( m * sizeof(double) );
	assert(pre_stimulus != NULL);
	
	for ( i = 0; i < J; ++i)
	{
		//Copy the pre trial results from each frequency block into pre_stimulus.
		memcpy(pre_stimulus, data + i*num_of_samples, sizeof(double) * num_of_samples);
		
		//Calculate mean and SD
		mean = gsl_stats_mean(pre_stimulus, stride, m);
    	sDeviation = gsl_stats_sd_m(pre_stimulus, stride, m, mean);

    	//Remove the Baseline
	    for ( k = 0; k < num_of_samples; ++k)
	    {
	    	value = data[i * num_of_samples + k] * data[i * num_of_samples + k];
	        output[i * num_of_samples + k] = (fabs(value) - mean) / sDeviation;
	    }
	}

	free(pre_stimulus);
	return(div_by_zero);
}

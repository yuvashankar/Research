//Functions to calculate ERSP
//Using the results from Grandchamp et al. to calculate the ESRP
//Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials
//Romain Grandchamp1,2* and Arnaud Delorme
//Implemented by Vinay Yuvashanakar
//McMaster University

#include "processEEG.h"
#include <gsl/gsl_statistics.h>

int ERSP (double * raw_data, double* scales, int sampling_frequency, int n, int J, int trials, 
	double * output)
{
	//Array Inits
	double * pre_stimulus, *wavelet_out, *baseline_out;
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;
	fftw_plan plan_forward, plan_backward;

	int i, j, x;

	//Calculate the necessary constants for the Continuous Wavelet Transform.
    const int PADDED_SIZE = CalculatePaddingSize(n, 1);
    const int m = PRE_EVENT_TIME * sampling_frequency;
    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SEC in Hz

    //Necessary Arrays for the CWT
    wavelet_out  = (double*) malloc( n * J * sizeof(double) );
    baseline_out = (double*) malloc( n * J * sizeof(double) );
    pre_stimulus = (double*) malloc( m     * sizeof(double) );

    //FFTW Memory Allocations
    data_in  =           (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data =           (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fftw_result  = 		 (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	//Populate the Data array
	PopulateDataArray(raw_data, n, PADDED_SIZE, data_in);

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	//Preapre for the plan backwards
	plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	/*Begin ERSP*/
	for ( x = 0; x < trials; ++x)
	{
		/*Begin Wavelet Analysis*/
		for (i = 0; i < J; ++i)
		{
			FrequencyMultiply(fft_data, PADDED_SIZE, scales[i], dw, 
				filter_convolution);

			//Take the inverse FFT and store it in fftw_result
			fftw_execute(plan_backward);

			//Calculate the power and store it in result this may need to be changed to accomodate for phase
			for (j = 0; j < n; ++j)
			{
				wavelet_out[i * n + j] = MAGNITUDE(fftw_result[j][0], fftw_result[j][1]);
			}
		}
		/*End Wavelet Analysis*/

		//Remove the baseline
		RemoveBaseline(pre_stimulus, wavelet_out,
			n, J, m,
			baseline_out);

		for ( i = 0; i < n * J; ++i)
		{
			output[i] += fabs(baseline_out[i]);
		}
	}
	
	for (int i = 0; i < n * J; ++i)
	{
		output[i] = output[i] / trials;
	}
	/*End ERSP*/

	//Sanitation Engineering
	fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
	fftw_free(data_in); fftw_free(fft_data); fftw_free(filter_convolution); fftw_free(fftw_result);
	free(pre_stimulus); free(baseline_out); free(wavelet_out);
	return(0);
}

int RemoveBaseline(double* pre_stimulus, double* pre_baseline_array, 
	const int n, const int J, const int m,
	double* output)
{
	double value;
	double mean, sDeviation;
	const int stride = 1;
	int i, j;
	// const int m = PRE_EVENT_TIME * sampling_frequency;

	for ( i = 0; i < J; ++i)
	{
		//Copy the pre trial results from each frequency block into pre_stimulus.
		memcpy(pre_stimulus, pre_baseline_array + i*n, sizeof(double) * n);
		
		//Calculate mean and SD
		mean = gsl_stats_mean(pre_stimulus, stride, m);
    	sDeviation = gsl_stats_sd_m(pre_stimulus, stride, m, mean);

    	//Remove the Baseline
	    for ( j = 0; j < n; ++j)
	    {
	    	value = pre_baseline_array[i * n + j] * pre_baseline_array[i * n + j];
	        output[i * n + j] = (fabs(value) - mean) / sDeviation;
	    }
	}

	return(0);
}

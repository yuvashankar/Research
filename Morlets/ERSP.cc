//Functions to calculate ERSP
//Using the results from Grandchamp et al. to calculate the ESRP
//Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials
//Romain Grandchamp1,2* and Arnaud Delorme
//Implemented by Vinay Yuvashanakar
//McMaster University

#include "processEEG.h"
// #include "wavelet.h"
#include <gsl/gsl_statistics.h>
#include <float.h>


// #define ZERO_TEST 0.00001


int ERSP (double * raw_data, double* scales, int sampling_frequency, int n, int J, int trials, 
	double * output)
{
	//Array Inits
	double * pre_stimulus, *wavelet_out, *baseline_out;
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;
	fftw_plan plan_forward, plan_backward;

	double value, mean, sDeviation;
	int i, j, k;
	int stride = 1;
	/*Begin Wavelet Analysis*/

	//Calculate Padding Required
    const int PADDED_SIZE = CalculatePaddingSize(n, 1);
    const int m = PRE_EVENT_TIME * sampling_frequency;

    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SEC in Hz

    wavelet_out  = (double*) malloc( n * J * sizeof(double) );
    baseline_out = (double*) malloc( n * J * sizeof(double) );
    pre_stimulus = (double*) malloc( m     * sizeof(double) );

    data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	
	filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );
	fftw_result  = 		 (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );

	PopulateDataArray(raw_data, data_in, n, PADDED_SIZE);

	// //populate the FFTW data vector. 
	// for (i = 0; i < n; ++i)
 //    {
 //    	data_in[i][0] = raw_data[i];
 //    	data_in[i][1] = 0.0;
 //    }

 //    //Force the rest of the data vector to zero just in case
 //    for (i = n; i < PADDED_SIZE; ++i)
 //    {
 //    	data_in[i][0] = 0.0;
 //    	data_in[i][1] = 0.0;
 //    }

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	
	
	//Preapre for the plan backwards
	plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
									 FFTW_BACKWARD, FFTW_ESTIMATE);

	for (int x = 0; x < trials; ++x)
	{
		for (i = 0; i < J; ++i)
		{
			//Calculate the corrosponding frequency to the scale
			// period[i] = (W_0)/(scales[i] * 2 * M_PI);

			//Compute the Fourier Morlet at 0 and N/2
			value = CompleteFourierMorlet(0.0, scales[i]);

			filter_convolution[0][0] = fft_data[0][0] * value;
			filter_convolution[0][1] = fft_data[0][1] * value;
			
			filter_convolution[PADDED_SIZE/2][0] = 0.0;
			filter_convolution[PADDED_SIZE/2][1] = 0.0;

			//Compute the Fourier Morlet Convolution in between
			for (j = 1; j < PADDED_SIZE/2 - 1; ++j)
			{
				value = CompleteFourierMorlet( j * dw , scales[i]);
				filter_convolution[j][0] = fft_data[j][0] * value;
				filter_convolution[j][1] = fft_data[j][1] * value;

				filter_convolution[PADDED_SIZE- j][0] = 0.0;
				filter_convolution[PADDED_SIZE- j][1] = 0.0;
			}

			//Take the inverse FFT. 
			fftw_execute(plan_backward);
		    
			//Calculate the power and store it in result
			for (j = 0; j < n; ++j)
			{
				wavelet_out[i * n + j] = MAGNITUDE(fftw_result[j][0], fftw_result[j][1]);
			}
		}
		/*End Wavelet Analysis*/

		/*Begin Baseline Removal*/
		for ( i = 0; i < J; ++i)
		{
			//Copy the pre trial results from each frequency block into pre_stimulus.
			// memcpy(pre_stimulus, wavelet_out + i*num_of_samples, sizeof(double) * num_of_samples);
			for (int j = 0; j < m; ++j)
			{
				pre_stimulus[j] = wavelet_out[i * n + j]; 
			}
			
			//Calculate mean and SD
			mean = gsl_stats_mean(pre_stimulus, stride, m);
	    	sDeviation = gsl_stats_sd_m(pre_stimulus, stride, m, mean);

	    	//Remove the Baseline
		    for ( k = 0; k < n; ++k)
		    {
		    	value = wavelet_out[i * n + k] * wavelet_out[i * n + k];
		        baseline_out[i * n + k] = (fabs(value) - mean) / sDeviation;
		    }
		}
		for ( i = 0; i < n * J; ++i)
		{
			output[i] += fabs(baseline_out[i]);
				// printf("Naan Alert!\n");
				// output[i] = DBL_MAX;
		}
	/*End Baseline Removal*/
	}

	/*Begin ERSP*/
	for (int i = 0; i < n * J; ++i)
	{
		output[i] = output[i] / trials;
	}
	/*End ERSP*/


	fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
	fftw_free(data_in); fftw_free(fft_data); fftw_free(filter_convolution); fftw_free(fftw_result);
	free(pre_stimulus); free(wavelet_out); free(baseline_out);
	return(0);
}





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

//Functions to calculate ERSP
//Using the results from Grandchamp et al. to calculate the ESRP
//Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials
//Romain Grandchamp1,2* and Arnaud Delorme
//Implemented by Vinay Yuvashanakar
//McMaster University

#include "processEEG.h"
#include <gsl/gsl_statistics.h>

int RemoveBaseline(double* pre_stimulus, double* pre_baseline_array, 
	const int n, const int J, const int sampling_frequency,
	double* output);

int FrequencyMultiply(const fftw_complex* fft_data, 
	const int data_size, const double scale, const double dw,
	 fftw_complex* filter_convolution);

int PopulateDataArray(double* input_data, const int data_size, const int padded_size, const int PAD_FLAG,
	fftw_complex* output_data);

int ERSP (double * raw_data, double* scales, int sampling_frequency, int n, int J, int trials, 
	double * output)
{
	//Array Inits
	double * pre_stimulus, *wavelet_out, *baseline_out;
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;
	fftw_plan plan_forward, plan_backward;

	int i, j, x;

	//Calculate the necessary constants for the Continuous Wavelet Transform.
    const int    PADDED_SIZE = CalculatePaddingSize(n, 1);
    const int    m           = PRE_EVENT_TIME * sampling_frequency;
    const double dw          = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SEC in Hz

    //Memory Allocations
    wavelet_out  = (double*) malloc( n * J * sizeof(double) );
    baseline_out = (double*) malloc( n * J * sizeof(double) );
    pre_stimulus = (double*) malloc( m     * sizeof(double) );

    //FFTW Memory Allocations
    data_in  =           (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data =           (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fftw_result  = 		 (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	//Populate the Data array
	PopulateDataArray(raw_data, n, PADDED_SIZE, 1, data_in);

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

	for ( i = 0; i < J; ++i)
	{
		//Copy the pre trial results from each frequency block into pre_stimulus.
		for ( j = 0; j < m; ++j)		
		{		
			pre_stimulus[j] = pre_baseline_array[i * n + j]; 		
		}

		//Calculate mean and standard deviation
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

int FrequencyMultiply(const fftw_complex* fft_data, 
	const int data_size, const double scale, const double dw,
	 fftw_complex* filter_convolution)
{
	int j; 
	double value; 
	//Compute the Fourier Morlet at 0 and N/2
	value = CompleteFourierMorlet(0.0, scale);

	filter_convolution[0][0] = fft_data[0][0] * value;
	filter_convolution[0][1] = fft_data[0][1] * value;
	
	filter_convolution[data_size/2][0] = 0.0;
	filter_convolution[data_size/2][1] = 0.0;

	//Compute the Fourier Morlet Convolution in between
	for (j = 1; j < data_size/2 - 1; ++j)
	{
		value = CompleteFourierMorlet( j * dw , scale);
		filter_convolution[j][0] = fft_data[j][0] * value;
		filter_convolution[j][1] = fft_data[j][1] * value;

		filter_convolution[data_size- j][0] = 0.0;
		filter_convolution[data_size- j][1] = 0.0;
	}

	return(0);
}
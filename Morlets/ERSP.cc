/** \file ERSP.cc
	
	\brief This file contains all of the function required to generate the Event Related Spectral Pertubation of EEG signals.
*/
#include "processEEG.h"
#include <gsl/gsl_statistics.h>

int RemoveBaseline(double* pre_stimulus, double* pre_baseline_array, 
	const int n, const int J, const int sampling_frequency,
	double* output);

int FrequencyMultiply(const fftw_complex* fft_data, 
	const int data_size, const double scale, const double dw,
	 fftw_complex* filter_convolution);

int PopulateDataArray(double* input_data, const int data_size, const int padded_size, const int padding_type,
	fftw_complex* output_data);

/**
	\fn int ERSP (double * raw_data, double* scales, int sampling_frequency, int n, int J, int trials, 
			double * output)

	\param raw_data A trials * n array containing the data to be analyzed
	\param scales A 1 x J array of the scales that the wavelets will be analyzed in
	\param sampling_frequency The frequency that the data was sampled in
	\param n The numer of samples in each data set
	\param J The number of scales to be analyzed
	\param trials The number of trials conducted for the ERSP
	\param output A n x J array with the resultant ERSP from all of the trials.

	This function conducts the Event Related Spectral Pertubation 
*/

int ERSP (double * raw_data, double* scales, const int sampling_frequency, const int n, 
	const int J, int const trials, const int padding_type, 
	double * output)
{
	//Array Inits
	double * pre_stimulus, *wavelet_out, *baseline_out;
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;
	fftw_plan plan_forward, plan_backward;

	int i, j, x;

	//Calculate the necessary constants for the Continuous Wavelet Transform.
    const int    PADDED_SIZE = CalculatePaddingSize(n, padding_type);
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
	PopulateDataArray(raw_data, n, PADDED_SIZE, padding_type, data_in);

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

/**
	\fn int RemoveBaseline(double* pre_stimulus, double* pre_baseline_array, 
			const int n, const int J, const int m,
			double* output)
	\brief A function that removes the pre stimulus noise found in EEG signals.

	\param pre_stimulus A 1 x m array to store the pre stimulus data
	\param pre_baseline_array An n x J array of the data that must be modified
	\param n The number of samples in the entire data array
	\param J The number of scales that were used
	\param m The size of the array before the stimulus
	\param output An n x J array that the function stores the result in. 

	This function follows the method outlined in the paper "Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials".

	The function will remove the baseline observed in in the pre stimulus by computing the z score on only the information before the stimulus. 
	The variable \a m is the number of samples before the stimulus was introduced. 

	All arrays must be pre allocated.
*/

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

/**
	\fn int FrequencyMultiply(const fftw_complex* fft_data, 
			const int data_size, const double scale, const double dw,
			fftw_complex* filter_convolution)
	\brief Multiples the signal with the wavelet at a specific scale in the frequency domain. 
	\param fft_data A fftw_complex * data_size array with the signal data in the frequency domain. 
	\param data_size The size of the data array
	\param scale THe scale of the wavelet that will be multipled with the signal array
	\param dw THe discrete increment in the frequency domain for the wavelet
	\param filter_convolution A fftw_complex * data_size array with the resulted multiplication

	This function mutliples the contents of fft_data with with the wavelet specified by the variable \a scale. 
	It stores the result in filter_convolution.

	All arrays must be pre allocated.
*/
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

/**
	\fn int PopulateDataArray(double* input_data, const int data_size, const int padded_size, const int padding_type,
			fftw_complex* output_data)
	\brief A function that copies and stores the input data in an array that is padded and friendly for FFTW.

	\param input_data A 1 x n array of the signal data
	\param data_size The size of the data array
	\param padded_size The size that the padded array needs to be
	\param PAD_FLAG The type of padding specified
	\param output_data  An fftw_complex * padded_size array that the data is stored in for FFTW to compute.

	This function takes the signal data from input_data and stores the result in an fftw_complex data array output_data.

	It returns the padding type
*/
int PopulateDataArray(double* input_data, const int data_size, const int padded_size, const int padding_type,
	fftw_complex* output_data)
{
	const double ramp = 2.0/data_size; // = 1.0/ n / 2
	double gain; 
	int i;

	int output_counter = 0;
	int input_counter = 0;

	switch(padding_type)
	{
		case 0: //No Padding what so ever
			//populate the FFTW data vector. 
			for (i = 0; i < data_size; ++i)
		    {
		    	output_data[i][0] = input_data[i];
		    	output_data[i][1] = 0.0;
		    }
		    break;
		
		case 1: //Zero - Padding
			//populate the FFTW data vector. 
			for (i = 0; i < data_size; ++i)
		    {
		    	output_data[i][0] = input_data[i];
		    	output_data[i][1] = 0.0;
		    }

		    //Force the rest of the data vector to zero just in case
		    for (i = data_size; i < padded_size; ++i)
		    {
		    	output_data[i][0] = 0.0;
		    	output_data[i][1] = 0.0;
		    }
		    break;
		
		case 2: //Duplicate array and ramp up and ramp down output
			for (i = 0; i < data_size; ++i)
		    {
		    	output_data[i][0] = input_data[i];
		    	output_data[i][1] = 0.0;
		    }

		    for (i = 0; i < (int) (0.5* data_size); ++i)
		    {
		    	output_counter = data_size + i;
		    	input_counter = (int) ((data_size/2) + i);
		    	gain = i * ramp;

		    	output_data[output_counter][0] = (1.0 - gain) * input_data[input_counter];
		    	output_data[output_counter][1] = 0.0;

		    	output_data[output_counter + (int) (0.5 * data_size)][0] = gain * input_data[i];
		    	output_data[output_counter + (int) (0.5 * data_size)][1] = 0.0;
		    }
			break;
		
		default: //Return just the array no padding. 
			//populate the FFTW data vector. 
			for (i = 0; i < data_size; ++i)
		    {
		    	output_data[i][0] = input_data[i];
		    	output_data[i][1] = 0.0;
		    }
		    break;
	}
    return(padding_type);
}
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

int PopulateDataArray(double* input_data, const int data_size, const int trial_number,
	const int padded_size, const int padding_type,
	fftw_complex* output_data);

/**
	\fn int ERSP (double * raw_data, double* scales, const int sampling_frequency, const int n, 
	const int J, int const trials, const int padding_type, 
	double * output)

	\param raw_data A trials * n array containing the data to be analyzed
	\param scales A 1 x J array of the scales that the wavelets will be analyzed in
	\param sampling_frequency The frequency that the data was sampled in
	\param n The numer of samples in each data set
	\param J The number of scales to be analyzed
	\param trials The number of trials conducted for the ERSP
	\param output A n x J array with the resultant ERSP from all of the trials.

	\return 0

	This function conducts the Event Related Spectral Pertubation of the given data set \a raw_data.
	It follows the method outlined by the paper: "Single-trial normalization for event-related spectral decomposition reduces sensitivity to noisy trials".

	This function uses the Continuous Wavelet Transform to generate the multi-resolution analysis of the given data. 

	This function is multi-threaded. 

	This function deals a lot with Fast Fourier Transforms, and can be optimized by using the Generate_FFTW_Wisedom() function. If no wisdom is provided, an approximate FFT algorithm will be used. 

	The variable raw_data must contain all of the data for each trial.

	raw_data, scales, and output must be pre-allocated. 
*/

int ERSP (double * raw_data, double* scales, const int sampling_frequency, const int n, 
	const int J, int const trials, const int padding_type, 
	double * output)
{
	int i, j, x;
	int number_of_threads = 1;
	//Calculate the necessary constants for the Continuous Wavelet Transform.
    const int    PADDED_SIZE = CalculatePaddingSize(n, padding_type);
    const int    m           = PRE_EVENT_TIME * sampling_frequency;
    const double dw          = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SEC in Hz

    fftw_init_threads();
    #pragma omp parallel num_threads(number_of_threads) private(i, j, x) shared(raw_data, output, scales, number_of_threads) default(none)
    {
    	//Array Inits
    	double * pre_stimulus, *wavelet_out, *baseline_out;
    	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;
    	fftw_plan plan_forward, plan_backward;

	    //Memory Allocations
	    wavelet_out  = (double*) malloc( n * J * sizeof(double) );
	    baseline_out = (double*) malloc( n * J * sizeof(double) );
	    pre_stimulus = (double*) malloc( m     * sizeof(double) );

	    //FFTW Memory Allocations
	    data_in            = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
		fft_data           = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
		filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
		fftw_result        = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

		#pragma omp critical (make_plan)
		{
			fftw_plan_with_nthreads(1);
			if (fftw_import_wisdom_from_filename("FFTW_Plan.wise") == 0)
			{
				printf("No FFTW Plan, using an unoptimized method\n");
				plan_forward  = fftw_plan_dft_1d(PADDED_SIZE, data_in,            fft_data,    FFTW_FORWARD,  FFTW_ESTIMATE);
				plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, FFTW_BACKWARD, FFTW_ESTIMATE);
			}
			else
			{
				plan_forward  = fftw_plan_dft_1d(PADDED_SIZE, data_in,            fft_data,    FFTW_FORWARD,  FFTW_EXHAUSTIVE);
				plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, FFTW_BACKWARD, FFTW_EXHAUSTIVE);	
			}
		}

		/*Begin ERSP*/
		#pragma omp for
		for ( x = 0; x < trials; ++x)
		{
			// memset(wavelet_out,        0.0, sizeof( double ) * n * J);
			// memset(baseline_out,       0.0, sizeof( double ) * n * J);
			// memset(pre_stimulus,       0.0, sizeof( double ) * m);

			memset(data_in,            0.0, sizeof( fftw_complex ) * PADDED_SIZE);
			memset(fft_data,           0.0, sizeof( fftw_complex ) * PADDED_SIZE);
			memset(filter_convolution, 0.0, sizeof( fftw_complex ) * PADDED_SIZE);
			memset(fftw_result,        0.0, sizeof( fftw_complex ) * PADDED_SIZE);

			/*Begin Wavelet Analysis*/
			PopulateDataArray(raw_data, n, x, 
							  PADDED_SIZE, padding_type, data_in);
			fftw_execute(plan_forward);

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
					wavelet_out[i * n + j] = fabs( wavelet_out[i * n + j] * wavelet_out[i * n + j] );
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

		#pragma omp for simd
		for ( i = 0; i < n * J; ++i)
		{
			output[i] /= trials;
		}
		/*End ERSP*/

		//Sanitation Engineering
		fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
		fftw_free(data_in); fftw_free(fft_data); fftw_free(filter_convolution); fftw_free(fftw_result);
		free(pre_stimulus); free(baseline_out); free(wavelet_out);

    }/*End of OpenMP*/

	return(0);
}

/**
	\fn int Generate_FFTW_Wisdom(int padded_size)
	
	\brief Analyzes the size of the FFTW arrays and generates the optimal plan. 
	
	\param padded_size The size of the FFT arrays. 

	\return 0 If successful
	\return 1 if unsuccessful

	This function can be used to optimize FFTW. This function will try to find the fastest FFT method based on the size of the array, and will store this information as "FFTW_plan.wise". 

	This function does not need to be used, but it can significantly improve performance if it is.
*/

int Generate_FFTW_Wisdom(int padded_size)
{
	int success_flag = 1;
	
	//Array Inits
	fftw_complex *data_in, *fft_data, *filter_convolution, *fftw_result;
	fftw_plan plan_forward, plan_backward;
	
	fftw_init_threads();
	
	//FFTW Memory Allocations
    data_in            = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * padded_size );
	fft_data           = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * padded_size );
	filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * padded_size );
	fftw_result        = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * padded_size );

	//FFTW Planning
	fftw_plan_with_nthreads(1);
	printf("Generating an Exhaustive FFTW Plan\n");
	plan_forward  = fftw_plan_dft_1d(padded_size, data_in,            fft_data,    FFTW_FORWARD,  FFTW_EXHAUSTIVE);
	plan_backward = fftw_plan_dft_1d(padded_size, filter_convolution, fftw_result, FFTW_BACKWARD, FFTW_EXHAUSTIVE);
	
	printf("Writing FFTW plan to FFTW_Plan.wise\n");
	if (fftw_export_wisdom_to_filename("FFTW_Plan.wise") != 0)
	{
		success_flag = 0;
	}

	fftw_destroy_plan(plan_forward); fftw_destroy_plan(plan_backward);
	fftw_free(data_in); fftw_free(fft_data); fftw_free(filter_convolution); fftw_free(fftw_result);

	return(success_flag);
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

	\return 0

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
		memset(pre_stimulus, 0.0, sizeof(double) * m);
		//Copy the pre trial results from each frequency block into pre_stimulus.
		for ( j = 0; j < m; ++j)		
		{		
			pre_stimulus[j] = pre_baseline_array[i * n + j]; 		
		}

		//Calculate mean and standard deviation
		mean       = gsl_stats_mean(pre_stimulus, stride, m);
    	sDeviation = gsl_stats_sd_m(pre_stimulus, stride, m, mean);

    	//Remove the Baseline
	    for ( j = 0; j < n; ++j)
	    {
	    	// value = pre_baseline_array[i * n + j] * pre_baseline_array[i * n + j];
	        output[i * n + j] = (pre_baseline_array[i * n + j] - mean) / sDeviation;
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
	\param scale THe scale of the wavelet that will be multiplied with the signal array
	\param dw THe discrete increment in the frequency domain for the wavelet
	\param filter_convolution A fftw_complex * data_size array with the resulted multiplication

	\return 0

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

	filter_convolution[0][0] = (fft_data[0][0]/data_size) * value;
	filter_convolution[0][1] = (fft_data[0][1]/data_size) * value;

	filter_convolution[data_size/2][0] = 0.0;
	filter_convolution[data_size/2][1] = 0.0;

	//Compute the Fourier Morlet Convolution in between
	for (j = 1; j < data_size/2 - 1; ++j)
	{
		value = CompleteFourierMorlet( j * dw , scale);
		filter_convolution[j][0] = (fft_data[j][0]/data_size) * value;
		filter_convolution[j][1] = (fft_data[j][1]/data_size) * value;

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

	\return padding_type

	This function takes the signal data from input_data and stores the result in an fftw_complex data array output_data.

	It returns the padding type
*/
int PopulateDataArray(double* input_data, const int data_size, const int trial_number,
	const int padded_size, const int padding_type,
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
		    	output_data[i][0] = input_data[trial_number * data_size + i];
		    	output_data[i][1] = 0.0;
		    }
		    break;
		
		case 1: //Zero - Padding
			//populate the FFTW data vector. 
			for (i = 0; i < data_size; ++i)
		    {
		    	output_data[i][0] = input_data[trial_number * data_size + i];
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
		    	output_data[i][0] = input_data[trial_number * data_size + i];
		    	output_data[i][1] = 0.0;
		    }

		    for (i = 0; i < (int) (0.5* data_size); ++i)
		    {
		    	output_counter = data_size + i;
		    	input_counter = (int) (trial_number * data_size + (0.5*data_size) + i);
		    	gain = i * ramp;

		    	output_data[output_counter][0] = (1.0 - gain) * input_data[input_counter];
		    	output_data[output_counter][1] = 0.0;

		    	output_data[output_counter + (int) (0.5 * data_size)][0] = gain * input_data[trial_number * data_size + i];
		    	output_data[output_counter + (int) (0.5 * data_size)][1] = 0.0;
		    }
			break;

		case 3: //Duplicate the signal once. 
			for ( i = 0; i < data_size; ++i)
			{
				output_data[i            ][0] = input_data[trial_number * data_size + i];
				output_data[i + data_size][0] = input_data[trial_number * data_size + i];

				output_data[i            ][1] = 0.0;
				output_data[i + data_size][1] = 0.0;
			}
			break;
		default: //Return just the array no padding. 
			//populate the FFTW data vector. 
			for (i = 0; i < data_size; ++i)
		    {
		    	output_data[i][0] = input_data[trial_number * data_size + i];
		    	output_data[i][1] = 0.0;
		    }
		    break;
	}
    return(padding_type);
}
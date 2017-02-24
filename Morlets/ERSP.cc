/** \file ERSP.cc
	
	\brief This file contains all of the function required to generate the Event Related Spectral Pertubation of EEG signals.
*/
#include "processEEG.h"
#include <gsl/gsl_statistics.h>

int ERSP (double * raw_data, double* scales, const int sampling_frequency, const int n, 
	const int J, int const trials, const int padding_type, 
	double * output)
{
	int i, j, x;
	int number_of_threads = 2;
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
					// wavelet_out[i * n + j] = fabs(wavelet_out[i * n + j] * wavelet_out[i * n + j]);
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
		mean       = gsl_stats_mean(pre_stimulus, stride, m);
    	sDeviation = gsl_stats_sd_m(pre_stimulus, stride, m, mean);

    	//Remove the Baseline
	    for ( j = 0; j < n; ++j)
	    {
	    	value = pre_baseline_array[i * n + j] * pre_baseline_array[i * n + j];
	        output[i * n + j] = (fabs(value) - mean) / sDeviation;
	        // output[i * n + j] = ( pre_baseline_array[i * n + j] - mean ) / sDeviation;
	    }
	}

	return(0);
}

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

		    //Horse the rest of the data vector to zero just in case
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
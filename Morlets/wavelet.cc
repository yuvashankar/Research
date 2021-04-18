 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

/**
	\file "Wavelet.cc"
	\brief This file contains all of the functions that support the ERSP and CWT functions
*/
#include "wavelet.h"
#include <omp.h>
#include <assert.h>
#include <gsl/gsl_math.h>	
#include <float.h>	
#include <gsl/gsl_statistics.h>

#define NORMALIZATION_FACTOR 9.1861911573363042e-02


int Wavelet(double* raw_data, double* scales, 
	double sampling_frequency, int n, int J, int padding_type,
	double* result)
{
	//Variable Declarations
	int i, j;
	fftw_complex *data_in, *fft_data;
	fftw_plan plan_forward;

	//Calculate Padding Required
    const int PADDED_SIZE = CalculatePaddingSize(n, padding_type);
    // const int PADDED_SIZE = n;

    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE);

    data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	PopulateDataArray(raw_data, n, 0, 
							  PADDED_SIZE, padding_type, data_in);

	fftw_init_threads();
	
	#pragma omp parallel private(i, j) shared (result, sampling_frequency, J, n, scales,  fft_data, data_in, plan_forward) default(none)
	{
		// double value;

		fftw_plan plan_backward;
		fftw_complex *filter_convolution, *fftw_result;

		filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
		fftw_result        = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

		#pragma omp critical (make_plan)
		{
			// FILE* wisdom = fopen("FFTW_Plan.wise", "r");
			// fftw_import_wisdom_from_filename("FFTW_Plan.wise");
			// Generate_FFTW_Wisdom(PADDED_SIZE);
			fftw_plan_with_nthreads(1);
			plan_forward  = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, FFTW_FORWARD, FFTW_ESTIMATE);
			plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
				FFTW_BACKWARD, FFTW_ESTIMATE);
			// fclose(wisdom);

			fftw_execute(plan_forward);
		}
		
	    #pragma omp for
		for (i = 0; i < J; ++i)
		{
			//Force the arrays to zero
			memset(filter_convolution, 0.0, sizeof( fftw_complex ) * PADDED_SIZE);
			memset(fftw_result,        0.0, sizeof( fftw_complex ) * PADDED_SIZE);

			//Multiply in the Frequency domain
			FrequencyMultiply(fft_data, PADDED_SIZE, scales[i], dw, 
								  filter_convolution);
			
			//Take the inverse FFT. 
			fftw_execute(plan_backward);
		    
			//Calculate the power and store it in result
			for (j = 0; j < n; ++j)
			{
				result[i * n + j] = MAGNITUDE(fftw_result[j][0], fftw_result[j][1]) / NORMALIZATION_FACTOR;
			}

		}
		
		//FFTW sanitation engineering. 
		fftw_destroy_plan(plan_backward);
	    fftw_free(fftw_result);
	    fftw_free(filter_convolution);
	}
	
	//Sanitation Engineering
	fftw_destroy_plan(plan_forward); 
	fftw_free(fft_data); fftw_free(data_in);
    return(0);
} /*Wavelet */

double* ShortTimeFourierTransform(double * raw_data, int n, int window_size)
{
	//Variable Declarations
	int i, j;
	fftw_complex *data_in, *fft_data;
	fftw_plan plan_forward;
	double* result;

	// FILE * STFT_FILE = fopen("STFT_Result.log", "w");

	int num_windows = ceil((double) n / window_size); 
	assert(num_windows > 0);

	// const int PADDED_SIZE = CalculatePaddingSize(window_size, 1);
	const int PADDED_SIZE = window_size;


	result =   (double*)             malloc(num_windows * PADDED_SIZE * sizeof(double));
	data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);

	for (i = 0; i < num_windows; ++i)
	{
		memset(data_in,  0.0, sizeof( fftw_complex ) * PADDED_SIZE);
		memset(fft_data, 0.0, sizeof( fftw_complex ) * PADDED_SIZE);

		//Fill Data into FFT Array
		for (j = 0; j < window_size; ++j)
		{
			if (i * window_size + j < n)
			{
				data_in[j][0] = raw_data[i * window_size + j];
				data_in[j][1] = 0.0;
			}
			else
			{
				data_in[j][0] = 0.0;
				data_in[j][1] = 0.0;
			}
		}
		
		fftw_execute(plan_forward);

		for (j = 0; j < window_size / 2; ++j)
		{
			result[j * num_windows + i] = MAGNITUDE(fft_data[j][0], fft_data[j][1]);
			// fprintf(STFT_FILE, "%d\t%f\t%d\n", j, result[j * num_windows + i], j * num_windows + i);
		}
	}

	// fclose(STFT_FILE);
	fftw_destroy_plan(plan_forward); 
	fftw_free(data_in); fftw_free(fft_data);

	return(result);
}

double* FFT(double * raw_data, int n)
{
	fftw_complex *data_in, *fft_data;
	const int PADDED_SIZE = n;
	
	double * result =   (double*)    malloc( PADDED_SIZE * sizeof(double) );
	data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	//Fill Data Array
	for (int i = 0; i < n; ++i)
	{
		data_in[i][0] = raw_data[i];
		data_in[i][1] = 0.0;
	}

	//Calculate the FFT of the data and store it in fft_data
	fftw_plan plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_execute(plan_forward);

	for (int i = 0; i < n/2; ++i)
	{
		result[i] = MAGNITUDE(fft_data[i][0], fft_data[i][1]);
	}

	return(result);
}



int CalculatePaddingSize(const int array_size, const int pad_flag)
{
	const int pad = ceil(log2(array_size));
	int out = array_size;
	switch(pad_flag)
	{
		case 0: //No Padding what so ever. 
			out = array_size;
			break;
		case 1: //Zero - Padding and preforming a Radix-2 Operation
			out = (int) pow(2, pad + 2);
    		break;
    	case 2: //Duplicate array and ramp up and ramp down output
    		out = 2 * array_size;
    		break;

    	case 3: //Repeat the array once. 
    		out = 2 * array_size;
    		break;

    	default: //Else return the array size
    		out = array_size;
    		break;
	}
	return(out);
}

double* GenerateScales(const double minimum_frequency, const double maximum_frequency, const double s_0)
{
	int min_i = freq_to_scale(maximum_frequency, s_0) + 1;
	int max_i = freq_to_scale(minimum_frequency, s_0);
	assert(min_i > 0); assert(max_i > 0);
	// printf("max_i = %d, min_i = %d\n", max_i, min_i);

	double * scales = (double*) malloc ( (max_i - min_i) * sizeof(double) );
	int count = ( max_i - min_i ) + 1;

	//Populate the scales array
	for (int i = 0; i < count; ++i)
	{
		int counterVariable = min_i + i;
		scales[i] = s_0 * pow(2, counterVariable * D_J);
	}
	return(scales);
}

int freq_to_scale(double frequency, double s_0)
{
	int out = floor( ( log2( (W_0) / (s_0 * 2 * M_PI * frequency) ) )/D_J);
	return(out);
}


double* IdentifyFrequencies(double* scales, int count)
{
	double * frequency = (double*) malloc( count * sizeof(double) );

	for (int i = 0; i < count; ++i)
	{
		frequency[i] = (W_0)/(scales[i] * 2 * M_PI);
	}
	return(frequency);
}

double CompleteFourierMorlet(double w, const double scale)
{
	// double norm = 1.0/sqrt(scale);
	double norm = sqrt(scale);
	w = w * scale; 
	double out = exp( -0.5 * ( W_0 - w ) * (W_0 - w) )
                                       - K_SIGMA * (exp ( -0.5 * w * w));

	out = norm * C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}

void TestCases(double *data, const int flag, double freq, double sampling_frequency, int data_size)
{
	// Fit a freq signal at two points
	double dt = 1./sampling_frequency;
	double fsig = freq/sampling_frequency;
	double dw = 2 * M_PI * fsig;
	double dampning_ratio = 0.014;

	double dwn = (2 * M_PI * dw)/ (sqrt(1 - dampning_ratio * dampning_ratio));
	// double dwn = sqrt(1 - dampning_ratio * dampning_ratio);
	// printf("dw = %f, dampning_ratio = %f, dwn = %f\n", dw, dampning_ratio, dwn);


	double w0 =  M_PI/4; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int) 1./fsig;
	double a = 10.990189448e-5;

	double frequency_increment = (MAX_FREQUENCY - MIN_FREQUENCY)/ 3.0; //3.0 seconds. 
	
	int t = (int) 2.0 * sampling_frequency; //At 2 seconds. 
	printf("t = %d\n", t);

	switch(flag)
	{
		//Impulse at T = 2 seconds
		case 1:		
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = (1.0/(a * sqrt(M_PI))) * exp( - pow( (i*dt - 2.0)/a, 2) );
			}

			// for (int i = 0; i < data_size; ++i)
			// {
			// 	data[i] = 0.0;
			// }

			// data[t] = 1.0;
			break;
		
		//Multiple Sines at t = 1500
		case 2:
			for (int i = data_size/2; i < data_size/2 + 2*one_peri; ++i)
			{
				data[i] = sin((i - data_size/2)* dw + w0) + sin((i - data_size/2)* 2* dw + w0);
			}
			break;

		//Multiple Sines at all times
		case 3:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = sin(i*dw + w0) + sin(i*2*dw + w0);
			}
			break;

		//Single sine at t = 1.0s;
		case 4: 
			for (int i = data_size/3; i < data_size/3 + 2 * one_peri; ++i)
			{
				data[i] = sin( (i - data_size/2) * dw);
			}
			break;

		
		case 5:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = cos( i * dw + w0 );
				
				if (i >= data_size/2 && i <= 2 * (data_size)/3)
				{
					data[i] = cos(i * dw + w0) + cos( 2.0*(i * dw+ w0) );

				}
			}
			break;

		case 6:
			for (int i = 0; i < data_size; ++i)
			{
				// data[i] = cos(i * dw + w0 );
				data[i] = 0.0;
				if (i >= data_size/3 && i <= data_size/2)
				{
					data[i] = cos(i * dw) + cos( 2.0 * (i * dw) );
				}
			}
			break;

		//Frequency Sweep
		case 7:
			for (int i = 0; i < data_size; ++i)
			{
				// fsig = frequency/sampling_frequency;
				// dw = 2*M_PI*fsig;

				data[i] = sin(w0 + 2 * M_PI * (MIN_FREQUENCY + (frequency_increment/2) * pow(i/sampling_frequency, 2)) );

				// frequency += frequency_increment; 
			}
			break;
			
		//Single sine all the way through. 
		case 8:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = cos(i*dw);
			}
			break;
		//Dampened Cosine
		case 9:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = exp( -dampning_ratio * i * dwn) * cos(i*dw + w0);
			}
			break;
	}
}

int WriteDebug(const double *data, const int length, const int sampling_frequency,
	const char* filename)
{
	FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    double t = 0.0;
    double dt = 1.0/sampling_frequency;

	for (int i = 0; i < length; ++i)
    {
    	// double value = (double) i/length;
    	fprintf(out_file, "%f\t%.16e\n", t, data[i]);
    	t += dt;
    }
    
    fclose(out_file);
    return 0;
}

int GetFileSize(char filename[])
{
	FILE* signalFile = fopen(filename, "r");
	assert(signalFile != NULL);
	// obtain file size:
	fseek (signalFile , 0 , SEEK_END);
	long lSize = ftell (signalFile);
	rewind (signalFile);

	char * buffer = (char*) malloc(sizeof(char)*lSize);
	assert(buffer != NULL);

	int result = fread (buffer, 1, lSize, signalFile);
	assert(result == lSize);

	char * token = strtok(buffer, "\n");
	
    //Get input from text.
	int counterVariable = 0;
	while (token !=NULL)
    {
    	// data[counterVariable] = atof(token);
    	counterVariable++;
        token = strtok (NULL, "\n");

    }
    fclose(signalFile);

    return (counterVariable);
}


int ReadFile(double data[], char filename[])
{
	FILE* signalFile = fopen(filename, "r");
	assert(signalFile != NULL);
	// obtain file size:
	fseek (signalFile , 0 , SEEK_END);
	long lSize = ftell (signalFile);
	rewind (signalFile);

	char * buffer = (char*) malloc(sizeof(char)*lSize);
	assert(buffer != NULL);

	int result = fread (buffer, 1, lSize, signalFile);
	assert(result == lSize);
	// puts(buffer);


	char * token = strtok(buffer, "\n");
	
    //Get input from text.
	int counterVariable = 0;
	while (token !=NULL)
    {
    	data[counterVariable] = atof(token);
    	counterVariable++;
        token = strtok (NULL, "\n");

    }
    fclose(signalFile);

    return (counterVariable);
}


double CWT_Cosine_Real(double time, double scale)
{
    double norm = 0.5 * sqrt(scale);
    
    double c_sigma = 1.0/ sqrt((1.0 + exp(- W_0 * W_0) - 2.0 * exp(-0.75 * W_0 * W_0)));
    double k_sigma = exp(-0.5 * W_0 * W_0);
    double quad_root_pi = pow(M_PI, -0.25);

    double w_o = FREQ * 2.0  * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + k_sigma * exp( -0.5 * w_o * w_o );
    double cosine = cos(w_o * time);

    double out = (c_sigma * quad_root_pi * norm * mor * cosine);
    return(out);
}

double CWT_Cosine_Complex(double time, double scale)
{
    double norm = 0.5 * sqrt(scale);
    double c_sigma = 1.0/ sqrt((1.0 + exp(- W_0 * W_0) - 2.0 * exp(-0.75 * W_0 * W_0)));
    double k_sigma = exp(-0.5 * W_0 * W_0);
    double quad_root_pi = pow(M_PI, -0.25);

    double w_o = FREQ * 2.0 * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + k_sigma * exp( -0.5 * w_o * w_o );
    double sine = sin(w_o * time);

    double out = (c_sigma * quad_root_pi * norm * mor * sine);
    return(out);
}

double CWT_Dirac_Real(double time, double scale)
{
    double norm = 1.0/sqrt(scale);
	// double norm = 1.0;

    double impulse = 2.0;
    time = (impulse - time )/scale;
    
    double out = exp( - 0.5 * time * time ) * ( cos( W_0 * time ) - K_SIGMA );

    out = (norm * C_SIGMA * QUAD_ROOT_PI * out);
    return(out);
}

double CWT_Dirac_Complex(double time, double scale)
{
    double norm = 1.0/sqrt(scale);
	// double norm = 1.0;

    double impulse = 2.0;
    time = (impulse - time )/scale;
    
    double out = exp( - 0.5 * time * time ) * ( sin( W_0 * time ) - K_SIGMA );

    out = (norm * C_SIGMA * QUAD_ROOT_PI * out);
    return(out);
}

void Numerical_Experiment(int n, int J, double dt, double fs)
{
	double* analytic_cwt_dirac, *analytic_cwt_cosine, *computed_cwt_dirac, *computed_cwt_cosine,
		*result_cwt_dirac, *result_cwt_cosine;

	double *dirac_signal, *cosine_signal;

	double* scales, *frequency;

	//Allocate All of the 2D Arrays
	analytic_cwt_dirac  =  (double*) malloc(n * J * sizeof(double));
	analytic_cwt_cosine =  (double*) malloc(n * J * sizeof(double));

	computed_cwt_dirac  =  (double*) malloc(n * J * sizeof(double));
	computed_cwt_cosine =  (double*) malloc(n * J * sizeof(double));

	result_cwt_dirac    =  (double*) malloc(n * J * sizeof(double));
	result_cwt_cosine   =  (double*) malloc(n * J * sizeof(double));


	//Allocate all of the 1D Arrays
	dirac_signal  = (double*) malloc(n * sizeof(double));
	cosine_signal = (double*) malloc(n * sizeof(double));

	//Things needed for Wavelet()
	scales    = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, dt);
  	frequency = IdentifyFrequencies(scales, J);

	//Populate the 1D Array
	TestCases( dirac_signal, 1, FREQ , fs, n);

	WriteDebug(dirac_signal, n, fs,
	"bruh.log");
	TestCases( cosine_signal, 8, FREQ, fs, n);

	//Compute the Wavelet Transform algorithmicallly
	Wavelet(dirac_signal, scales, 
		fs, n, J, 1,
		computed_cwt_dirac);

	Wavelet(cosine_signal, scales, 
		fs, n, J, 0,
		computed_cwt_cosine);

	//Compute the wavelet transform analytically
    for (int i = 0; i < J; ++i)
    {
        double tau = 0.0;

        for (int j = 0; j < n; ++j)
        {

            double dirac_realVal = CWT_Dirac_Real(tau, scales[i]);
            double dirac_compVal = - CWT_Dirac_Complex(tau, scales[i]);

            double cosine_realVal = CWT_Cosine_Real(tau, scales[i]);
            double cosine_compVal = CWT_Cosine_Complex(tau, scales[i]);

            analytic_cwt_dirac[i * n + j]  = MAGNITUDE(dirac_realVal, dirac_compVal);
            analytic_cwt_cosine[i * n + j] = MAGNITUDE(cosine_realVal, cosine_compVal);

            tau += dt;
        }
    }

    //Get the absolute error
    for (int i = 0; i < n * J; ++i)
    {
    	result_cwt_dirac[i] = fabs(computed_cwt_dirac[i] - analytic_cwt_dirac[i]);
    	result_cwt_cosine[i] =  fabs(computed_cwt_cosine[i] - analytic_cwt_cosine[i]);
    }

    ARRAY_DATA dirac_error, cosine_error, computed_dirac, analytic_dirac;

    dirac_error  = Max(result_cwt_dirac,  n * J);
    cosine_error = Max(result_cwt_cosine, n * J);
    computed_dirac = Max(computed_cwt_dirac, n * J);
    analytic_dirac = Max(analytic_cwt_dirac, n * J);


    printf("Computed Dirac Max = %.16g\n", computed_dirac.value);
    printf("Analytic Dirac Max = %.16g\n", analytic_dirac.value);


    printf("Maximum Absolute Error of Impulse is: %.16g\n", dirac_error.value);
    printf("Maximum Absolute Error of Cosine is: %.16g\n", cosine_error.value);

	// Plot(result_cwt_dirac, frequency, J, n, 0, fs,
	// 		"result_cwt_dirac",
	// 		"result_cwt_dirac");

	// Plot(result_cwt_cosine, frequency, J, n, 0, fs,
	// 		"result_cwt_cosine",
	// 		"result_cwt_cosine");

	// Plot(computed_cwt_dirac, frequency, J, n, 0, fs,
	// 		"computed_cwt_dirac",
	// 		"computed_cwt_dirac");

	// Plot(analytic_cwt_dirac, frequency, J, n, 0, fs,
	// 		"analytic_cwt_dirac",
	// 		"analytic_cwt_dirac");

	// Plot(computed_cwt_cosine, frequency, J, n, 0, fs,
	// 		"Computed CWT of a Cosine Function",
	// 		"computed_cwt_cosine");


	// Plot(analytic_cwt_cosine, frequency, J, n, 0, fs,
	// 		"Analytic CWT of a Cosine Function",
	// 		"analytic_cwt_cosine");

	free(dirac_signal); free(cosine_signal);
	free(result_cwt_dirac); free(result_cwt_cosine);
	free(computed_cwt_dirac); free(computed_cwt_cosine);
	free(analytic_cwt_dirac); free(analytic_cwt_cosine);

	free(scales); free(frequency);

}

void STFT_Timing_Test()
{
	//Start the timer!
    double t = omp_get_wtime();

    double * data;
    double *scales, *frequency;
    double* stft_result;

    const int sampling_frequency = 2048;

    const int n = 3 * sampling_frequency;
    const double d_t = 1.0/ (double) sampling_frequency;

    const int J = (int) freq_to_scale(MIN_FREQUENCY, d_t) - freq_to_scale(MAX_FREQUENCY, d_t);
    assert (J > 0);

    data    =  (double*) malloc(n *     sizeof(double));
    stft_result = (double*) malloc( WINDOW_SIZE * n * sizeof(double));

    scales    = GenerateScales(MIN_FREQUENCY, MAX_FREQUENCY, d_t);
    frequency = IdentifyFrequencies(scales, J);
    assert(scales != NULL); assert(frequency!= NULL);

    TestCases( data, 9, 16.0 , sampling_frequency, n);

    for (int i = 0; i < 10000; ++i)
    {
        stft_result = ShortTimeFourierTransform(data, n, WINDOW_SIZE);
    }

    WriteSTFTFile(stft_result, WINDOW_SIZE/2, ceil( (double) n / WINDOW_SIZE), sampling_frequency, "STFT_Result.log");

    free(scales); free(frequency);

    free(stft_result); free(data);

    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %e\n", t/10000.0);
}
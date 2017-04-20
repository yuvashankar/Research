/**
	\file "Wavelet.cc"
	\brief This file contains all of the functions that support the ERSP and CWT functions
*/
#include "wavelet.h"
#include <omp.h>
#include <assert.h>
#include <gsl/gsl_statistics.h>

#define TEST 0.00001

int Wavelet(double* raw_data, double* scales, 
	double sampling_frequency, int n, int J,
	double* result)
{
	//Variable Declarations
	int i, j;
	fftw_complex *data_in, *fft_data;
	fftw_plan plan_forward;


	//Calculate Padding Required
    const int PADDED_SIZE = CalculatePaddingSize(n, 1);
    // const int PADDED_SIZE = n;

    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SAMPLE in RAD/SEC

    data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );


	//populate the FFTW data vector
	for (i = 0; i < n; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	// data_in[i + n][0] = raw_data[i];

    	// data_in[i + n][1] = 0.0;
    	data_in[i    ][1] = 0.0;
    }

    //Force the rest of the data vector to zero just in case
    for (int i = n; i < PADDED_SIZE; ++i)
    {
    	data_in[i][0] = 0.0;
    	data_in[i][1] = 0.0;
    }

    double *temp = (double*) malloc(n * sizeof(double));
    FILE* debug_file = fopen("debug.log", "w");

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	// #pragma omp parallel num_threads(1) private(i, j) shared (result, sampling_frequency, J, n, scales,  fft_data) default(none)
	// {
		double value;

		fftw_plan plan_backward;
		fftw_complex *filter_convolution, *fftw_result;

		filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );
		fftw_result  = 		 (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );

		// #pragma omp critical (make_plan)
		// {
			//Preapre for the plan backwards
			plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
				FFTW_BACKWARD, FFTW_ESTIMATE);
		// }
		
	    // #pragma omp for
		for (i = 0; i < J; ++i)
		{
			//Force the arrays to zero
			memset(filter_convolution, 0.0, sizeof( fftw_complex ) * PADDED_SIZE);
			memset(fftw_result,        0.0, sizeof( fftw_complex ) * PADDED_SIZE);

			//Compute the Fourier Morlet at 0 and N/2
			value = CompleteFourierMorlet(0.0, scales[i]);

			filter_convolution[0][0] = ( fft_data[0][0] / PADDED_SIZE ) * value;
			filter_convolution[0][1] = ( fft_data[0][1] / PADDED_SIZE ) * value;
			
			filter_convolution[PADDED_SIZE/2][0] = 0.0;
			filter_convolution[PADDED_SIZE/2][1] = 0.0;

			//Compute the Fourier Morlet Convolution in between
			for (j = 1; j < PADDED_SIZE/2 - 1; ++j)
			{
				value = CompleteFourierMorlet( j * dw , scales[i]);
				filter_convolution[j][0] = ( fft_data[j][0] / PADDED_SIZE ) * value;
				filter_convolution[j][1] = ( fft_data[j][1] / PADDED_SIZE ) * value;

				// filter_convolution[j][0] = (fft_data[j][0]) * value;
				// filter_convolution[j][1] = (fft_data[j][1]) * value;
				filter_convolution[PADDED_SIZE- j][0] = 0.0;
				filter_convolution[PADDED_SIZE- j][1] = 0.0;
			}

			//Take the inverse FFT. 
			fftw_execute(plan_backward);
		    
			//Calculate the power and store it in result
			for (j = 0; j < n; ++j)
			{
				result[i * n + j] = MAGNITUDE(fftw_result[j][0], fftw_result[j][1]);
				temp[j] = result[i * n + j];
			}
			
			// double total = 0;
			// for (int j = 0; j < n; ++j)
			// {
			// 	total += temp[j];
			// }

			// double variance = gsl_stats_variance(temp, 1, n);
			// fprintf(debug_file, "%.16f\t%.16f\t%.16f\n", SCALE_TO_FREQ(scales[i]), variance, total);

		}

		free(temp);
		fclose(debug_file);

		//FFTW sanitation engineering. 
		fftw_destroy_plan(plan_backward);
	    fftw_free(fftw_result);
	    fftw_free(filter_convolution);
	// }
	//Sanitation Engineering
	fftw_destroy_plan(plan_forward); 
	fftw_free(fft_data); fftw_free(data_in);
    return(0);
} /*Wavelet */

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
			out = (int) pow(2, pad + 1);
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
	int min_i = FREQ_TO_SCALE(maximum_frequency) + 1;

	int max_i = FREQ_TO_SCALE(minimum_frequency);
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
	// double DT = 1./sampling_frequency;
	double fsig = freq/sampling_frequency;
	double dw = 2 * M_PI * fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;

	double frequency_increment = (MAX_FREQUENCY - MIN_FREQUENCY)/ 3.0; //3.0 seconds. 
	
	int t = 2 * sampling_frequency; //At 2 seconds. 

	switch(flag)
	{
		//Impulse at T = 2 seconds
		case 1:		
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = 0.0;
			}

			data[t] = 1.0;
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
				data[i] = sin( (i - data_size/2) * dw + w0 );
			}
			break;

		
		case 5:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = cos( i * dw + w0 );
				
				if (i >= data_size/2 && i <= 2 * (data_size)/3)
				{
					data[i] = 2 * cos(i * dw + w0);
				}
			}
			break;

		case 6:
			for (int i = 0; i < data_size; ++i)
			{
				data[i] = cos(i * dw + w0 );
				if (i >= data_size/3 && i <= data_size/2)
				{
					data[i] = cos(i * (dw - 0.005) + w0);
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
				data[i] = sin(i*dw + w0);
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

void FillData(double * data)
{
	// Fit a FREQ signal at two points
	// double DT = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2d   Pitch %.f   Discrete Period = %d \n",FS,FREQ,one_peri);

	for (int i = 0; i < DATA_SIZE; ++i)
	{
		data[i] = 0.0;
	}

	// //Impulse Sample
	// data[2000] = 1.0;
	int i;
	// double t=0;
	for(i=0;i<DATA_SIZE;i++){
		// data[i] = sin(i*dw) + sin(i*dw*4);
		data[i]=0.;
		// if((i>200)&(i<400))data[i]=sin( (i-200)*dw+w0);
		if((i>0.25*DATA_SIZE)&(i<0.25*DATA_SIZE+one_peri)) data[i]=sin( (i-200)*dw+w0);
		if((i>0.5*DATA_SIZE)&(i<0.5*DATA_SIZE+2*one_peri))data[i]=sin( (i-1000)*dw+w0);
		if((i>0.75*DATA_SIZE)&(i<0.75*DATA_SIZE+3*one_peri))data[i]=sin( (i-2000)*dw+w0);
	}
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
    double norm = sqrt(scale);
    double w_o = FREQ * 2 * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + K_SIGMA * exp( -0.5 * w_o * w_o );
    double cosine = cos(w_o * time);

    double out = 0.5 * C_SIGMA * QUAD_ROOT_PI * norm * mor * cosine;
    return(out);
}

double CWT_Cosine_Complex(double time, double scale)
{
    double norm = sqrt(scale);
    double w_o = FREQ * 2 * M_PI;
    w_o *= scale;
    double mor = exp( - 0.5 * (W_0 - w_o) * (W_0 - w_o) ) + K_SIGMA * exp( -0.5 * w_o * w_o );
    double sine = sin(w_o * time);

    double out = 0.5 * C_SIGMA * QUAD_ROOT_PI * norm * mor * sine;
    return(out);
}

double CWT_Dirac_Real(double time, double scale)
{
    double impulse = 2.0;
    time = (impulse - time )/scale;
    double norm = 1.0/sqrt(scale);
    double out = exp( - 0.5 * time * time ) * ( cos( W_0 * time ) - K_SIGMA );

    out = norm * C_SIGMA * QUAD_ROOT_PI * out;
    return(out);
}

double CWT_Dirac_Complex(double time, double scale)
{
    double impulse = 2.0;
    time = (impulse - time )/scale;
    double norm = 1.0/sqrt(scale);
    double out = exp( - 0.5 * time * time ) * ( sin( W_0 * time ) - K_SIGMA );
    out = norm * C_SIGMA * QUAD_ROOT_PI * out;
    return(out);
}
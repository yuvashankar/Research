/**
	\file "Wavelet.cc"
	\brief This file contains all of the functions that support the ERSP and CWT functions
*/
#include "wavelet.h"
#include <omp.h>

#define TEST 0.00001

/**
	\fn int Wavelet(double* raw_data, double* scales, 
			double sampling_frequency, int n, int J,
			double* result)
	\brief A function that computes the Continuous Wavelet Transform for the data given in \a raw_data
	\param raw_data A 1 x n array with the data required
	\param scales A 1 x J array with all of the scales for generating the wavelets
	\param sampling_frequency The sampling frequency of the given data
	\param n The size of the input data
	\param J The number of scales that is provided
	\param result An n x J array of contiguous memory that stores the result

	This function preforms the Continuous Wavelet Transform using Morlet Wavelets on the data given in raw_data. 
	It stores the result in the result array.

	This function only modifies the result array. The arrays must be pre allocated for this function to work. 

	You can provide the function with scales of your choosing, or one can generate dyadic scales with the GenerateScales() function.

	This function is optimized using openmp to allow for multi threading. 

*/
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

    const double dw = (2 * M_PI * sampling_frequency)/(PADDED_SIZE); //NOT IN RAD/SEC in Hz

    data_in  = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );
	fft_data = (fftw_complex *) fftw_malloc( sizeof( fftw_complex ) * PADDED_SIZE );

	//populate the FFTW data vector. 
	for (i = 0; i < n; ++i)
    {
    	data_in[i][0] = raw_data[i];
    	data_in[i][1] = 0.0;
    }

    //Force the rest of the data vector to zero just in case
    for (int i = n; i < PADDED_SIZE; ++i)
    {
    	data_in[i][0] = 0.0;
    	data_in[i][1] = 0.0;
    }

	//Calculate the FFT of the data and store it in fft_data
	plan_forward = fftw_plan_dft_1d(PADDED_SIZE, data_in, fft_data, 
									FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	#pragma omp parallel num_threads(2) private(i, j) shared (result, sampling_frequency, J, n, scales,  fft_data) default(none)
	{
		double value;

		fftw_plan plan_backward;
		fftw_complex *filter_convolution, *fftw_result;

		filter_convolution = (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );
		fftw_result  = 		 (fftw_complex *) fftw_malloc( sizeof( fftw_complex )* PADDED_SIZE );

		#pragma omp critical (make_plan)
		{
			//Preapre for the plan backwards
			plan_backward = fftw_plan_dft_1d(PADDED_SIZE, filter_convolution, fftw_result, 
				FFTW_BACKWARD, FFTW_ESTIMATE);
		}
		
	    #pragma omp for
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
				result[i * n + j] = MAGNITUDE(fftw_result[j][0], fftw_result[j][1]);
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

/**
	\fn int CalculatePaddingSize(int array_size, int pad_flag)
	
	\brief Calculates the size that the padded array should be. 

	\param array_size The cardinal or size of the signal sample
	\param pad_flag The type of padding required: see PAD_FLAG

	\return paadded_size The cardinal or size that the padded array should be

	This function computes the size of the padded array depending on the type of padding specified. 
	It takes the size of the data array, and type of pad, and returns how large the padded array should be.
*/
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
    	default: //Else return the array size
    		out = array_size;
    		break;
	}
	return(out);
}

/**
	\fn double* GenerateScales(const double minimum_frequency, const double maximum_frequency, const double s_0)
	
	\brief This function generates the scales that will be used in the Continuous Wavelet Transform.
	
	\param minimum_frequency The lowest frequency that must be observed
	\param maximum_frequency The higest frequency that must be observed
	\param s_0 The smallest scale usually it is \f$ 2 * \delta t \f$

	\return scales An 1 x n array with the dyadic scales.

	This function computes the dyadic scales to be generated to accurately compute the multi resolution analysis of a signal. 
	Given the minimum frequency and the maximum frequency, the function will generate a 1 x n array with the scales necessary scale factors for the Continuous Wavelet Transform

	The Scales array will be allocated in this function, so it is wise to deallocate this array after it is used. 
*/

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

/**
	\fn double* IdentifyFrequencies(double* scales, int count)
	
	\brief Compute the corrosponding frequency to scales used.

	\param scales A 1 x count array of the scales used
	\param count The cardinal of the scales array

	\return frequency The corrosponding frequency for each scale in the scale array

	This function computes the corrosponding frequency for each scale provided in the scales array. 

	It allocates memory and returns the allocated array

	It is wise to dealloate this array after use with the free() function. 

*/
double* IdentifyFrequencies(double* scales, int count)
{

	double * frequency = (double*) malloc( count * sizeof(double) );

	for (int i = 0; i < count; ++i)
	{
		frequency[i] = (W_0)/(scales[i] * 2 * M_PI);
	}
	return(frequency);

}

// double FourierMorlet(double w, double scale, double normal)
// {
// 	double exponent = -0.5 * (scale * w - W_0) * (scale * w - W_0);
// 	double out = QUAD_ROOT_PI* normal * exp(exponent);
// 	return(out);
// }

/**
	\fn double CompleteFourierMorlet(const double w, const double scale)
	
	\brief Computes the Morlet Wavelet in the frequency domain. 

	\param w
	\param scale

	\return morlet

	This function generates the Morlet Wavelet in the frequency domain normalized by the scale. 

	The formula computed is 
	\f[
		\hat{\Psi}_\sigma(\omega) = c_\sigma \pi^{-\frac{1}{4}}(e^{-\frac{1}{2}(\sigma - \omega)^2} - \kappa_\sigma e^{-\frac{1}{2} \omega^2})
	\f]
*/
double CompleteFourierMorlet(const double w, const double scale)
{
	// double W_0 = 6.0;
	// double C_SIGMA = pow( (1.0 + exp(-W_0_2) - 2.0 * exp(-0.75 * W_0_2)), -0.5 );
	// double K_SIGMA = exp( -0.5 * W_0_2 );

	// double norm = pow(scale, -0.5);
	double norm = 1.0/sqrt(scale);
	double out = exp( -0.5 * ( W_0 - scale * w ) * (W_0 - scale * w) ) 
					- K_SIGMA * (exp ( -0.5 * scale * w * w));
	out = norm * C_SIGMA * QUAD_ROOT_PI * out;
	return(out);
}

/**
	\fn void TestCases(double *data, const int flag)
	\brief Generates a suite of test case data for wavelet analysis

	\param data The 1 x n data array to be populated
	\param flag The type of test data to be generated

	This function populates the data array with 3 seconds of sample data. The \a flag parameter specifies the type of test data that will be generated

	<table>
	<caption id="multi_row">TestCases Flags</caption>
	<tr><th> Flag Type 			<th> Output
	<tr><td> 1         			<td> Impulse at T = 2 seconds
	<tr><td> 2         			<td> 2 Sine waves at t = 1.5 seconds at FREQ and 2 * FREQ
	<tr><td> 3 					<td> 2 sine waves at FREQ and 2 * FREQ from t = 0 to 3 s
	<tr><td> 4 					<td> Single sine wave at t = 1.0 s
	<tr><td> 5 					<td> sin(x) from t = 0.0 to 3.0 and 2*sin(x) from t = 1.5 - 2.0 s
	<tr><td> 6 					<td> sin(x) from t = 0.0 - 3.0 and sin(x - w0) where w0 = 0.005 from t = 1.0 s - 1.5 s
	<tr><td> 7 					<td> Frequency Sweet from MIN_FREQUENCY to MAX_FREQUENCY
	</table>
*/

void TestCases(double *data, const int flag)
{

	// Fit a FREQ signal at two points
	// double DT = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	// printf("FS  %.2d   Pitch %.f   Discrete Priode = %d \n",FS,FREQ,one_peri);

	// double frequency = MIN_FREQUENCY;
	double frequency_increment = (MAX_FREQUENCY - MIN_FREQUENCY)/ 3.0; //3.0 seconds. 
	
	int t = 2 * FS; //At 2 seconds. 

	switch(flag)
	{
		//Impulse at T = 2 seconds
		case 1:
			
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = 0.0;
			}

			data[t] = 1.0;
			break;
		
		//Multiple Sines at t = 1500
		case 2:
			for (int i = DATA_SIZE/2; i < DATA_SIZE/2 + 2*one_peri; ++i)
			{
				data[i] = sin((i - DATA_SIZE/2)* dw + w0) + sin((i - DATA_SIZE/2)* 2* dw + w0);
			}
			break;
		//Multiple Sines at all times
		case 3:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = sin(i*dw + w0) + sin(i*2*dw + w0);
			}
			break;
		//Single sine at t = 1.0s;
		case 4: 
			for (int i = DATA_SIZE/3; i < DATA_SIZE/3 + 2 * one_peri; ++i)
			{
				data[i] = sin( (i - DATA_SIZE/2) * dw + w0 );
			}
			break;
		//Single sine all the way through. 
		case 5:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = cos( i * dw + w0 );
				if (i >= DATA_SIZE/2 && i <= 2 * (DATA_SIZE)/3)
				{
					data[i] = 2 * cos(i * dw + w0);
				}
			}
			break;
		case 6:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				data[i] = cos(i * dw + w0 );
				if (i >= DATA_SIZE/3 && i <= DATA_SIZE/2)
				{
					data[i] = cos(i * (dw - 0.005) + w0);
				}
			}
			break;
		//Frequency Sweep
		case 7:
			for (int i = 0; i < DATA_SIZE; ++i)
			{
				// fsig = frequency/FS;
				// dw = 2*M_PI*fsig;

				data[i] = sin(w0 + 2 * M_PI * (MIN_FREQUENCY + (frequency_increment/2) * pow(i/FS, 2)) );

				// frequency += frequency_increment; 
			}
			break;
	}
}

/**
	\fn int WriteFile(const double *data, const double *period, const int x, const int y, const char* filename)

	\brief A function that writes the Wavelet Results to the disk. 

	\param data A x x y array with the  data that is going to be written 
	\param period A 1 x y array with the frequencies that were analyzed
	\param x The number of samples in the signal
	\param y The number of frequencies analyzed
	\param filename The name of the file that will be written

	\return 0 if successful
	\return -1 if unsuccessful

	This function will write the resultant data computed by Wavelet() and ERSP() into the disk so that it can be graphed by Gnuplot. 
	One can plot the output of this function using the matrix.gplot file. 
*/

int WriteFile(const double *data, const double *frequency, const int x, const int y, 
	const char* filename)
{

    FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    //Xticks
    fprintf(out_file, "%d\t", x);
    for (int i = 0; i < y; ++i)
    {
    	fprintf(out_file, "%f\t", (double) i/FS);
    }
    fprintf(out_file, "\n");

    double small_eps = 0.00001; //Add a small eps so that logs of zero don't happen. 
	for (int i = 0; i < x; ++i)
    {
    	//Feed Frequency
    	fprintf(out_file, "%f\t", frequency[i]);

    	//Feed Data
        for (int j = 0; j < y; ++j)
        {
            // value = Magnitude(result[i*n + j], result[i*n + j]);
            fprintf(out_file, "%.5f\t", data[i*y + j] + small_eps);
        }
        //Ready for the next line.
        fprintf(out_file, "\n");

    }

    fclose(out_file);
    return(0);
}

/**
	\fn int WriteDebug(const double *data, const int length, const char* filename)
	
	\brief A function that writes a 1 - d matrix into a log file

	\param data A 1 - dimentional data array containing the data to be written
	\param length The size of the data array
	\param filename THe name of the file to be written

	\return 0 if successful
	\return -1 if unsuccessful

	This function writes a 1 - dimentional array to the disk, it's useful when trying to quickly get the results from an array. 
*/


int WriteDebug(const double *data, const int length, const char* filename)
{
	FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

	for (int i = 0; i < length; ++i)
    {
    	double value = (double) i/length;
    	fprintf(out_file, "%f\t%f\n", value, data[i]);
    }
    
    fclose(out_file);
    return 0;
}

/**
    \fn void FillData(double * data)
	\brief Populates the input data array with a 3 sparse sine waves. 
	\param data A 1 - dimentional block of memory that will be overwritten. 

	Similar to TestCases()
*/
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

/**
	\fn int ReadFile(double data[], char filename[])

	
	\param A pre allocated 1 dimentional array
	\param filename The name and location of the file to be opened

	\return array_size The number of elements that was read

	This function opens a file, and reads the input assuming that the file is stored with one value at every line.

*/
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

/*! \file wavelet.h
    \brief The supporting header file for generating the Continuous Wavelet Transform
*/
//Includes
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#include <omp.h>

typedef struct 
{
	double value;
	int    index;
} ARRAY_DATA;

//Global Constants

/*! \def QUAD_ROOT_PI 
\brief the quad root of pi (\f$ \pi^{-0.25} \f$) computed to machine precision
*/
#define QUAD_ROOT_PI 0.75112554446494248

/** 	
	\var C_SIGMA 
 	\brief A constant needed to compute the Morlet Wavelet precalculated to machine precision

 	\f[ C_\sigma = (1 + e^{-\sigma^2} - 2e^{-\frac{3\sigma^2}{4}})^{-\frac{1}{2}} \f]
*/
#define C_SIGMA 1.0000000000018794


/**
	\var K_SIGMA
	\brief A constant needed to compute the Morlet Wavelet precalculated to machine eps

	\f[ \kappa_\sigma = e^{-\frac{\sigma^2}{2}} \f]
*/
#define K_SIGMA 1.522997974471262843e-08


/**
	\var W_0 
	\brief The fundamental frequency of the Morlet Wavelet
*/
#define W_0 6.0


/**
	\var W_0_2
	\brief The fundamental frequency of the Morlet Wavelet squared
*/
#define W_0_2 36.0


/**
	\var D_J
	The amount of "sub octaves" or sub scales inbetween the major scales that will be used. 
	The lower the number, the higher the resolution of the result.
*/
#define D_J 0.03125

/**
	\var PAD_FLAG
	\brief The type of padding specified for the Continuous Wavelet Transform. 

	The Continuous Wavelet Transform uses the Fast Fourier Transform. It is sometimes efficient to add additional values to the edge of the data array to improve the speed of the FFT. This method is commonly called padding the data array. 
	This variable determins the type of padding that will be used to assist the Fast Fourier Transform. THe padding options are:

	0 - No Padding
	- The array will be analyzed with no padding. 

	1 - Zero Padding
	- The size of the array will be enlarged to the closest power of two and zeros will be added to the end.

	2 - Ramp Padding
	- The array will be doubled in size, and the signal will be ramped up and ramped down to gradually. 

	3 - Repetion Padding
	- The signal will be repeated once, in effect the signal length will be doubled. 

	If none of these are specified, the array is not padded by default.
*/
#define PAD_FLAG 1

//Signal Constants

/*!
	\var FS
	\brief Used by TestCases() to generate sample data
*/
#define FS 500

/**
	\var DT 
	\brief \f$ \delta t = \frac{1}{f_s} \f$
*/
#define DT 1.0/FS

/**
	\var S0 
	\brief the lowest scale that can be used to compute the CWT \f$ s_0 = 2 \delta t \f$
*/
#define S0 DT

#define FREQ 128
#define DATA_SIZE 6144

//Plotting Constants
/**
	\var MAX_FREQUENCY
	\brief The maximum frequency that will be analyzed
*/
#define MAX_FREQUENCY FS/2

/**
	\var MIN_FREQUENCY
	\brief The minimum frequency that will be analyzed
*/
#define MIN_FREQUENCY 0.5
	
#define MAX_DATA_SIZE 10000000

#define MIN_I FREQ_TO_SCALE(MAX_FREQUENCY)
#define MAX_I FREQ_TO_SCALE(MIN_FREQUENCY)

//Macros
/*!
\def FREQ_TO_SCALE(x)
\brief Converts a given frequency \a x to a scale, handy for debugging. Note the scale is divided into sub octaves
*/
#define FREQ_TO_SCALE(x) floor( ( log2( (W_0) / (S0 * 2 * M_PI * x) ) )/D_J)

/*!
    \def SCALE_TO_FREQ(x)
	\brief Converts a given scale \a x to its corrosponding frequency. 
*/
#define SCALE_TO_FREQ(x) (W_0)/(x * 2 * M_PI)

/*!
	\def MAGNITUDE(x, y)
	\brief Computes the 2- norm or the x ^ 2 + y ^ 2, of \a x and \a y
*/
#define MAGNITUDE(x,y) sqrt( (x * x) + (y * y) )

/**
    \fn void FillData(double * data)
	\brief Populates the input data array with a 3 sparse sine waves. 
	\param data A 1 - dimentional block of memory that will be overwritten. 

	Similar to TestCases()
*/
// void FillData(double * data);

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
void TestCases(double *data, const int flag, double freq, double sampling_frequency, int data_size);

/**
	\fn int ReadFile(double data[], char filename[])

	\param A pre allocated 1 dimentional array
	\param filename The name and location of the file to be opened

	\return array_size The number of elements that was read

	This function opens a file, and reads the input assuming that the file is stored with one value at every line.

*/
int  ReadFile(double data[], char filename[]);

int GetFileSize(char filename[]);

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

int  WriteFile(const double *data, const double *frequency, const int x, const int y, int sampling_frequency,
	const char filename[]);

/**
	\fn int WriteGnuplotScript(const char graph_title[], const char filename[])

	\param graph_title The title that the graph should be in
	\param filename The name of the file that will be graphed

	This function will generate a file called "script.gplot" that will can be used to 
	plot data generated using WriteFile()
	\returns 0 on success
*/
int WriteGnuplotScript(const char graph_title[], const char filename[]);

/**
	\fn int Plot_PNG(double * data, double * periods, int num_x, int num_y, char graph_title[], 
					const char filename[])
	\brief Uses PNGWriter to plot the results of the Continuous Wavelet Transform.

	\param data An num_x * num_y size data array that will be plotted
	\param periods A num_y x 1 array of the frequency used to analyze the signal data
	\param num_x The size of the data
	\param num_y The number of frequencies analyzed
	\param graph_title The title of the graph
	\param filename The name of the file that the data will be written to. 

	This function plots a two dimentional array into a png using PNGWriter. It utilizes GetRange() and GetColour() to function.
	In order to reduce the horizontal length of the image, a stride variable was introduced that will skip a certain number of elements in the data array. THe higher this number the slimmer the resultant graph will be. 
	
	The lines_size plots the same pixel in a number of vertical columns. THis allows the vertical scaling to be increased and decreased. 
*/
int Plot_PNG(double * data, double * periods, int num_x, int num_y, char graph_title[], 
	const char filename[]);

int Plot(double * data, double * frequency, int num_x, int num_y, int plot_type, int sampling_frequency,
	char graph_title[],
	char filename[]);

/**
	\fn int WriteDebug(const double *data, const int length, const int sampling_frequency,
	const char* filename)
	
	\brief A function that writes a 1 - d matrix into a log file

	\param data A 1 - dimentional data array containing the data to be written
	\param length The size of the data array
	\param sampling_frequency How often the data array was sampled
	\param filename THe name of the file to be written

	\return 0 if successful
	\return -1 if unsuccessful

	This function writes a 1 - dimentional array to the disk, it's useful when trying to quickly get the results from an array. 
*/
int WriteDebug(const double *data, const int length, const int sampling_frequency,
	const char* filename);

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
	double * output);

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
double CompleteFourierMorlet(double w, const double scale);

double CompleteRealMorlet (double x, double scale);
double CompleteComplexMorlet(double x, double scale);


double CWT_Cosine_Real(double time, double scale);
double CWT_Cosine_Complex(double time, double scale);
double CWT_Dirac_Real(double time, double scale);
double CWT_Dirac_Complex(double time, double scale);

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
	double sampling_frequency, int n, int J, int padding_type,
	double* result);

/**
    \fn void CleanData(double * data, double n)

    \param data An 1 x n array with the data to be cleaned
    \param n The size of the data array.

Takes a 1 x n array and preforms the Z-Score Calculation
The array data will be rewritten
*/
void CleanData(double * data, double n);

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
double* GenerateScales(const double minimum_frequency, const double maximum_frequency, const double s_0);

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
double* IdentifyFrequencies(double* scales, int count);

void Convolute(double *data, double *conWindow, double * complexWindow, int data_size, int conSize,
	double* realResult, double* complexResult);

int CWT_Convolution(double *data, double * scales, int data_size, int num_of_scales, 
	double* result);

/**
	\fn int CalculatePaddingSize(int array_size, int pad_flag)
	
	\brief Calculates the size that the padded array should be. 

	\param array_size The cardinal or size of the signal sample
	\param pad_flag The type of padding required: see PAD_FLAG

	\return paadded_size The cardinal or size that the padded array should be

	This function computes the size of the padded array depending on the type of padding specified. 
	It takes the size of the data array, and type of pad, and returns how large the padded array should be.
*/
int CalculatePaddingSize(const int array_size, const int pad_flag);

/**
	\fn int Generate_FFTW_Wisdom(int padded_size)
	
	\brief Analyzes the size of the FFTW arrays and generates the optimal plan. 
	
	\param padded_size The size of the FFT arrays. 

	\return 0 If successful
	\return 1 if unsuccessful

	This function can be used to optimize FFTW. This function will try to find the fastest FFT method based on the size of the array, and will store this information as "FFTW_plan.wise". 

	This function does not need to be used, but it can significantly improve performance if it is.
*/
int Generate_FFTW_Wisdom(int padded_size);

/**
	\fn double Max(double* array, int size)

	\brief A function that finds the Maximum of a given array
	\param array The array to be analyzed
	\param size The size of the array
*/
ARRAY_DATA Max(double * array, int size);

/**
	\fn double Min(double* array, int size)

	\brief A function that finds the minimum of a given array
	\param array The array to be analyzed
	\param size The size of the array
*/
double Min(double* array, int size);

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
	const int n, const int J, const int sampling_frequency,
	double* output);

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
	 fftw_complex* filter_convolution);

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
	fftw_complex* output_data);


int Find_Peaks(double* array, double* frequency, int sampling_frequency, int n, int J);

double* ShortTimeFourierTransform(double * raw_data, int n, int window_size);


#define WINDOW_SIZE 500
int WriteSTFTFile(const double *data, const int x, const int y, int sampling_frequency, const char filename[]);

double* FFT(double * raw_data, int n);

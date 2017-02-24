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

//Global Constants

/*! \def QUAD_ROOT_PI 
	\brief the quad root of pi (\f$ \pi^{-0.25} \f$) computed to machine precision
*/
#define QUAD_ROOT_PI 0.7511255444649425

/** 	
	\var C_SIGMA 
 	\brief A constant needed to compute the Morlet Wavelet precalculated to machine precision

 	\f[ C_\sigma = (1 + e^{-\sigma^2} - 2e^{-\frac{3\sigma^2}{4}})^{-\frac{1}{2}} \f]
*/
#define C_SIGMA 1.00000000000187941284040286940
/**
	\var K_SIGMA
	\brief A constant needed to compute the Morlet Wavelet precalculated to machine eps

	\f[ \kappa_\sigma = e^{-\frac{\sigma^2}{2}} \f]
*/
#define K_SIGMA 1.52299797447126284361366292335e-08

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
#define D_J 0.0625

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
#define FS 2048

/**
	\var DT 
	\brief \f$ \delta t = \frac{1}{f_s} \f$
*/
#define DT 1.0/FS

/**
	\var S0 
	\brief the lowest scale that can be used to compute the CWT \f$ s_0 = 2 \delta t \f$
*/
#define S0 2.0 * DT

#define FREQ 16.0
#define DATA_SIZE 6144 

//Plotting Constants
/**
	\var MAX_FREQUENCY
	\brief The maximum frequency that will be analyzed
*/
#define MAX_FREQUENCY 60.0

/**
	\var MIN_FREQUENCY
	\brief The minimum frequency that will be analyzed
*/
#define MIN_FREQUENCY 2.0
	
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
<<<<<<< HEAD
=======
// #define MAGNITUDE(x,y) (x * x) + (y * y)
>>>>>>> ERSP
#define MAGNITUDE(x,y) sqrt( (x * x) + (y * y) )

void FillData(double * data);
void TestCases(double *data, const int flag);

int  ReadFile(double data[], char filename[]);

int  WriteFile(const double *data, const double *frequency, const int x, const int y, 
	const char filename[]);

int WriteGnuplotScript(const char graph_title[], const char filename[]);

int Plot_PNG(double * data, double * periods, int num_x, int num_y, char graph_title[], 
	const char filename[]);

int Plot(double * data, double * frequency, int num_x, int num_y, int plot_type,
	char graph_title[],
	char filename[]);

int WriteDebug(const double *data, const int length, const int sampling_frequency,
	const char* filename);

int ERSP (double * raw_data, double* scales, const int sampling_frequency, const int n, 
	const int J, int const trials, const int padding_type, 
	double * output);

double CompleteFourierMorlet(double w, const double scale);
double CompleteRealMorlet (double x, double scale);
double CompleteComplexMorlet(double x, double scale);


double CWT_Cosine_Real(double time, double scale);
double CWT_Cosine_Complex(double time, double scale);
double CWT_Dirac_Real(double time, double scale);
double CWT_Dirac_Complex(double time, double scale);

int Wavelet(double* raw_data, double* scales, 
	double sampling_frequency, int n, int J,
	double* result);

void CleanData(double * data, double n);

double* GenerateScales(const double minimum_frequency, const double maximum_frequency, const double s_0);
double* IdentifyFrequencies(double* scales, int count);

void Convolute(double *data, double *conWindow, double * complexWindow, int data_size, int conSize,
	double* realResult, double* complexResult);

int CWT_Convolution(double *data, double * scales, int data_size, int num_of_scales, 
	double* result);

int CalculatePaddingSize(const int array_size, const int pad_flag);
int Generate_FFTW_Wisdom(int padded_size);

double Max(double * array, int size);
double Min(double* array, int size);


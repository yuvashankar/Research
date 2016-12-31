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

/*! \var QUAD_ROOT_PI 
	\brief the quad root of pi (pi^-0.25) computed to machine precision
*/
#define QUAD_ROOT_PI 0.7511255444649425

/** 	\var C_SIGMA 
 	\brief A constant needed to compute the Morlet Wavelet precalculated to machine precision
*/
#define C_SIGMA 1.0000000000018794
/**
	\var K_SIGMA
	\brief A constant needed to compute the Morlet Wavelet precalculated to machine eps
*/
#define K_SIGMA 1.5229979744712628e-08

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

#define D_J 0.125

//Signal Constants
#define FS 2048
#define DT 1.0/FS
#define S0 2.0 * DT
#define FREQ 16.0
#define DATA_SIZE 6144

//Plotting Constants
#define MAX_FREQUENCY 512.0
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
#define MAGNITUDE(x,y) (x * x) + (y * y)

void FillData(double * data);
void TestCases(double *data, int flag);

int ReadFile(double data[], char* filename);
int WriteFile(double *data, double *frequency, int x, int y, const char* filename);
int WriteDebug(double *data, int length, const char* filename);

int ERSP (double * data, double* scales, int sampling_frequency, int n, int J, int trials, 
	double * output);

void Plot(double * data, double * periods, int num_x, int num_y);

double CompleteFourierMorlet(double w, double scale);

int Wavelet(double* raw_data, double* scales, 
	double sampling_frequency, int n, int J,
	double* result);

void CleanData(double * data, double n);

double* GenerateScales(double minimum_frequency, double maximum_frequency);
double* IdentifyFrequencies(double* scales, int count);

void Convolute(double *data, double *conWindow, double * complexWindow, double conSize,
	double* result, double* complexResult);

int CalculatePaddingSize(int array_size, int FLAG);

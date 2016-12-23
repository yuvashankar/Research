
//Header file for Wavelet.c
//By: Vinay Yuvashankar
//Email: yuvashv@mcmaster.ca

//Includes
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#include <omp.h>

//Global Constants
#define QUAD_ROOT_PI 0.7511255444649425 //Precalculated to machine precision Pi^-0.25
#define C_SIGMA 1.0000000000018794 //pow( (1.0 + exp(-W_0_2) - 2.0 * exp(-0.75 * W_0_2)), -0.5 )
#define K_SIGMA 1.5229979744712628e-08
#define FOURIER_WAVELENGTH_FACTOR (4.0 * M_PI)/(W_0 + sqrt(2.0 + W_0_2))

#define W_0 6.0
#define W_0_2 36.0
#define D_J 0.125

//Sample Rate
#define FS 2048.0
#define DT 1.0/FS
#define S0 2.0 * DT

#define MAX_FREQUENCY 1024.0
#define MIN_FREQUENCY 0.2

#define MIN_I FREQ_TO_SCALE(MAX_FREQUENCY)
#define MAX_I FREQ_TO_SCALE(MIN_FREQUENCY)

//Measuring Frequency
#define FREQ 30.0

//Plotting Constants
#define PLOT_OY 100
#define PLOT_OX 100

#define MAX_DATA_SIZE 10000000

#define DATA_SIZE 6144

//Macros
#define FREQ_TO_SCALE(x) floor( ( log2( (W_0) / (S0 * 2 * M_PI * x) ) )/D_J)
#define SCALE_TO_FREQ(x) (W_0)/(x * 2 * M_PI)
#define MAGNITUDE(x,y) (x * x) + (y * y)

void FillData(double * data);
void TestCases(double *data, int flag);

int ReadFile(double data[], char filename[]);
int WriteFile(double *data, double *frequency, int x, int y, char filename[]);

void Plot(double * data,int num_x,int num_y);
double max(double * array, int size);
double min(double * array, int size);

int WriteTestCases(double *data, int length, char filename[]);

double FourierMorlet(double w, double scale, double normal);

double CompleteFourierMorlet(double w, double scale);

int Wavelet(double* raw_data,  double* period, double* scales, 
	double sampling_frequency, int n, int J,
	double* result);

void CleanData(double * data, double n);

int RemoveBaseline(double* data, int num_of_samples, int J, 
	int trials, double sampling_frequency, 
	double* output);

double* GenerateScales(double minimum_frequency, double maximum_frequency);

double CompleteRealMorlet (double x, double scale);
double CompleteComplexMorlet(double x, double scale);

void Convolute(double *data, double *conWindow, double * complexWindow, double conSize,
	double* result, double* complexResult);

int CalculatePaddingSize(int array_size, int FLAG);

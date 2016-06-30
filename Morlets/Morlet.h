//Header file for Morlet.c
//By: Vinay Yuvashankar
//Email: yuvashv@mcmaster.ca

//Includes
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_statistics.h>


#define W_0 5.0
#define W_0_2 25.0

//Global Constants
#define QUAD_ROOT_PI 0.7511255444649425 //Precalculated to machine precision


#define FOURIER_WAVELENGTH_FACTOR 1.0
//(8 * M_PI)/W_0

//Sample Rate
#define FS 1000


#define MAX_FREQUENCY 128.0
#define MIN_FREQUENCY 0.5

//Measuring Frequency
#define FREQ 4.0

#define MAX_DATA_SIZE 10000000   

#define DATA_SIZE 3000

//Morlet Functions
void FillData(double * data);
void TestCases(double *data, int flag);

int ReadFile(double data[], char filename[]);
int WriteFile(double *data, double *frequency, int x, int y, char filename[]);
int WriteTestCases(double *data, int length, char filename[]);

double FourierMorlet(double w, double scale, double k, double cSigma,
	double normal);

int Wavelet(double* raw_data, int sampling_frequency, int n, double dj, double s0, int J, double* result, double* frequency);

double Magnitude (double x, double y);

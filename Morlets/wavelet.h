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
// # define M_PI        3.14159265358979323846    /* pi */
//Global Constants
#define QUAD_ROOT_PI 0.7511255444649425 //Precalculated to machine precision
#define FOURIER_WAVELENGTH_FACTOR (4.0 * M_PI)/(W_0 + sqrt(2.0 + W_0_2))

#define W_0 5.0
#define W_0_2 25.0

//Sample Rate
#define FS 2048.0

#define MAX_FREQUENCY 128.0
#define MIN_FREQUENCY 0.5

//Measuring Frequency
#define FREQ 20.0

#define MAX_DATA_SIZE 10000000 

#define DATA_SIZE 6144

void FillData(double * data);
void TestCases(double *data, int flag);

int AllocateMemory(double *data, double *result, double *frequency, 
	double n, double sampling_frequency, 
	double max_frequency, double dj, double s0);

int ReadFile(double data[], char filename[]);
int WriteFile(double *data, double *frequency, int x, int y, char filename[]);
int WriteTestCases(double *data, int length, char filename[]);

double FourierMorlet(double w, double scale, double normal);

int Wavelet(double* raw_data,  double* frequency, 
	double sampling_frequency, int n, double dj, double s0, int J, double maximum_frequency,
	double* result);

double Magnitude (double x, double y);

void CleanData(double * data, double n);

int RemoveBaseline(double* data, int num_of_samples, int J, 
	int trials, double sampling_frequency);

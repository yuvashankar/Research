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
#include <gsl/gsl_math.h>

//Global Constants
#define quadRootPi 0.7511255444649425 //Precalculated to machine precision

//Sample Rate
#define FS 1000.0

//Measuring Frequency
#define FREQ 19.0

#define DATA_SIZE 504
#define MAX_SCALES 50 * 4

#define MAX_CONV_SIZE 512

//Morlet Functions
void fillData(double * data);
int FillDataComplex(fftw_complex * data);

double Morlet(double x, double w0, double scale);
double FourierMorlet(double w, double w0, double scale);

double ComplexMorlet(double x, double w0, double scale);

double Magnitude (double x, double y);

int createFilter(double* conWindow, double* complexWindow, double frequency);
int CreateComplexFilter(double* conWindow, double frequency);

void convolute(double* data, int conSize, double* conWindow, 
	double* complexWindow, double* result, double* complexResult);
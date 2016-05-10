#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Morlet.h"

void fillData(double * data)
{
	// Fit a FREQ signal at two points
	// double dt = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Priode = %d \n",FS,FREQ,one_peri);
	// double t=0;
	int i;
	for(i=0;i<DATA_SIZE;i++){
		data[i]=0.;
		if((i>200)&(i<400))data[i]=sin( (i-200)*dw+w0);
		//if((i>200)&(i<200+one_peri)) data[i]=sin( (i-200)*dw+w0);
		if((i>1000)&(i<1000+2*one_peri))data[i]=sin( (i-1000)*dw+w0);
		if((i>2000)&(i<2000+3*one_peri))data[i]=sin( (i-2000)*dw+w0);
	}
}

void FillDataComplex(fftw_complex * data)
{
	// Fit a FREQ signal at two points
	// double dt = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Priode = %d \n",FS,FREQ,one_peri);
	// double t=0;
	int i;
	for(i=0;i<DATA_SIZE;i++)
	{
		// data[i][0] = sin(i * dw + w0);
		// data[i][1] = 0.0;
		data[i][0]= 0.; data[i][1] = 0.;
		if((i>200)&(i<400))
		{
			data[i][0]=sin( (i-200)*dw+w0);
			data[i][1] = 0.0;
		}
		if((i>1000)&(i<1000+2*one_peri))
		{
			data[i][0]=sin( (i-200)*dw+w0);
			data[i][1] = 0.0;
		}
		if((i>2500)&(i<2500+3*one_peri))
		{
			data[i][0]=sin( (i-200)*dw+w0);
			data[i][1] = 0.0;
		}
	}
}

double Morlet(double x, double w0, double scale)
{
    const double w02 = w0 * w0;
    // const double sqPi = pow( M_PI, -.25 );
    const double k = exp( -.5 * w02 );

    double normal = 1./sqrt(scale);

    x = x * scale;
    // Now we take the real part of it only !!
    double more =  quadRootPi * exp (-.5 * x*x) * (cos (w0*x) - k) * normal;
    
    return(more);
}

double FourierMorlet(double w, double w0, double scale)
{
	// const double w02 = w0 * w0;
	// const double k = exp(-0.5 * w02);
	// const double cSigma = sqrt((1. + exp(-w02) - 2*exp(-0.75*w02)));

	// double out = exp( -0.5 * (w0 - w)*(w0 - w)) - k * exp(-0.5 * w*w);
	// out = cSigma * out;
	// return(out);

	const double exponent = -0.5 * (scale * w - w0)*(scale * w - w0);
	const double normal = sqrt(scale * w) * quadRootPi * sqrt(DATA_SIZE);
	double out = normal * exp(exponent);
	return(out);
}

double ComplexMorlet(double x, double w0, double scale)
{
	const double w02 = w0 * w0;
    const double sqPi = pow( M_PI, -.25 );
    const double k=exp( -.5 * w02 );

    double normal = 1./sqrt(scale);

    x = x * scale;

    // Now we take the complex part.
    double more =  sqPi * exp(-.5 * x * x) * (sin( w0 * x ) - k ) * normal;
    
    return(more);
}

double Magnitude (double x, double y)
{
	double output = x * x + y * y;
	output = sqrt(output);
	return (output);
}

int createFilter(double* conWindow, double* complexWindow, double frequency)
{
	double signalFrequency = frequency/FS;
	// double dw = 2 * M_PI * signalFrequency;

	int conSize = (int) 1./signalFrequency;
	conSize *=4;
	
	double dt = 1.0/FS;

	double t = 0.0;
	double scale;

	double normal; 
	for (int i = 0; i < MAX_SCALES; ++i)
	{
		t = 0.0; //Gotta reset the time to zero for every scale. 
		scale = pow(2, i);
		// printf("Scale is: %f\n", scale);
		normal = sqrt(2*M_PI*scale/dt);


		for (int j = 0; j < conSize; ++j)
		{
			conWindow[i * MAX_CONV_SIZE + j] = normal * Morlet (t, 5.0 , scale);
			complexWindow[i * MAX_CONV_SIZE + j] = ComplexMorlet (t, 5.0, scale);
			t += dt;
		}
	}

	return(conSize);
}

int CreateComplexFilter(double* conWindow, double frequency)
{
	double scale = 22.0;
	double dt = 1.0/FS;
	// double normal = sqrt (2 * M_PI * scale / dt);
	double df = 1./DATA_SIZE/dt;
	// printf("Dt = %f, Df = %f\n", dt, df);
	printf("Starting df: %f\n", df);

	double value;
	for (int i = 0; i < DATA_SIZE; ++i)
	{
		value = FourierMorlet(i*df, 5.0, scale);
		// value = value * normal;
		
		conWindow[i] = value;
	}
	printf("Ending df: %f\n", DATA_SIZE * df);
	return df;
}

void convolute(double* data, int conSize, double* conWindow, double* complexWindow, double* result, double* complexResult)
{

	for (int i = 0; i < MAX_SCALES; ++i)
	{
		for (int j = 0; j < DATA_SIZE; ++j) //For every element in the data file
		{
			result[i*DATA_SIZE + j] = 0.0;
			complexResult[i*DATA_SIZE + j] = 0.0;

			for (int k = -conSize; k < conSize; ++k) 
			{
				if ( (j - k) > 0)
				{
					if (k >= 0)
					{
						result[i * DATA_SIZE + j] += data[j - k] * conWindow[i * MAX_CONV_SIZE + k];
						complexResult[i * DATA_SIZE + j] += data[j - k] * complexWindow[i * MAX_CONV_SIZE + k];
					}
						

					if (k < 0)
					{
						result[i * DATA_SIZE + j] += data[j - k] * conWindow[i * MAX_CONV_SIZE - k];
						complexResult[i * DATA_SIZE + j] -= data[j - k] * complexWindow[i * MAX_CONV_SIZE - k];
					}
						
				}
			}
		}
	}
}
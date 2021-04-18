 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

//Sample Rate
#define FS 1000.0
//Measuring Frequency
#define FREQ 19.

#define MAX_SCALE 10

#define DATA_SIZE 3000

#define MAX_CONV_SIZE 512

#define MORLET_CENT_FREQ 0.8125 //This was taken from the matlab centfreq command not sure if it still pertains to our situation
#define LOWER_FREQ_BOUND 10.0
#define UPPER_FREQ_BOUND 30.0

double data[DATA_SIZE];

double result[DATA_SIZE * MAX_SCALE];
double complexResult[DATA_SIZE * MAX_SCALE];

double oldResult[DATA_SIZE];
double oldComplex[DATA_SIZE];

int conSize;
double conWindow[MAX_CONV_SIZE];
double complexWindow[MAX_CONV_SIZE];


void fillData(void)
{
	// // double signal[SIGNAL_SIZE];

	// FILE* signalFile = fopen("bic.txt", "r");
	// assert(signalFile != NULL);

	// // obtain file size:
	// fseek (signalFile , 0 , SEEK_END);
	// long lSize = ftell (signalFile);
	// rewind (signalFile);

	// char * buffer = (char*) malloc(sizeof(char)*lSize);
	// assert(buffer);

	// size_t result = fread (buffer, 1, lSize, signalFile);
	// assert(result == lSize);
	// // puts(buffer);


	// char * token = strtok(buffer, "\n");
	
 //    //Get input from txt.
	// int counterVariable = 0;
	// while (token !=NULL)
 //    {
 //    	data[counterVariable] = atof(token);
 //    	double difference = data[counterVariable] - atof(token);
 //    	// printf("signal: %f, difference: %f\n", signal[counterVariable], difference);
 //    	counterVariable++;
 //        token = strtok (NULL, "\n");
 //    }
 //    fclose(signalFile);

	// Fit a FREQ signal at two points
	double dt = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Priode = %d, dw: %f \n",FS,FREQ,one_peri, dw);
	double t=0;
	int i;
	for(i=0;i<DATA_SIZE;i++)
	{
		data[i]=0.;
		if((i>200)&(i<200+one_peri)) 
			data[i]=cos( (i-200)*dw+w0);
		if((i>1000)&(i<1000+2*one_peri))
			data[i]=cos( (i-1000)*0.75 * dw+w0);
		if((i>2000)&(i<2000+3*one_peri))
			data[i]=cos( (i-2000) * 1.5 * dw + w0);
	}
}

double Morlet(double x, double w0, double scale)
{
    const double w02=w0*w0;
    const double sqPi = pow(M_PI,-.25);
    const double k=exp(-.5*w02);

    double normal = 1./sqrt(scale);

    x=x*scale;
    // Now we take the real part of it only !!
    double more =  sqPi*exp(-.5*x*x) * (cos(  w0 * x) -k)*normal;
    
    return(more);
}

double ComplexMorlet(double x, double w0, double scale)
{
	const double w02 = w0 * w0;
    const double sqPi = pow( M_PI, -.25 );
    const double k=exp( -.5 * w02 );

    double normal = 1./sqrt(scale);

    x = x * scale;

    // Now we take the real part of it only !!
    double more =  sqPi * exp(-.5 * x * x) * (sin( w0 * x ) -k ) * normal;
    
    return(more);
}

void createFilter(double frequency, double scale)
{
	double signalFrequency = frequency/FS;
	double dw = 2 * M_PI * signalFrequency;

	conSize = (int) 1./signalFrequency;
	conSize *=4;

	// double dt = 4.0/conSize;
	double dt = 1.0/FS;


	double t = 0;
	for (int i = 0; i < conSize; ++i)
	{
		// double scale = 22.0;
		// //Use Scale of 0.5
		conWindow[i] = Morlet (t, 5.0 , scale);
		complexWindow[i] = ComplexMorlet (t, 5.0, scale);
		// conWindow[i] = value;
		t += dt;
	}
	// printf("time after loop: %f\n", t);
}

double magnitude (double x, double y)
{
	double value = x * x +
		y * y;
	value = sqrt (value);

	return (value);
}

void convolute(void)
{
	double scales[MAX_SCALE] = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560};

	double lowerScale = MORLET_CENT_FREQ / (LOWER_FREQ_BOUND * FS);
	double upperScale = MORLET_CENT_FREQ / (UPPER_FREQ_BOUND * FS);
	//Divide the scale amounts into 50 steps. 
	double step = (upperScale - lowerScale) / MAX_SCALE;
	double scale = lowerScale;
	// printf("Lower Scale Bound: %f, Upper Scale Bound: %f\n", lowerScale, upperScale);
	for (int i = 0; i < MAX_SCALE; ++i)
	{
		// double scale = pow(2, i);
		createFilter(FREQ, scales[i]);
		scale += step;
		for (int j = 0; j < DATA_SIZE; ++j) //For every element in the data file
		{
			result[j] = 0.0;
			for (int k = -conSize; k < conSize; ++k) 
			{
				if ( (j - k) > 0)
				{
					if (k >= 0)
					{
						result[i*j + j] += data[j - k] * conWindow[k];
						complexResult[i*j + j] += data[j - k] * conWindow[k];
					}
						
					if (j < 0)
					{
						result[i*j + j] += data[j - k] * conWindow[-k];
						complexResult[i*j + j] -= data[j - k] * conWindow[-k];
					}
						
				}
			}
		}
	}
}

void oldConvolute(void)
{
	for (int i = 0; i < DATA_SIZE; ++i) //For every element in the data file
	{
		result[i] = 0.0;
		for (int j = -conSize; j < conSize; ++j) 
		{
			if ( (i - j) > 0)
			{
				if (j >= 0)
				{
					oldResult[i] += data[i - j] * conWindow[j];
					oldComplex[i] += data[i - j] * complexWindow[j];
				}
					

				if (j < 0)
				{
					oldResult[i] += data[i - j] * conWindow[-j];
					oldComplex[i] -= data[i - j] * complexWindow[i];
				}
					
			}
		}
	}
}


int main(void)
{
	fillData();
	// createFilter(FREQ);
	convolute();

	// PLOT ALL INTO ONE FILE, EVEN THE FILTER WHICH IS SHORTER
	// FOR THE VALUE IS JUST USE THE ABS
	FILE* out_file=fopen("DATA.log","w");
	
	double value;


	for (int i = 0; i < MAX_SCALE; ++i)
	{
		for (int j = 0; j < DATA_SIZE; ++j)
		{
			value = magnitude(result[i*j + j], result[i*j + j]);

			result[i*j + j] = value;
		}
	}

	// for (int i = 0; i < MAX_SCALE; ++i)
 //    {
 //    	for (int j = 0; j < DATA_SIZE; ++j)
 //    	{
 //    		// value = result[i*j + j] * result[i*j + j] + 
 //    		// 	complexResult[i*j + j] * complexResult[i*j + j];
 //    		// value = sqrt(value);
 //    		value = result[i*j + j];

 //    		fprintf(out_file, "%f\t", value);
 //    	}
    	
 //    	fprintf(out_file, "\n");
 //    }

	for (int i = 0; i < DATA_SIZE; ++i)
	{
		fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i, data[i],
			result[0*i + i] + 0., result[1*i + i] + 5., result[2*i + i] + 10, result[3*i + i] + 15,
			result[4*i + i] + 20, result[5*i + i] + 25, result[6*i + i] + 30, result[7*i + i] + 35,
			result[8*i + i] + 40, result[9*i + i] + 45);
	}

    fclose(out_file);

    createFilter(FREQ, 22.);
    oldConvolute();

    FILE* data_output = fopen("Output.log", "w");
    for(int i=0;i<DATA_SIZE;i++)
    {

    	value = magnitude (oldResult[i], oldComplex[i]);
		if(i<conSize)
			fprintf(data_output,"%d\t%f\t%f\t%f\t%f\t%f\n",i,data[i], fabs(oldResult[i]),conWindow[i], complexWindow[i], value);
		else
			fprintf(data_output,"%d\t%f\t%f\t%f\t%f\t%f\n",i,data[i], fabs(oldResult[i]),0. , 0., value);
	}

	return 0;
}
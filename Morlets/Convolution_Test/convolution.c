#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Sample Rate
#define FS 1000.0

//Measuring Frequency
#define FREQ 19.

#define DATA_SIZE 3000

#define MAX_CONV_SIZE 512

double data[DATA_SIZE];
double result[DATA_SIZE];

int conSize; // THE SIZE OF THE FILTER THAT WE SLIDE
double conWindow[MAX_CONV_SIZE];

void fillData(void)
{

	// Fit a FREQ signal at two points
	double dt = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Priode = %d \n",FS,FREQ,one_peri);
	double t=0;

	for (int i = 0; i < DATA_SIZE; i++)
	{
		data[i]=0.;
		if((i>200)&(i<200+one_peri)) 
			data[i]=cos( (i-200)*dw+w0);

		if((i>1000)&(i<1000+2*one_peri))
			data[i]=cos( (i-1000)*dw+w0);

		if((i>2000)&(i<2000+3*one_peri))
			data[i]=cos( (i-2000)*dw+w0);
	}
}

double Morlet(double x, double w0, double scale)
{
    const double w02 = w0 * w0;
    const double sqPi = pow( M_PI, -.25 );
    const double k=exp( -.5 * w02 );

    double normal = 1./sqrt(scale);

    x = x * scale;
    // Now we take the real part of it only !!
    double more =  sqPi * exp(-.5 * x * x) * (cos( w0 * x ) -k ) * normal;
    
    return(more);
}

void createFilter(double frequency)
{
	double signalFrequency = frequency/FS;
	double dw = 2 * M_PI * signalFrequency;

	conSize = (int) 1./signalFrequency;

	// double start = 0;

	double dt = 4.0/conSize;
	// double t = -start;
	// double value;
	// conSize = (int) 512;
	// printf("dt = %f\n", dt);
	double t = 0;

	for (int i = 0; i < conSize; ++i)
	{
		// //Use Scale of 0.5
		conWindow[i] = Morlet( t, 5 , 0.5);
		// conWindow[i] = value;
		t += dt;

	}
	// printf("time after loop: %f\n", t);
}

void convolute(void)
{
	for (int i = 0; i < DATA_SIZE; ++i) //For every element in the data file
	{
		result[i] = 0.0;
		for (int j = -conSize; j < conSize; ++j) 
		{
			if ( (i - j) > 0)
			{
				if (j >= 0)
					result[i] += data[i - j] * conWindow[j];
				if (j < 0)
					result[i] += data[i - j] * conWindow[-j];
			}
		}
	}
}

int main(void)
{
	fillData();
	createFilter(FREQ);
	convolute();

	// PLOT ALL INTO ONE FILE, EVEN THE FILTER WHICH IS SHORTER
	// FOR THE VALUE IS JUST USE THE ABS
	FILE* out_file=fopen("DATA","w");
	
	
    for (int i = 0;i < DATA_SIZE;i++)
    {
		if(i<conSize)
			fprintf(out_file,"%d\t%f\t%f\t%f\n",i,data[i], fabs(result[i]),conWindow[i]);
		else
			fprintf(out_file,"%d\t%f\t%f\t%f\n",i,data[i], fabs(result[i]),0.);
	}
	fclose(out_file);

	return 0;
}
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Sample Rate
#define FS 1000.0

//Measuring Frequency
#define FREQ 19.0

#define DATA_SIZE 3000
#define MAX_SCALES 10

#define MAX_CONV_SIZE 512

double data[DATA_SIZE];

double result[DATA_SIZE * MAX_SCALES];
double complexResult[DATA_SIZE * MAX_SCALES];



int conSize; // THE SIZE OF THE FILTER THAT WE SLIDE

double conWindow[MAX_CONV_SIZE * MAX_SCALES];
double complexWindow[MAX_CONV_SIZE * MAX_SCALES];

double test[MAX_CONV_SIZE];

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
	int i;
	for(i=0;i<DATA_SIZE;i++){
		data[i]=0.;
		if((i>200)&(i<400))data[i]=sin( (i-200)*dw+w0);
		//if((i>200)&(i<200+one_peri)) data[i]=sin( (i-200)*dw+w0);
		if((i>1000)&(i<1000+2*one_peri))data[i]=sin( (i-1000)*dw+w0);
		if((i>2000)&(i<2000+3*one_peri))data[i]=sin( (i-2000)*dw+w0);
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

double ComplexMorlet(double x, double w0, double scale)
{
	const double w02 = w0 * w0;
    const double sqPi = pow( M_PI, -.25 );
    const double k=exp( -.5 * w02 );

    double normal = 1./sqrt(scale);

    x = x * scale;

    // Now we take the complex part.
    double more =  sqPi * exp(-.5 * x * x) * (sin( w0 * x ) -k ) * normal;
    
    return(more);
}

double Magnitude (double x, double y)
{
	double output = x * x + y * y;
	output = sqrt(output);
	return (output);
}

void createFilter(double frequency)
{
	double signalFrequency = frequency/FS;
	double dw = 2 * M_PI * signalFrequency;

	conSize = (int) 1./signalFrequency;
	conSize *=4;
	
	double dt = 1.0/FS;

	// double scales[MAX_SCALES] = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560};
	double scales[MAX_SCALES];

	double t = 0.0;
	double temp;
	for (int i = 0; i < MAX_SCALES; ++i) //Run ten times.
	{
		t = 0; //Gotta reset the time to zero for every scale. 

		for (int j = 0; j < conSize; ++j)
		{
			temp = Morlet (t, 5.0 , 22.0);
			conWindow[i*j + j] = temp;
			test[i] = temp;

			temp = ComplexMorlet (t, 5.0, 22.0);
			complexWindow[i*j + j] = temp;
			// conWindow[i*j + j] = Morlet (t, 5.0 , 22.0);
			// complexWindow[i*j + j] = ComplexMorlet (t, 5.0, 22.0);

			t += dt;
		}
	}
}

void convolute(void)
{

	for (int i = 0; i < MAX_SCALES; ++i)
	{
		for (int j = 0; j < DATA_SIZE; ++j) //For every element in the data file
		{
			result[i*j + j] = 0.0;
			complexResult[i*j + j] = 0.0;

			for (int k = -conSize; k < conSize; ++k) 
			{
				if ( (j - k) > 0)
				{
					if (k >= 0)
					{
						result[i*j + j] += data[j - k] * conWindow[i*k + k];
						complexResult[i*j + j] += data[j - k] * complexWindow[i*k + k];
					}
						

					if (j < 0)
					{
						result[i*j + j] += data[j - k] * conWindow[-(i*k + k)];
						complexResult[i*j + j] -= data[j - k] * complexWindow[-(i*k + k)];
					}
						
				}
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
	FILE* out_file=fopen("DATA.log","w");
	
	double value;
 //    for (int i = 0;i < DATA_SIZE;i++)
 //    {
 //    	// printf("Data[%d]: %f\n", i, data[i]);
 //    	//a^2 + b^2 = c^2... you should have tried this hours ago.
 //    	const double value = Magnitude(result[i], complexResult[i]);

 //    	// sqrt(result[i] * result[i] + complexResult[i] * complexResult[i]);

	// 	if(i<conSize)
	// 		fprintf (out_file, "%d\t%f\t%f\t%f\t%f\t%f\n", i, data[i], result[i],
	// 			complexResult[i], value, conWindow[i]);
	// 	else
	// 		fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\n", i, data[i], result[i],
	// 			complexResult[i], value, 0.0);
	// }
	//Find the magnitude. 
	for (int i = 0; i < MAX_SCALES; ++i)
	{
		for (int j = 0; j < DATA_SIZE; ++j)
		{
			value = Magnitude(result[i*j + j], complexResult[i*j + j]);
			result[i*j + j] = value;
		}
	}
	//Print into a file. 
	for (int i = 0; i < MAX_CONV_SIZE; ++i)
	{
		value = (conWindow[i] - test[i]);
		fprintf(out_file, "%d\t%f\t%f\n", i, conWindow[i],
			value);
		// fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i, data[i],
		// 	result[0*i + i] + 0., result[1*i + i] + 5., result[2*i + i] + 10, result[3*i + i] + 15,
		// 	result[4*i + i] + 20, result[5*i + i] + 25, result[6*i + i] + 30, result[7*i + i] + 35,
		// 	result[8*i + i] + 40, result[9*i + i] + 45);
	}

	fclose(out_file);

	return 0;
}
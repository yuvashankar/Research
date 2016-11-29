#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//Sample Rate
#define FS 2048.0

//Measuring Frequency
#define FREQ 19.0

#define DATA_SIZE 3000
#define MAX_SCALES 11

#define MORLET_CENT_FREQ 0.8125 //This was taken from the matlab centfreq command not sure if it still pertains to our situation
#define LOWER_FREQ_BOUND 5.0
#define UPPER_FREQ_BOUND 25.0

#define MAX_CONV_SIZE 512

void fillData(double * data)
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
		if(i == 1000) data[i] = 1;
		// if((i>200)&(i<400))data[i]=sin( (i-200)*dw+w0);
		// //if((i>200)&(i<200+one_peri)) data[i]=sin( (i-200)*dw+w0);
		// if((i>1000)&(i<1000+2*one_peri))data[i]=sin( (i-1000)*dw+w0);
		// if((i>2000)&(i<2000+3*one_peri))data[i]=sin( (i-2000)*dw+w0);
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
	double dw = 2 * M_PI * signalFrequency;

	int conSize = (int) 1./signalFrequency;
	conSize *=4;
	
	double dt = 1.0/FS;

	// double lowerScale = (5.0 * M_PI) / LOWER_FREQ_BOUND;
	// double upperScale = (5.0 * M_PI) / UPPER_FREQ_BOUND;
	
	double lowerScale = (MORLET_CENT_FREQ * FS) / LOWER_FREQ_BOUND;
	double upperScale = (MORLET_CENT_FREQ * FS) / UPPER_FREQ_BOUND;
	//Divide the scale amounts into 10 steps. 
	// double step = (upperScale - lowerScale) / MAX_SCALES;
	// double scale = lowerScale;

	// double step = -10.0;
	// double scale = 100.0;

	double t = 0.0;
	double scale = 2 * dt; 

	for (int i = 0; i < MAX_SCALES; ++i) //Run ten times.
	{
		t = 0.0; //Gotta reset the time to zero for every scale. 

		for (int j = 0; j < conSize; ++j)
		{
			conWindow[i * MAX_CONV_SIZE + j] = Morlet (t, 5.0 , scale);
			complexWindow[i * MAX_CONV_SIZE + j] = ComplexMorlet (t, 5.0, scale);
			t += dt;

		}
		scale = (2 * dt) * pow(2, i + 1);
		// scale += step;
		
		double w0 = (MORLET_CENT_FREQ * FS) / scale;
		printf("Scale: %f, w0 at scale: %f\n", scale, w0);
	}
	// printf("Lower Scale: %f, Upper Scale: %f, Scale Step: %f\n", lowerScale, upperScale, step);
	return(conSize);
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

int main(void)
{
	//Allocate Memory for the necessary arrays.
	double *data = malloc(DATA_SIZE * sizeof(double));
	assert(data != NULL);
	
	double *result = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
	double *complexResult = malloc (DATA_SIZE*MAX_SCALES * sizeof(double));
	assert(result != NULL);
	assert(complexResult != NULL);

	double *conWindow = malloc(MAX_CONV_SIZE * MAX_SCALES * sizeof(double));
	double *complexWindow = malloc(MAX_CONV_SIZE * MAX_SCALES * sizeof(double));
	assert(conWindow != NULL);
	assert(complexWindow!= NULL);

	double * test = malloc(MAX_CONV_SIZE* sizeof(double));
	assert (test != NULL);

	int conSize;

	fillData(data);
	
	conSize = createFilter(conWindow, complexWindow, FREQ);

	convolute(data, conSize, conWindow, complexWindow, result, complexResult);

	// PLOT ALL INTO ONE FILE, EVEN THE FILTER WHICH IS SHORTER
	FILE* out_file=fopen("DATA.log","w");
	
	double value;

	for (int i = 0; i < MAX_SCALES; ++i)
	{
		for (int j = 0; j < DATA_SIZE; ++j)
		{
			value = Magnitude(result[i * DATA_SIZE + j], complexResult[i * DATA_SIZE + j]);
			result[i * DATA_SIZE + j] = value;

			// fprintf(out_file, "%f\t", result[i*DATA_SIZE + j]);
		}
		// fprintf(out_file,"\n");
	}

	// Print into a file. 
	for (int i = 0; i < DATA_SIZE; ++i)
	{
		fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i, data[i],
			result[0*DATA_SIZE + i] + 0., result[1*DATA_SIZE + i] + 5., result[2*DATA_SIZE + i] + 10, result[3*DATA_SIZE + i] + 15,
			result[4*DATA_SIZE + i] + 20, result[5*DATA_SIZE + i] + 25, result[6*DATA_SIZE + i] + 30, result[7*DATA_SIZE + i] + 35,
			result[8*DATA_SIZE + i] + 40, result[9*DATA_SIZE + i] + 45);

		// fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i, data[i],
		// 	conWindow[0*MAX_CONV_SIZE + i], conWindow[1*MAX_CONV_SIZE + i], conWindow[2*MAX_CONV_SIZE + i], conWindow[3*MAX_CONV_SIZE + i],
		// 	conWindow[4*MAX_CONV_SIZE + i], conWindow[5*MAX_CONV_SIZE + i], conWindow[6*MAX_CONV_SIZE + i], conWindow[7*MAX_CONV_SIZE + i],
		// 	conWindow[8*MAX_CONV_SIZE + i], conWindow[9*MAX_CONV_SIZE + i]);
	}
	fclose(out_file);

	//Sanitation Engineering
	free(data);
	free(result);
	free(complexResult);
	free(conWindow);
	free(complexWindow);

	return 0;
}
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



#include "Morlet.h"


void fillData(double * data)
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

int FillDataComplex(fftw_complex * data)
{


	// FILE* signalFile = fopen("sst_nino3.dat", "r");
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
	// double sum = 0.0;
	// while (token !=NULL)
 //    {
 //    	data[counterVariable][0] = atof(token);
 //    	data[counterVariable][1] = 0.0;

 //    	counterVariable++;
 //        token = strtok (NULL, "\n");
 //    }
 //    fclose(signalFile);

 //    return (counterVariable);

    // printf("counterVariable: %d\n", counterVariable);

    // Sample Sine Wave.
	// Fit a FREQ signal at two points
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
		if((i>1500)&(i<1500+3*one_peri))
		{
			data[i][0]=sin( (i-200)*dw+w0);
			data[i][1] = 0.0;
		}
	}
	return(DATA_SIZE);
}

int ReadFile(double data[], char filename[])
{
	FILE* signalFile = fopen(filename, "r");
	assert(signalFile != NULL);
	// obtain file size:
	fseek (signalFile , 0 , SEEK_END);
	long lSize = ftell (signalFile);
	rewind (signalFile);

	char * buffer = (char*) malloc(sizeof(char)*lSize);
	assert(buffer != NULL);

	size_t result = fread (buffer, 1, lSize, signalFile);
	assert(result == lSize);
	// puts(buffer);


	char * token = strtok(buffer, "\n");
	
    //Get input from text.
	int counterVariable = 0;
	double sum = 0.0;
	while (token !=NULL)
    {
    	data[counterVariable] = atof(token);
    	counterVariable++;
        token = strtok (NULL, "\n");

    }
    fclose(signalFile);

    return (counterVariable);
}

int WriteFile(double *data, int x, int y, char filename[])
{

    FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    // assert(out_file != NULL);

	for (int i = 0; i < x; ++i)
    {
        for (int j = 0; j < y; ++j)
        {
            // value = Magnitude(result[i*n + j], result[i*n + j]);
            fprintf(out_file, "%f\t", data[i*y + j]);
        }

        fprintf(out_file, "\n");
    }
    fclose(out_file);
    return(0);
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

// double FourierMorlet(double w, double w0, double scale)
// {
// 	// const double w02 = w0 * w0;
// 	// const double k = exp(-0.5 * w02);
// 	// const double cSigma = sqrt((1. + exp(-w02) - 2*exp(-0.75*w02)));

// 	// double out = exp( -0.5 * (w0 - w)*(w0 - w)) - k * exp(-0.5 * w*w);
// 	// out = cSigma * out;
// 	// return(out);
// 	double dw = 2 * M_PI / (DATA_SIZE * 1);
// 	int kplus = 0;
// 	if (w > 0) kplus = 1;

// 	const double exponent = -0.5 * (scale * w - w0)*(scale * w - w0) * kplus;
// 	const double normal = sqrt(scale * dw) * quadRootPi * sqrt(DATA_SIZE);
// 	double out = normal * exp(exponent); 
// 	out = out * kplus; //Heaviside Step Function.

// 	return(out);
// }

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

int CreateComplexFilter(double* conWindow)
{
	// double scale;
	// double value;

	// double dt = 1.0/FS;
	// double df = 1./DATA_SIZE/dt;
	// double s0 = 2 * dt;
	// double dj = 0.25;
	
	// int J = (int) floor((log(DATA_SIZE * dt / s0) / log(2)) / dj);
	
	// double k[DATA_SIZE];
 //    for (int i = 0; i < DATA_SIZE/2 + 1; ++i)
 //    {
 //        k[i] = (i * 2 * M_PI / (DATA_SIZE * 1));
 //        // printf("%f\n", k[i]);
 //    }

 //    int counterVariable = DATA_SIZE/2 - 1;
 //    for (int i = DATA_SIZE/2 + 1; i < DATA_SIZE; ++i)
 //    {
 //        k[i] = -k[counterVariable];
 //        counterVariable -- ;
 //        // printf("%f\n", k[i]);
 //    }

	// for (int i = 0; i < J; ++i)
	// {
	// 	scale = s0 * pow(2, i * dj);

	// 	for (int j = 0; j < DATA_SIZE; ++j)
	// 	{
	// 		value = FourierMorlet(k[j], 5.0, scale);
	// 		conWindow[i * DATA_SIZE + j] = value;
	// 	}
	// }
	
	
	// printf("Ending df: %f\n", DATA_SIZE * df);
	// return (J);
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
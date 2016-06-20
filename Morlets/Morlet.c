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

	for (int i = 0; i < DATA_SIZE; ++i)
	{
		data[i] = 0.0;
	}

	// //Impulse Sample
	// data[2000] = 1.0;

	///Sine Wave Sample
	int i;
	// double t=0;
	for(i=0;i<DATA_SIZE;i++){
		data[i] = sin(i*dw) + sin(i*dw*4);
		// data[i]=0.;
		// if((i>200)&(i<400))data[i]=sin( (i-200)*dw+w0);
		// //if((i>200)&(i<200+one_peri)) data[i]=sin( (i-200)*dw+w0);
		// if((i>1000)&(i<1000+2*one_peri))data[i]=sin( (i-1000)*dw+w0);
		// if((i>2000)&(i<2000+3*one_peri))data[i]=sin( (i-2000)*dw+w0);
	}
}

void TestCases(double *data, int flag)
{

	// Fit a FREQ signal at two points
	// double dt = 1./FS;
	double fsig = FREQ/FS;
	double dw = 2*M_PI*fsig;
	double w0 =  0.01; // A SMALL PHASE SHIFT SO ITS NOT ALL INTERGER ALIGNED
	int one_peri = (int)1./fsig;
	printf("FS  %.2f   Pitch %.f   Discrete Priode = %d \n",FS,FREQ,one_peri);
	
	switch(flag)
	{
		//Impulse
		case 1:
			data[1500] = 1.0;
			break;
		
		//Multiple Sines
		case 2:
			for (int i = 1500; i < 1500 + 2*one_peri; ++i)
			{
				data[i] = sin((i - 1500)* dw + w0) + sin((i - 1500)* 2* dw + w0);
			}
			break;
	}

}

int FillDataComplex(fftw_complex * data)
{
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
	while (token !=NULL)
    {
    	data[counterVariable] = atof(token);
    	counterVariable++;
        token = strtok (NULL, "\n");

    }
    fclose(signalFile);

    return (counterVariable);
}

int WriteFile(double *data, double *frequency, int x, int y, char filename[])
{

    FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

    //Xticks
    fprintf(out_file, "%d\t", x);
    for (int i = 0; i < y; ++i)
    {
    	fprintf(out_file, "%d\t", i);
    }
    fprintf(out_file, "\n");


	for (int i = 0; i < x; ++i)
    {
    	fprintf(out_file, "%f\t", frequency[i]);
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

int WriteTestCases(double *data, int length, char filename[])
{
	FILE* out_file=fopen(filename,"w");
    if (out_file == NULL) return -1;

	for (int i = 0; i < length; ++i)
    {
    	fprintf(out_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", i,
    		data[length*0 + i] + 0., data[length*1 + i] + 5., data[length*2 + i] + 10, 
    		data[length*3 + i] + 15, data[length*4 + i] + 20, data[length*5 + i] + 25, 
    		data[length*6 + i] + 30, data[length*7 + i] + 35, data[length*8 + i] + 40, 
    		data[length*9 + i] + 45, data[length*10+ i] + 50);
    }
    
    fclose(out_file);
    return 0;

}

double Morlet(double x, double w0, double scale)
{
    const double w02 = w0 * w0;
    // const double sqPi = pow( M_PI, -.25 );
    const double k = exp( -.5 * w02 );

    double normal = 1./sqrt(scale);

    x = x * scale;
    // Now we take the real part of it only !!
    double more =  QUAD_ROOT_PI * exp (-.5 * x*x) * (cos (w0*x) - k) * normal;
    
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
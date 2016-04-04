#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define quadRootPi 0.7511255444649425 //Precalculated to machine precision
#define samplingFrequency 500

double Morlet(double t, double w0, double a);

int main(int argc, char const *argv[])
{
	FILE * signalFile;

	long lSize;
	char * buffer;
	char *token;
	size_t result;

	double normal; 

	double signal[500];

	signalFile = fopen("signal.txt", "r");
	assert(signalFile != NULL);

	// obtain file size:
	fseek (signalFile , 0 , SEEK_END);
	lSize = ftell (signalFile);
	rewind (signalFile);

	buffer = (char*) malloc(sizeof(char)*lSize);
	assert(buffer);

	result = fread (buffer, 1, lSize, signalFile);
	assert(result == lSize);
	// puts(buffer);


	token = strtok(buffer, "\n");
	

	int counterVariable = 0;
	while (token !=NULL)
    {
    	signal[counterVariable] = atof(token);
    	double difference = signal[counterVariable] - atof(token);
    	// printf("signal: %f, difference: %f\n", signal[counterVariable], difference);
    	counterVariable++;
        token = strtok (NULL, "\n");
    }

    // for (int tau = 0; i < 500; ++tau) //from the beginning of the signa till the end
    // {
    // 	for (int s = 0; i < count; ++s)
    // 	{
    // 		normal = 1.0/sqrt(s);
    // 		// double value = signal[tau] * Morlet()
    // 	}
    // }

	fclose(signalFile);
	return 0;
}

double Morlet(double t, double w0, double a)
{
	double innerBracket = cos(w0 * t) - exp(-0.5 * w0 * w0);
	double outerExp = exp(-0.5 * t * t);
	double powerNormal = 1.0/pow(a, -0.5);

	double morletWave = powerNormal * quadRootPi * outerExp * innerBracket;
	return(morletWave);
}
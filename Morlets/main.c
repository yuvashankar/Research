#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define quadRootPi 0.7511255444649425 //Precalculated to machine precision
#define samplingFrequency 1000
#define SIGNAL_SIZE 4001
#define SCALE_SIZE 50

double Morlet(double t, double w0, double a);

int main(int argc, char const *argv[])
{
	FILE * signalFile;
	FILE * outputFile;

	long lSize;
	char * buffer;
	char *token;
	size_t result;

	double * output;

	double normal; 

	double signal[SIGNAL_SIZE];

	signalFile = fopen("bic.txt", "r");
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
    
    double dt = 1.0/SIGNAL_SIZE;
    double w0 = samplingFrequency;
    double t = 0.0;

    output = malloc( SIGNAL_SIZE * SCALE_SIZE * sizeof(double));
    assert (output != NULL);

    for (int tau = 0; tau < SIGNAL_SIZE; ++tau) //from the beginning of the signa till the end
    {
    	for (int j = 0; j < SCALE_SIZE; ++j)
    	{
    		//Using a Dydatic Scale for the scales. 
    		double s = pow(2, -j);
    		// double multiplier = pow(2, -j/2);
    		// double mor = Morlet(s*t, w0, )	

    		// double psi = 
    		double value = signal[tau] * Morlet(t, w0, s);
    		output[tau*j + j] = value;
    	}
    	//Increment t
    	t += dt;
    }

    outputFile = fopen("output.txt", "w");
    assert (outputFile != NULL);

    for (int i = 0; i < SIGNAL_SIZE; ++i)
    {
    	for (int j = 0; j < SCALE_SIZE; ++j)
    	{
    		fprintf(outputFile, "%f\t", output[i*j + j]);
    	}
    	fprintf(outputFile, "\n");
    }

    free(buffer);
    free(output);
	fclose(signalFile);
	fclose(outputFile);
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
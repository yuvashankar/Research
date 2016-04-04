#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define quadRootPi 0.7511255444649425 //Precalculated to machine precision

int main(int argc, char const *argv[])
{
	FILE * signalFile;

	long lSize;
	char * buffer;
	char *token;
	size_t result;

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

	fclose(signalFile);
	return 0;
}

double Morlet(double t, double w0, double a)
{

}
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

int main(int argc, char const *argv[])
{
	FILE * signalFile;
	char buffer[10000];
	char *token;
	double signal[500];

	signalFile = fopen("signal.csv", "r");
	assert(signalFile != NULL);

	fgets(buffer, 10000, signalFile);
	assert(buffer);
	// puts(buffer);

	token = strtok(buffer, ",");
	

	int counterVariable = 0;
	while (token !=NULL)
    {
    	signal[counterVariable] = atof(token);
    	counterVariable++;
    	
        token = strtok (NULL, ",");
    }

	fclose(signalFile);
	return 0;
}
/*
FindTriggers.c
This function should take an array input and return the rising and falling edges of the triggers. 

Inputs: 
	StatusInput: The Status Channel Input from the BDF or EDF flie. use the edfread_digital_samples
	numberOfElements: The size of statusInput
	outputBuffer: a 1 x 2 * MAXIMUM_TRIGGERS long long array that is to be populated. 

Outputs:
	outputBuffer: a 1 x 2 * MAXIMUM_TRIGGERS long long array with the odd entries being the
				rising edge and the even entries being the falling edges. 
*/

#include <stdio.h>

void FindTriggers(const int * statusInput, const long long numberOfElements,
						long long * outputBuffer)
{
	long long counterVariable = 0;
	//int i needs to be long long because we are recording it into a long long arrray. 
	for (long long i = 0; i < numberOfElements; ++i)
	{
		if ( (statusInput[i] & 1) != (statusInput[i-1] & 1) )
        {
        	//Record position when you see a change of state. 
        	outputBuffer[counterVariable] = i;
        	counterVariable++;
        }
	}
}
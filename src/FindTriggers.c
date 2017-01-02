/**
    \file "FindTriggers.c"
*/

/**
    \fn long long FindTriggers(const int * statusInput, const long long numberOfRecords,
                        long long * outputBuffer)

    \brief This function should take an array input and return the rising and falling edges of the triggers. 

    \param statusInput: The Status Channel Input from the BDF or EDF flie. use the edfread_digital_samples
	\param numberOfRecords: The size of statusInput
	\param outputBuffer: a 1 x 2 * MAXIMUM_TRIGGERS long long array with the odd entries being the
				rising edge and the even entries being the falling edges. 

    \return counterVariable The number of triggers that were found.                
*/

#include <stdio.h>
#define MAXIMUM_TRIGGERS 1000000

long long FindTriggers(const int * statusInput, const long long numberOfRecords,
						long long * outputBuffer)
{
	long long counterVariable = 0;
	//int i needs to be long long because we are recording it into a long long arrray. 
	int edge = 0;
    for (long long i = 0; i < numberOfRecords; ++i)
    {
    	//Bit and the lower 16 bits to see if any of the triggers have been triggered to on. 
        if ( ((statusInput[i] & 0x0000FFFF) > 0) && (edge == 0) && (statusInput[i-1] != statusInput[i]) ) //Rising Edge Detected.
        {
            outputBuffer[counterVariable] = i;
            counterVariable++;
            edge = 1;
        }

        if (statusInput[i-1] != statusInput[i] && edge == 1) //Falling Edge Detected.
        {
            edge = 0;
        }

        if (counterVariable > MAXIMUM_TRIGGERS)
        	return -1;
    }

    return counterVariable;
}
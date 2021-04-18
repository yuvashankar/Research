 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

/**
    \file "FindTriggers.c"
*/

/**
    \fn int64_t FindTriggers(const int * statusInput, const int64_t numberOfRecords,
                        int64_t * outputBuffer)

    \brief This function should take an array input and return the rising and falling edges of the triggers. 

    \param statusInput: The Status Channel Input from the BDF or EDF flie. use the edfread_digital_samples
	\param numberOfRecords: The size of statusInput
	\param outputBuffer: a 1 x 2 * MAXIMUM_TRIGGERS int64_t array with the odd entries being the
				rising edge and the even entries being the falling edges. 

    \return counterVariable The number of triggers that were found.                
*/

#include <stdio.h>
#include <inttypes.h>
#define MAXIMUM_TRIGGERS 1000000

int FindTriggers(const int * statusInput, const int64_t numberOfRecords,
						int64_t * outputBuffer)
{
	int64_t counterVariable = 0;
    int counterVariable_int = 0;
	//int i needs to be int64_t because we are recording it into a int64_t arrray. 
	int edge = 0;
    for (int64_t i = 0; i < numberOfRecords; ++i)
    {
    	//Bit and the lower 16 bits to see if any of the triggers have been triggered to on. 
        if ( ((statusInput[i] & 0x0000FFFF) > 0) && (edge == 0) && (statusInput[i-1] != statusInput[i]) ) //Rising Edge Detected.
        {
            outputBuffer[counterVariable] = i;
            
            counterVariable++;
            counterVariable_int++;
            edge = 1;
        }

        if (statusInput[i-1] != statusInput[i] && edge == 1) //Falling Edge Detected.
        {
            edge = 0;
        }

        if (counterVariable > MAXIMUM_TRIGGERS)
        	return -1;
    }

    return counterVariable_int;
}
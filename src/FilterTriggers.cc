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
    \file "FilterTriggers.c"
*/
#include <inttypes.h>
#include "processEEG.h"
#include <wavelet.h>


/**
    \fn int FilterTriggers(const int code, const int button, const int numberOfRecords, 
            const int64_t * triggerList,
            const int * readBuffer, 
            int * outputBuffer)

    \brief Filteres the triggers coming in, and finds the specified events

    \param code The code of the trigger list that is needed 
                Possible Inputs: 1, or 2
    \param button The code of the button that is needed
                Possible Inputs 1, or 2
    \param triggerList The list of all of the possible triggers
    \param readBuffer The Status Channel input from the file
    
    \param outputBuffer The buffer that FilterTriggers will populate with the location of the location of the triggers that we're looking for

    \return counterVariable The number of triggers found.
*/
int FilterTriggers(const int code, 
    const int button, 
    const int numberOfRecords, 
    const int64_t * triggerList,
    const int * readBuffer, 
    int * outputBuffer)
{
    int readCode;
    int buttonCode;
    int counterVariable = 0;
	for (int i = 0; i < numberOfRecords; ++i)
    {
        readCode = readBuffer[i] & 255;
        buttonCode = (readBuffer[i] >> 8) & 3;

        if ( (readCode == code) && (buttonCode == button) )
        {
            outputBuffer[counterVariable] = triggerList[i];
            counterVariable++;
        }
    }
    return counterVariable;
}
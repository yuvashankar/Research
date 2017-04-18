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
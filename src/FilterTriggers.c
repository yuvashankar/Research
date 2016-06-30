/*
Filter Triggers: 
Filteres the triggers coming in, and finds the specified events

Inputs: 
        Code: The code of the trigger list that is needed 
                Possible Inputs: 1, or 2
        button: The code of the button that is needed
                Possible Inputs 1, or 2
        triggerList: The list of all of the possible triggers
        readBuffer: The Status Channel Input from the file

Outputs: 
        outputBuffer: The buffer that FilterTriggers will populate with the
            location of the location of the triggers that we're looking for

        counterVariable: Output of the number of triggers found.

*/

int FilterTriggers(const int code, const int button, const int numberOfRecords, 
    const long long * triggerList,
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
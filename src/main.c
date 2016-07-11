/*
 This program should find al the events in a BDF file and create individual
 arrays of the EEG data.

 The first argument must be the file that you are analyzing. 
 */

#include "edflib.h"
#include "processEEG.h"
#include <assert.h>


int main(int argc, char const *argv[])
{
    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag,
        channel,
        filteredTriggerNumber,
        readFlag;

    int32_t* rawStatus,
        *buffer,
        *filteredBuffer;

    double sampleFrequency,
        samplesToRead,
        offset2sec;

    long long numberOfRecords,
        numberOfTriggers,
        offset;

    double * data,
        *tempBuffer;

    long long * triggerList;
    
    struct edf_hdr_struct edfHeader;

    //Opening and reading the file into the edfHeader input the first argument is the input file.
    openFlag = OpenFile(argv[1], &edfHeader);
    assert( openFlag == 0 );

    //Values that I'll need for getting the EDF information.
    handle = edfHeader.handle;
    sampleFrequency = ((double)edfHeader.signalparam[1].smp_in_datarecord /
                       (double)edfHeader.datarecord_duration) * EDFLIB_TIME_DIMENSION;
    numberOfChannels = edfHeader.edfsignals;
    numberOfRecords = edfHeader.signalparam[numberOfChannels - 1].smp_in_file;
    channel = numberOfChannels - 1;
    samplesToRead = (PRE_EVENT_TIME + POST_EVENT_TIME) * sampleFrequency;

    //Allocate Necessary Memory
    rawStatus =         (int*) malloc( numberOfRecords*sizeof(int) );
    triggerList = (long long*) malloc( MAXIMUM_TRIGGERS * sizeof(long long) );
    tempBuffer =     (double*) malloc( samplesToRead * sizeof(double) );
    assert(rawStatus != NULL); assert(triggerList!= NULL); assert(tempBuffer!= NULL);

    //Read the status Signal
    readFlag = edfread_digital_samples(handle, channel, numberOfRecords, rawStatus);
    assert (readFlag != 1);

    //Find Parse the file and find the triggers.
    numberOfTriggers = FindTriggers(rawStatus, numberOfRecords, triggerList);
    assert (numberOfTriggers != -1);

    //Allocating Read Buffer This has to be done after I find all of the triggers. 
    buffer =         (int *) malloc( numberOfTriggers * sizeof(int) );
    filteredBuffer = (int *) malloc( numberOfTriggers * sizeof(int) );
    assert (buffer!= NULL); assert(filteredBuffer != NULL);

    //Read the status channel and put it in the buffer
    for (int i = 0; i < numberOfRecords; ++i)
    {
        edfseek(handle, channel, triggerList[i], EDFSEEK_SET);
        edfread_digital_samples(handle, channel, 1, &buffer[i]);
    }

    //Filter the Triggers to what you want.
    filteredTriggerNumber = FilterTriggers(1, 2, numberOfTriggers, triggerList,
        buffer, filteredBuffer);
    printf("Filtered Trigger Number: %d\n", filteredTriggerNumber);

    //Allocate the necessary memory to copy all of the data into a contigious directory. 
    printf("Will allocate %f Mbs of Memory\n", filteredTriggerNumber * samplesToRead * numberOfChannels * sizeof(double)/1048576);
    data = (double*) malloc( filteredTriggerNumber * samplesToRead * numberOfChannels * sizeof(double) );
    assert(data!= NULL);

    //Load Data from Buffer onto the Data File.
    int dataOffset = 0;
    for (long long i = 0; i < filteredTriggerNumber; ++i)
    {
        for (int j = 0; j <= (numberOfChannels-1) ; j++)
        {
            edfseek(handle, j, filteredBuffer[i], EDFSEEK_SET);
            readFlag = edfread_physical_samples(handle, j, samplesToRead, tempBuffer);
            assert(readFlag != -1);

            dataOffset = (int) ((numberOfChannels-1)*i + (i+j) ) * samplesToRead;
            memcpy(&data[dataOffset], tempBuffer, samplesToRead * sizeof(double));
        }
    }
    printf("File Read into Memory\n");

    // FILE *debug_file = fopen("debug.log", "w");
    // assert(debug_file != NULL);
    // for (int i = 0; i < samplesToRead; ++i)
    // {
    //     fprintf(debug_file, "%d\t%f\n", i, data[i]);
    // }
    // fclose(debug_file);

    //clean up and close up
    edfclose_file(handle);

    free(buffer);
    free(filteredBuffer);
    free(rawStatus);
    free(triggerList);
    free(data);
    free(tempBuffer);

    return 0;
}
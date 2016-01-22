/*
 This program should find al the events in a BDF file and create individual
 arrays of the EEG data.
 */

#include "edflib.h"
#include <assert.h>
#include <string.h>
#include <math.h>

#define PRE_EVENT_TIME 1.0
#define POST_EVENT_TIME 2.0

#define MAXIMUM_TRIGGERS 1000000

int OpenFile(const char* fileName, struct edf_hdr_struct *header);
long long FindTriggers(const int * statusInput, const long long numberOfElements, long long * outputBuffer);

int main(int argc, char const *argv[])
{
    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag,
        channel,
        readFlag;

    int32_t* rawStatus;

    double sampleFrequency,
        samplesToRead,
        offset2sec;

    long long numberOfRecords,
        numberofTriggers,
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
    rawStatus = (int*) malloc(numberOfRecords*sizeof(int));
    triggerList = (long long*) malloc(MAXIMUM_TRIGGERS * sizeof(long long));
    tempBuffer = (double*) malloc(samplesToRead * sizeof(double));
    assert(rawStatus);
    assert(triggerList);
    assert(tempBuffer);

    //Read the status Signal
    readFlag = edfread_digital_samples(handle, channel, numberOfRecords, rawStatus);
    assert (readFlag != 1);

    //Find Parse the file and find the triggers.
    numberofTriggers = FindTriggers(rawStatus, numberOfRecords, triggerList);
    assert (numberofTriggers != -1);

    // for (int i = 0; i < numberofTriggers; ++i)
    // {
    //     edfseek(handle, channel, triggerList[i], EDFSEEK_SET);
    //     edfread_digital_samples(handle, channel, 1, )
    // }

    //Allocate the necessary memory to copy all of the data into a contigious directory. 
    printf(" Will allocate %f Mbs of Memory\n", numberofTriggers*samplesToRead*numberOfChannels*sizeof(double)/1048576);
    data = (double*) malloc (numberofTriggers*samplesToRead*numberOfChannels*sizeof(double));
    assert(data);

    // //Load Data from Buffer onto the Data File.
    // int dataOffset = 0;
    // for (long long i = 0; i < numberofTriggers; ++i)
    // {
    //     for (int j = 0; j <= (numberOfChannels-1) ; j++)
    //     {
    //         edfseek(handle, j, triggerList[i], EDFSEEK_SET);
    //         readFlag = edfread_physical_samples(handle, j, samplesToRead, tempBuffer);
    //         assert(readFlag != -1);

    //         dataOffset = (int) ((numberOfChannels-1)*i + (i+j) ) * samplesToRead;
    //         memcpy(&data[dataOffset], tempBuffer, samplesToRead * sizeof(double));
    //     }
    // }


    //clean up and close up
    edfclose_file(handle);

    free(rawStatus);
    free(triggerList);
    free(data);
    free(tempBuffer);

    return 0;
}
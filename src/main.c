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

#define MAXIMUM_TRIGGERS 1000

int OpenFile(const char* fileName, struct edf_hdr_struct *header);

int main(int argc, char const *argv[])
{
    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag,
        channel,
        edge,
        readFlag;

    int32_t* rawStatus;

    double sampleFrequency;

    long long numberOfRecords,
        status_sample_duration,
        offset;
    
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

    //Allocate Necessary Memory
    rawStatus = (int*) malloc(numberOfRecords*sizeof(int));
    assert(rawStatus);

    //Read the status Signal
    readFlag = edfread_digital_samples(handle, channel, numberOfRecords, rawStatus);
    assert (readFlag != 1);

    edge = 0; 
    for (int i = 0; i < numberOfRecords; ++i)
    {
        offset2sec = ( ((double)i * 
            (double)edfHeader.file_duration) /
            ((double)EDFLIB_TIME_DIMENSION * 
            (double)edfHeader.signalparam[channel].smp_in_file));

        if ( ((rawStatus[i] & 0x0000FFFF) > 0) && (edge == 0) && (rawStatus[i-1] != rawStatus[i]) ) //Rising Edge Detected.
        {
            printf("%f,\t %x, \t %x\n", offset2sec, rawStatus[i], (rawStatus[i] & 0x0000FFFF));
            edge = 1;
        }

        if (rawStatus[i-1] != rawStatus[i] && edge == 1) //Falling Edge Detected.
        {
            edge = 0;
        }
    }

    //clean up and close up
    edfclose_file(handle);

    free(rawStatus);

    return 0;
}
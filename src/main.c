/*
 This program should find al the events in a BDF file and create individual
 arrays of the EEG data.
 */

#include "edflib.h"
#include <assert.h>
#include <string.h>
#include <math.h>

#include <limits.h>

#define PRE_EVENT_TIME 1.0
#define POST_EVENT_TIME 2.0

int OpenFile(const char* fileName, struct edf_hdr_struct *header);

unsigned int_to_int(int k) {
    if (k == 0) return 0;
    if (k == 1) return 1;                       /* optional */
    return (k % 2) + 10 * int_to_int(k / 2);
}

int main(int argc, char const *argv[])
{
    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag,
        readFlag;

    int* rawStatus;

    double sampleFrequency, 
        triggerTime;

    long long numberOfRecords;
    
    struct edf_hdr_struct edfHeader;

    //Opening and reading the file into the hdr input the first argument is the input file.
    openFlag = OpenFile(argv[1], &edfHeader);
    assert( openFlag == 0 );

    //Values that I'll need for getting the EDF information.
    handle = edfHeader.handle;
    sampleFrequency = ((double)edfHeader.signalparam[1].smp_in_datarecord /
                       (double)edfHeader.datarecord_duration) * EDFLIB_TIME_DIMENSION;
    numberOfChannels = edfHeader.edfsignals;
    numberOfRecords = edfHeader.signalparam[numberOfChannels - 2].smp_in_file;

    //Allocate Necessary Memory
    rawStatus = (int*) malloc(numberOfRecords*sizeof(int));
    assert(rawStatus);

    //Read the File
    readFlag = edfread_digital_samples(handle, numberOfChannels - 1, numberOfRecords, rawStatus);
    assert (readFlag != -1);
    for (int i = 0; i < numberOfRecords; ++i)
    {
        printf("%u\n", int_to_int(rawStatus[i]));

        // if ( (rawStatus[i] & 1) != (rawStatus[i-1] & 1) )
        // {
        //    triggerTime = ( ((double)i * 
        //                     (double)edfHeader.file_duration) /
        //                     ((double)EDFLIB_TIME_DIMENSION * 
        //                     (double)edfHeader.signalparam[numberOfChannels - 1].smp_in_file));
        //     printf("Trigger found at %d, at %f, rawStatus[i] = %u \n",i, triggerTime, int_to_int(rawStatus[i]));
        // }
    }

    //clean up and close up
    edfclose_file(handle);
    free(rawStatus);
    return 0;
}
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

int OpenFile(const char* fileName, struct edf_hdr_struct *header);
int AllocateMemory(double* A, int samplesToRead, int numberOfChannels);

int main(int argc, char const *argv[])
{
    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag;

    double sampleFrequency;

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
    numberOfRecords = edfHeader.signalparam[numberOfChannels - 2].smp_in_datarecord;
    printf("%lld\n", numberOfRecords);


    return 0;
}
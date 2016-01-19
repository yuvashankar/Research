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

 #define MAXIMUM_TRIGGERS 1000

int OpenFile(const char* fileName, struct edf_edfHeader_struct *header);

int main(int argc, char const *argv[])
{
    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag,
        bufferSize,
        jump_to_file,
        counterVariable,
        channel,
        readFlag;

    char* rawStatus;

    double* triggerTime;

    double sampleFrequency;

    int edgeDetected = 0;

    long long numberOfRecords,
        status_sample_duration,
        offset;
    
    struct edf_edfHeader_struct edfHeader;

    struct edf_annotation_struct *annotation;

    //Opening and reading the file into the edfHeader input the first argument is the input file.
    openFlag = OpenFile(argv[1], &edfHeader);
    assert( openFlag == 0 );
    edfclose_file(handle);

    //Values that I'll need for getting the EDF information.
    handle = edfHeader.handle;
    sampleFrequency = ((double)edfHeader.signalparam[1].smp_in_datarecord /
                       (double)edfHeader.datarecord_duration) * EDFLIB_TIME_DIMENSION;
    numberOfChannels = edfHeader.edfsignals;
    numberOfRecords = edfHeader.signalparam[numberOfChannels - 1].smp_in_file;
    channel = numberOfChannels - 1;

    //Allocate Necessary Memory
    jump_to_file = (numberOfChannels - 1) * sampleFrequency * 3;
    status_sample_duration = EDFLIB_TIME_DIMENSION / (long long)sampleFrequency;
    rawStatus = (char*) calloc(1, bufferSize);
    triggerTime = (double*) malloc( MAXIMUM_TRIGGERS * sizeof(double) );
    assert(rawStatus);
    assert(triggerTime);
    // printf("bufferSize: %d, numberOfRecords: %lld\n", bufferSize, numberOfRecords);

    //Read the File
    offset = edfHeader->edfHeadersize;
    offset += (edfHeader->edfparam[channel].sample_pntr / edfHeader->edfparam[channel].smp_per_record) * edfHeader->recordsize;
    offset += edfHeader->edfparam[channel].buf_offset;
    offset += ((edfHeader->edfparam[channel].sample_pntr % edfHeader->edfparam[channel].smp_per_record) * bytes_per_smpl);
    smp_per_record = hdr->edfparam[channel].smp_per_record;
    jump = edfHeader->recordsize - (smp_per_record * bytes_per_smpl);

    FILE *file = fopen(argv[1], "r");
    assert(file);

    fseeko(file, offset, SEEK_SET);

    for (long long i = 0; i < numberOfRecords; ++i)
    {
        counterVariable = 0;
        // fseeko(file, (long long)jump_to_file, SEEK_CUR);

        readFlag = fread(rawStatus, bufferSize, 1, file);
        assert(readFlag == 1);

        for (int j = 0; i < bufferSize; i+=3)
        {
            if(rawStatus[j + 2] & 1)
            {
                if (!edgeDetected) //Discovered a new edge
                {
                    annotation = (struct edf_annotation_struct *) calloc (1, sizeof(struct edf_annotation_struct));
                    assert(annotation);
                    triggerTime[counterVariable] = (double)(i * EDFLIB_TIME_DIMENSION) + ((i / 3) * status_sample_duration);
                    counterVariable++;
                    edgeDetected = 1;
                } 
                else
                {
                    edgeDetected = 0;
                }
            }

        }
    }

    //clean up and close up
    // edfclose_file(handle);
    fclose(file);
    free(rawStatus);
    free(triggerTime);
    return 0;
}
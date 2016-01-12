/*
 This program should find al the events in a BDF file and create individual
 arrays of the EEG data.
 */

#include "edflib.h"
#include <assert.h>
#include <string.h>

#define PRE_EVENT_TIME 1.0
#define POST_EVENT_TIME 2.0

int OpenFile(const char* fileName, struct edf_hdr_struct *header);
int AllocateMemory(double* A, int samplesToRead, int numberOfChannels);

int main(int argc, char const *argv[])
{
    
    int handle,
        channel,
        samplesToRead,
        readFlag,
        openFlag,
        numberOfChannels;
    
    double sampleFrequency,
        seekStartTime,
        seekEndTime,
        seekDuration,
        eventTime;
    
    long long offset;
    
    double *buffer,
            *tempBuffer;
    
    struct edf_hdr_struct edfHeader;
    
    
    //Opening and reading the file into the hdr input the first argument is the input file.
    openFlag = OpenFile(argv[1], &edfHeader);
    assert( openFlag == 0 );
    
    handle = edfHeader.handle;
    
    //Get the event time, this is not needed, but I'm implementing it because A10_C1 has no events.
    eventTime = atof(argv[2]);
    assert(eventTime > 0);
    
    //Calculate the sampeling frequency from the first channel. Change if necessary.
    sampleFrequency = ((double)edfHeader.signalparam[1].smp_in_datarecord /
                       (double)edfHeader.datarecord_duration) * EDFLIB_TIME_DIMENSION;
    
    numberOfChannels = edfHeader.edfsignals;
    
    seekStartTime = eventTime - PRE_EVENT_TIME;
    //This is needed to deal with the poossibility of negative time.
    if (seekStartTime < 0)
        seekStartTime = 0;
    
    seekEndTime = eventTime + POST_EVENT_TIME;
    //If the seek Duration is longer than the experiment itself, than it is automatically
    //reset to the length of the experiment.
    seekDuration = seekEndTime - seekStartTime;
    if (seekDuration > edfHeader.file_duration / EDFLIB_TIME_DIMENSION)
        seekDuration = edfHeader.file_duration / EDFLIB_TIME_DIMENSION;
    
    //the number of inputs to read is #of seconds * samples/sec
    samplesToRead = seekDuration * sampleFrequency;
    
    //Allocate the necessary memory.
    buffer = (double*) malloc( samplesToRead * numberOfChannels * sizeof(double));
    tempBuffer = (double*) malloc(samplesToRead * sizeof(double));
    assert(buffer);
    assert (tempBuffer);
    

    
    offset = ((((double)seekStartTime) /
               ((double)edfHeader.file_duration /
                (double)EDFLIB_TIME_DIMENSION)) * ((double)edfHeader.signalparam[channel].smp_in_file));
    

    
    for (int i = 0; i <= (numberOfChannels-1) ; i++)
    {
        edfseek(handle, i, offset, EDFSEEK_SET);
        readFlag = edfread_physical_samples(handle, i, samplesToRead, tempBuffer);
        assert(readFlag != -1);
        
        memcpy(&buffer[samplesToRead*i], tempBuffer, samplesToRead * sizeof(double));
    }
    
    
    printf("\nread %i samples, started at %f seconds from start of file: and ended at %f\n\n",
           samplesToRead, seekStartTime, seekEndTime);
    printf("\n\n");
    
//    //check if the buffer is populated
//    for(int i=0; i<200; i++)
//    {
//        printf("%.0f  \n", buffer[i]);
//    }
    
    //Close file and clean up.
    edfclose_file(handle);
    free(buffer);
    free(tempBuffer);
    
    
    return 0;
}
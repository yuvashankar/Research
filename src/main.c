/*
 This program should find al the events in a BDF file and create individual
 arrays of the EEG data.

 The first argument must be the file that you are analyzing. 
 */

#include "edflib.h"
#include "processEEG.h"
#include "wavelet.h"
#include <assert.h>

int main(int argc, char const *argv[])
{
    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag,
        channel,
        filteredTriggerNumber,
        readFlag,
        waveletFlag,
        writeFlag,
        J; //needed by Wavelet

    int32_t* rawStatus,
        *buffer,
        *filteredBuffer;

    double sampleFrequency,
        samplesToRead,
        dj, dt, s0;

    long long numberOfRecords,
        numberOfTriggers;

    double *tempBuffer,
        *result, //Needed by Wavelet
        *period;

    long long * triggerList;
    
    struct edf_hdr_struct edfHeader;

    //Debug File ... just in case. 
    FILE *debug_file = fopen("Main_debug.log", "w");
    assert(debug_file != NULL);

    //Opening and reading the file into the edfHeader input the first argument is the input file.
    openFlag = OpenFile(argv[1], &edfHeader);
    assert( openFlag == 0 );

    //File Information
    handle = edfHeader.handle;
    sampleFrequency = ((double)edfHeader.signalparam[1].smp_in_datarecord /
                       (double)edfHeader.datarecord_duration) * EDFLIB_TIME_DIMENSION;
    numberOfChannels = edfHeader.edfsignals;
    numberOfRecords = edfHeader.signalparam[numberOfChannels - 1].smp_in_file;
    channel = numberOfChannels - 1;
    samplesToRead = (PRE_EVENT_TIME + POST_EVENT_TIME) * sampleFrequency;

    //Allocate Necessary Memory
    rawStatus =         (int32_t*) malloc( numberOfRecords*sizeof(int32_t) );
    triggerList = (long long*) malloc( MAXIMUM_TRIGGERS * sizeof(long long) );
    
    tempBuffer =     (double*) malloc( samplesToRead * sizeof(double) );
    assert(rawStatus != NULL); assert(triggerList!= NULL); assert(tempBuffer!= NULL);

    //Read the status Signal --> Output to rawStatus
    readFlag = edfread_digital_samples(handle, channel, numberOfRecords, rawStatus);
    assert (readFlag != 1);

    //Find and parse the file and find the triggers -->  Output to TriggerList
    numberOfTriggers = FindTriggers(rawStatus, numberOfRecords, triggerList);
    assert (numberOfTriggers != -1);

    //Allocating Read Buffer
    buffer =         (int32_t *) malloc( numberOfTriggers * sizeof(int32_t) );
    filteredBuffer = (int32_t *) malloc( numberOfTriggers * sizeof(int32_t) );
    assert (buffer!= NULL); assert(filteredBuffer != NULL);

    //Read the status channel when there is a trigger and put it in the buffer
    for (int i = 0; i < numberOfRecords; ++i)
    {
        edfseek(handle, channel, triggerList[i], EDFSEEK_SET);
        readFlag = edfread_digital_samples(handle, channel, 1, &buffer[i]);
        assert(readFlag != -1);
    }

    //Filter the Triggers to what you want.
    filteredTriggerNumber = FilterTriggers(1, 2, numberOfTriggers, triggerList,
        buffer, filteredBuffer);
    printf("Filtered Trigger Number: %d\n", filteredTriggerNumber);

    //Allocate the necessary memory to copy all of the data into a contigious directory. 
    // printf("Will allocate %f Mbs of Memory\n", filteredTriggerNumber * samplesToRead * numberOfChannels * sizeof(double)/1048576);
    // data = (double*) malloc( filteredTriggerNumber * samplesToRead * numberOfChannels * sizeof(double) );
    // assert(data!= NULL);
    
    //Begin Wavelet Analysis
    dj = 0.25;
    dt = 1.0/sampleFrequency;
    s0 = 2 * dt;

    // J = (int) ceil(log2 ( 1.0/(s0 * MIN_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) )/dj);
    J = log2( (samplesToRead * dt)/s0 )/dj;

    //Wavelet Memory Allocations
    result = (double*) malloc(J * samplesToRead * sizeof(double));
    period = (double*) malloc(J * sizeof(double));
    assert(result != NULL); assert(period != NULL);

    // TestCases(tempBuffer, 1);
    edfseek(handle, 0, filteredBuffer[0], EDFSEEK_SET);
    readFlag = edfread_physical_samples(handle, 4, samplesToRead, tempBuffer);
    
    for (int i = 0; i < samplesToRead; ++i)
    {
        fprintf(debug_file, "%d\t%f\n", i, tempBuffer[i]);
    }

    waveletFlag = Wavelet(tempBuffer, period ,
        sampleFrequency, samplesToRead, dj, s0, J, MAX_FREQUENCY,
        result);
    assert(waveletFlag!= -1);

    writeFlag = WriteFile(result, period, J, samplesToRead, "DATA.log");
    assert(writeFlag != -1);

    printf("Wavelet Analysis Done!\n");


    // printf("Beginning Wavelet Analysis\n");
    // //Load Data from Buffer onto the Data File.
    // int dataOffset = 0;
    // // for (long long i = 0; i < filteredTriggerNumber; ++i)
    // for (long long i = 0; i < 1; ++i)
    // {
    //     // for (int j = 0; j <= (numberOfChannels-1) ; j++)
    //     for (int j = 0; j <= 1; ++j)
    //     {
    //         edfseek(handle, j, filteredBuffer[i], EDFSEEK_SET);
    //         readFlag = edfread_physical_samples(handle, j, samplesToRead, tempBuffer);
    //         assert(readFlag != -1);



    //         // dataOffset = (int) ((numberOfChannels-1)*i + (i+j) ) * samplesToRead;
    //         // memcpy(&data[dataOffset], tempBuffer, samplesToRead * sizeof(double));
    //     }
    // }








    //clean up and close up
    edfclose_file(handle);

    free(result);
    free(period);

    free(buffer);
    free(filteredBuffer);
    free(rawStatus);
    free(triggerList);
    // free(data);
    free(tempBuffer);

    fclose(debug_file);

    return 0;
}
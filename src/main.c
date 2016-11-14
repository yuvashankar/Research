/*
 This program should find al the events in a BDF file and create individual
 arrays of the EEG data.

 The first argument must be the file that you are analyzing. 
 */


#include "processEEG.h"
#include <assert.h>
#include <hdf5.h>
#include <float.h>

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
        samplesToRead,
        J; //needed by Wavelet

    int* rawStatus, 
        *buffer,
        *filteredBuffer;

    double sampleFrequency,
        dj, dt, s0;

    long long numberOfRecords,
        numberOfTriggers;

    double *tempBuffer,
        *result, //The Final Result
        *period, //The Corrosponding Frequencies
        *wavelet_result, //What Wavelet spits out
        *baseline_out; //What Baseline spits out

    long long * triggerList;
    
    struct edf_hdr_struct edfHeader;

    //Debug File ... just in case. 
    FILE *debug_file = fopen("Main_debug.log", "w");
    assert(debug_file != NULL);

    //Opening and reading the file into the edfHeader input the first argument is the input file.
    openFlag = OpenFile(argv[1], &edfHeader);
    assert( openFlag == 0 );

    //Get File Information
    handle = edfHeader.handle;
    sampleFrequency = ( ( double )edfHeader.signalparam[1].smp_in_datarecord /
                        ( double )edfHeader.datarecord_duration               ) * EDFLIB_TIME_DIMENSION;
    numberOfChannels = edfHeader.edfsignals;
    numberOfRecords = edfHeader.signalparam[numberOfChannels - 1].smp_in_file;
    channel = numberOfChannels - 1; //The status channel.
    samplesToRead = (PRE_EVENT_TIME + POST_EVENT_TIME) * sampleFrequency;

    //Allocate Necessary Memory...why yes I am OCD.
    rawStatus =         (int*) malloc( numberOfRecords  * sizeof      (int) );
    triggerList = (long long*) malloc( MAXIMUM_TRIGGERS * sizeof(long long) );
    tempBuffer =     (double*) malloc( samplesToRead    * sizeof   (double) );
    assert(rawStatus != NULL); assert(triggerList!= NULL); assert(tempBuffer!= NULL);

    //Read the status Signal --> Output to rawStatus
    readFlag = edfread_digital_samples(handle, channel, numberOfRecords, rawStatus);
    assert (readFlag != 1);

    //Find and parse the file and find the triggers -->  Output to TriggerList
    numberOfTriggers = FindTriggers(rawStatus, numberOfRecords, triggerList);
    assert (numberOfTriggers != -1);

    //Allocating Read Buffer
    buffer =         (int *) malloc( numberOfTriggers * sizeof(int) );
    filteredBuffer = (int *) malloc( numberOfTriggers * sizeof(int) );
    assert (buffer!= NULL); assert(filteredBuffer != NULL);
    printf("Malloc'd buffer and Filtered Buffer\n");



    //Read the status channel when there is a trigger and put it in the buffer
    for (int i = 0; i < numberOfRecords; ++i)
    {
        int value; 
        edfseek(handle, channel, triggerList[i], EDFSEEK_SET);
        
        //Sometimes this function will fail, I need to put a failsafe so that seg faults don't occur. 
        readFlag = edfread_digital_samples(handle, channel, 1, &value);
        assert(readFlag != -1);
        buffer[i] = value;
    }

    //Filter the Triggers to what you want.
    filteredTriggerNumber = FilterTriggers(1, 2, numberOfTriggers, triggerList,
        buffer, filteredBuffer);
    printf("Number of Filtered Triggers Found: %d\n", filteredTriggerNumber);
    
    //Begin Wavelet Analysis
    dj = 0.0625;
    dt = 1.0/sampleFrequency;
    s0 = 2 * dt;

    J = (int) ceil(log2 ( 1.0/(s0 * MIN_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) )/dj);
    int start = (int) floor( log2( 1.0/(s0 * MAX_FREQUENCY * FOURIER_WAVELENGTH_FACTOR) ) /dj);
    // J = log2( (samplesToRead * dt)/s0 )/dj;

    //Wavelet Memory Allocations
    result =         (double*) malloc(J * samplesToRead * sizeof(double));
    wavelet_result = (double*) malloc(J * samplesToRead * sizeof(double));
    baseline_out =   (double*) malloc(J * samplesToRead * sizeof(double));
    period =         (double*) malloc(J *                 sizeof(double));
    assert(result != NULL); assert(period != NULL); assert(wavelet_result != NULL);

    printf("Beginning Wavelet Analysis\n");

    for (int i = 0; i < filteredTriggerNumber; ++i)
    {
        edfseek(handle, 0, filteredBuffer[i], EDFSEEK_SET);
        
        // readFlag = edfread_physical_samples(handle, 4, samplesToRead, tempBuffer);
        TestCases(tempBuffer, 3);
        
        //Preform a Z-Score on the read data. 
        CleanData(tempBuffer, samplesToRead);

        //Preform the Wavelet Analysis
        waveletFlag = Wavelet(tempBuffer, period,
            sampleFrequency, samplesToRead, s0, J, MAX_FREQUENCY,
            wavelet_result);
        assert(waveletFlag!= -1);

        RemoveBaseline(wavelet_result, samplesToRead, J, filteredTriggerNumber, sampleFrequency, baseline_out);
        for (int j = start; j < J; ++j)
        {
            for (int k = 0; k < samplesToRead; ++k)
            {
                result[j * samplesToRead + k] += baseline_out[j * samplesToRead + k];
            }
        }   
    }

    for (int i = start; i < J; ++i)
    {
        for (int j = 0; j < samplesToRead; ++j)
        {
            result[i * samplesToRead + j] = result[i * samplesToRead + j] / filteredTriggerNumber;
        }
    }

    printf("Wavelet Analysis done\n");

    printf("Writing to File\n");

    writeFlag = WriteFile(result, period, J, samplesToRead, "DATA.log");
    assert(writeFlag != -1);

    printf("Freeing Memory and closing files\n");
    //clean up and close up
    edfclose_file(handle);

    free(period);
    free(result);
    free(wavelet_result);
    
    free(filteredBuffer);
    free(buffer);
    
    free(tempBuffer);
    free(triggerList);
    free(rawStatus);
    // free(data);
    fclose(debug_file);
    printf("Memory Cleaned and I'm done\n");

    return 0;
}
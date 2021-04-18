 // MIT License

 // Copyright (c) [2017] [Vinay Yuvashankar]

 // Permission is hereby granted, free of charge, to any person obtaining a copy
 // of this software and associated documentation files (the "Software"), to deal
 // in the Software without restriction, including without limitation the rights
 // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 // copies of the Software, and to permit persons to whom the Software is
 // furnished to do so, subject to the following conditions:

 // The above copyright notice and this permission notice shall be included in all
 // copies or substantial portions of the Software.

 // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 // SOFTWARE.

#include <assert.h>
#include "wavelet.h"
#include <inttypes.h>
#include <gsl/gsl_statistics.h>
#include <omp.h>

int main(int argc, char const *argv[])
{
    if (argc < 2)
    {
        printf("Not enough input arguments please provide a BDF file\n");
        return(-1);
    }

    double t = omp_get_wtime();

    //Variable Declarations
    int handle,
        numberOfChannels,
        openFlag,
        channel,
        filteredTriggerNumber,
        readFlag,
        writeFlag,
        samplesToRead,
        numberOfTriggers,
        read_channel,
        min_i,
        max_i,
        J; //needed by Wavelet

    int *rawStatus, 
        *buffer,
        *value,
        *filteredBuffer;

    double  sampleFrequency,
            dt;

    int64_t numberOfRecords;

    double  *tempBuffer,
            *data,
            *result, //The Final Result
            *scales,
            *frequency; //The Corrosponding Frequencies

    char    graph_title[256],
            file_name[256];

    int64_t * triggerList;
    
    struct edf_hdr_struct edfHeader;

    // Debug File ... just in case. 
    // FILE *debug_file = fopen("Main_debug.log", "w");
    // assert(debug_file != NULL);

    //Opening and reading the file into the edfHeader input the first argument is the input file.
    openFlag = OpenFile(argv[2], &edfHeader);
    read_channel = atoi(argv[1]);
    assert( openFlag == 0 );

    //Get File Information
    handle           = edfHeader.handle;
    sampleFrequency  = (( double ) edfHeader.signalparam[1].smp_in_datarecord /
                        ( double ) edfHeader.datarecord_duration              ) * EDFLIB_TIME_DIMENSION;
    numberOfChannels = edfHeader.edfsignals;
    numberOfRecords  = edfHeader.signalparam[numberOfChannels - 1].smp_in_file;
    channel          = numberOfChannels - 1; //The status channel.
    samplesToRead    = (PRE_EVENT_TIME + POST_EVENT_TIME) * sampleFrequency;

    //Generate Constants for the CWT
    dt               = 1.0/sampleFrequency;
    min_i            = floor( ( log2( (W_0) / (dt * 2 * M_PI * sampleFrequency/2.0) ) )/D_J);
    max_i            = floor( ( log2( (W_0) / (dt * 2 * M_PI * MIN_FREQUENCY) ) )      /D_J);
    assert(min_i > 0); assert(max_i > 0);

    J     = (int) max_i - min_i;

    /*Allocate Necessary Memory*/
    rawStatus   =    (int*)     malloc( numberOfRecords  * sizeof   (int)     );
    triggerList =    (int64_t*) malloc( MAXIMUM_TRIGGERS * sizeof   (int64_t) );
    tempBuffer  =    (double*)  malloc( samplesToRead    * sizeof   (double)  );
    assert(rawStatus != NULL); assert(triggerList!= NULL); assert(tempBuffer!= NULL);

    //Read the status Signal --> Output to rawStatus
    readFlag = edfread_digital_samples(handle, channel, numberOfRecords, rawStatus);
    assert (readFlag != 1);

    //Find and parse the file and find the triggers -->  Output to TriggerList
    numberOfTriggers = FindTriggers(rawStatus, numberOfRecords, triggerList);
    assert (numberOfTriggers != -1);

    //Allocating Read Buffer
    buffer         = (int *) malloc( numberOfTriggers * sizeof(int) );
    filteredBuffer = (int *) malloc( numberOfTriggers * sizeof(int) );
    assert (buffer!= NULL); assert(filteredBuffer != NULL);

    scales         = (double*) malloc ( J * sizeof(double) );
    frequency      = (double*) malloc ( J * sizeof(double) );
    assert(scales != NULL); assert(frequency != NULL);

    //Wavelet Memory Allocations
    data           = (double*) malloc( samplesToRead * filteredTriggerNumber * sizeof(double));
    result         = (double*) malloc( J * samplesToRead * sizeof(double));
    assert(result != NULL); assert(data != NULL);

    /*Main Program Begins*/
    GenerateScalesAndFrequency(min_i, max_i, dt, scales, frequency);
    
    //Read the status channel when there is a trigger and put it in the buffer
    value = (int*) malloc( sizeof(int) );
    for (int i = 0; i < numberOfTriggers; ++i)
    {
        edfseek(handle, channel, triggerList[i], EDFSEEK_SET);

        readFlag = edfread_digital_samples(handle, channel, 1, value);
        assert(readFlag != -1);
        buffer[i] = *value;
    }

    //Filter the Triggers to what you want.
    filteredTriggerNumber = FilterTriggers(1, 2, numberOfTriggers, triggerList,
                                            buffer, filteredBuffer);
    printf("Number of Filtered Triggers Found: %d\n", filteredTriggerNumber);


    printf("Beginning Wavelet Analysis\n");
    for (int i = 0; i < filteredTriggerNumber; ++i)
    {
        edfseek(handle, 0, filteredBuffer[i], EDFSEEK_SET);
        
        readFlag = edfread_physical_samples(handle, read_channel, samplesToRead, tempBuffer);
        // TestCases( tempBuffer, 5, 16.0 , sampleFrequency, samplesToRead);

        // Preform a Z-Score on the read data. 
        CleanData(tempBuffer, samplesToRead);

        for (int j = 0; j < samplesToRead; ++j)
        {
            data[i * samplesToRead + j] = tempBuffer[j];
        }
    }
    
    ERSP (data, scales, sampleFrequency, samplesToRead, 
                J, filteredTriggerNumber, 1, 
                result);

    printf("Writing to File\n");
    sprintf(graph_title, "Channel %d", read_channel);
    sprintf(file_name,   "Channel_%d.log", read_channel);

    writeFlag = WriteFile(result, frequency, J, samplesToRead, sampleFrequency,
                          file_name);
    WriteGnuplotScript(graph_title, file_name);

    assert(writeFlag != -1);

    printf("Freeing Memory and closing files\n");
    edfclose_file(handle);

    free(data);   free(result);
    free(scales); free(frequency);
    free(buffer); free(filteredBuffer);
    free(rawStatus); free(triggerList); free(tempBuffer);
    
    t = omp_get_wtime() - t;

    printf("Execution Time = %f\n", t);

    return 0;
}
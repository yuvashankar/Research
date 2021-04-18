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

/**
	\file "processEEG.h"
	\brief The files needed to open and process BDF Files.
*/

#include "edflib.h"
#include "wavelet.h"

//Constants
/**
	\val PRE_EVENT_TIME
	\brief The amount of seconds before the stimulus
*/
#define PRE_EVENT_TIME 1.0

/**
	\val POST_EVENT_TIME
	\brief The amount of seconds after the stimulus
*/
#define POST_EVENT_TIME 2.0

#define MAXIMUM_TRIGGERS 1000000

//Functions that are going to be used.
int OpenFile(const char* fileName, struct edf_hdr_struct *header);

int FindTriggers(const int * statusInput, const int64_t numberOfElements, 
	int64_t * outputBuffer);

int FilterTriggers(const int code, 
    const int button, 
    const int numberOfRecords, 
    int64_t * triggerList,
    const int * readBuffer, 
    int * outputBuffer);

// int RemoveBaseline(double* data, int num_of_samples, int J, 
// 	int trials, double sampling_frequency, 
// 	double* output);
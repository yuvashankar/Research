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
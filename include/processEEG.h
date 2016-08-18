#include "edflib.h"
#include "wavelet.h"

//Constants
#define PRE_EVENT_TIME 1.0
#define POST_EVENT_TIME 2.0

#define MAXIMUM_TRIGGERS 1000000

//Functions that are going to be used.
int OpenFile(const char* fileName, struct edf_hdr_struct *header);

long long FindTriggers(const int * statusInput, const long long numberOfElements, 
	long long * outputBuffer);

int FilterTriggers(const int code, 
    const int button, 
    const int numberOfRecords, 
    const long long * triggerList,
    const int * readBuffer, 
    int * outputBuffer);

void CleanData(double * data, double n);

void RemoveBaseline(double* data, double num_of_samples, int J, 
	int trials, double sampling_frequency);
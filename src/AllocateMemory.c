/*
Given the samples to read and the number of channels, this function will allocate an nxm matrix
*/

#include <assert.h>

int AllocateMemory(double *A, int samplesToRead, int numberOfChannels)
{
	A = (double*) malloc(samplesToRead * numberOfChannels);
	assert(A);
    return 0;
}
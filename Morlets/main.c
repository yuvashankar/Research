#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "Morlet.h"

int main(void)
{
    //Allocate Memory for the necessary arrays.
    double *data = malloc(DATA_SIZE * sizeof(double));
    assert(data != NULL);
    
    double *result = malloc(DATA_SIZE * MAX_SCALES * sizeof(double));
    double *complexResult = malloc (DATA_SIZE*MAX_SCALES * sizeof(double));
    assert(result != NULL);
    assert(complexResult != NULL);

    double *conWindow = malloc(MAX_CONV_SIZE * MAX_SCALES * sizeof(double));
    double *complexWindow = malloc(MAX_CONV_SIZE * MAX_SCALES * sizeof(double));
    assert(conWindow != NULL);
    assert(complexWindow!= NULL);

    FILE* out_file=fopen("DATA.log","w");

    //This used to be the convolution stuff commented out to test out the Fourier Convolution stuff. 
    // int conSize;
    // fillData(data);
    // conSize = createFilter(conWindow, complexWindow, FREQ);
    // convolute(data, conSize, conWindow, complexWindow, result, complexResult);
    // FILE* out_file=fopen("DATA.log","w");
    
    // double value;
    // for (int i = 0; i < MAX_SCALES; ++i)
    // {
    //     for (int j = 0; j < DATA_SIZE; ++j)
    //     {
    //         value = Magnitude(result[i * DATA_SIZE + j], complexResult[i * DATA_SIZE + j]);
    //         result[i * DATA_SIZE + j] = value;

    //         fprintf(out_file, "%f\t", result[i*DATA_SIZE + j]);
    //     }
    //     fprintf(out_file,"\n");
    // }
    fclose(out_file);

    //Sanitation Engineering
    free(data);
    free(result);
    free(complexResult);
    free(conWindow);
    free(complexWindow);

    return 0;
}
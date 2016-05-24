#include "Morlet.h"

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double * data, *result;

    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    double dj, dt, s0, J;
    dt = 1.0/FS;
    dj = 1;
    s0 = 2 * dt;
    // J = 7/dj;
    J = 10;

    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    assert(data != NULL); assert(result != NULL);

    //populate the data array
    fillData(data);
    

    int out  = Wavelet(data, dt, n, dj, s0, J, result);
    // printf("Wavelet Flag = %d\n", out);

    int writeFlag = WriteFile(result, J, n, "DATA.log");

    free(data); free(result);
    fclose(out_file);

    return 0;
}
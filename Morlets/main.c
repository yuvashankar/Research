#include "Morlet.h"

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double *data, *result, *frequency;

    //Open the Output file
    // FILE* out_file=fopen("DATA.log","w");
    // assert(out_file != NULL);

    double dj, dt, s0, J;
    dt = 1.0/FS;
    dj = 0.25;
    s0 = 2 * dt;
    J = 7/dj; //Where did this 7 come from? dunno m8.
    // printf("dt = %f, dj = %f, s0 = %f, J = %f\n", dt, dj, s0, J);

    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    frequency = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL);

    // populate the data array
    fillData(data);
    // TestCases(data, 1);
    // char *read_file_name = "sst_nino3.dat";
    // ReadFile(data, read_file_name);

    int out  = Wavelet(data, dt, n, dj, s0, J, result, frequency);
    int writeFlag = WriteFile(result, frequency, J, n, "DATA.log");

    free(data); free(result);
    // fclose(out_file);
    return 0;
}
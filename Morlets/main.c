#include "Morlet.h"
#include "math.h"

int main(void)
{
    //Size of Data
    int n = DATA_SIZE;
    double *data, *result, *frequency;

    //Open the Output file
    FILE* out_file=fopen("DATA.log","w");
    assert(out_file != NULL);

    double dj, dt, s0, J, maxScale, minScale;
    dt = 1.0/FS;
    dj = 0.25;
    s0 = 2 * dt;


    J = (log2(n * dt)/s0)/dj;

    maxScale = MAX_FREQUENCY * (W_0 + sqrt(2 + W_0_2)) / 4 * M_PI;
    minScale = MIN_FREQUENCY * (W_0 + sqrt(2 + W_0_2)) / 4 * M_PI;


    printf("dt = %f, dj = %f, s0 = %f, J = %f\n", dt, dj, s0, J);
    printf("maxScale = %f, minScale = %f\n", maxScale, minScale);

    data = malloc(n * sizeof(double));
    result = malloc(J * n * sizeof(double));
    frequency = malloc(J * sizeof(double));
    assert(data != NULL); assert(result != NULL); assert(frequency != NULL);

    // populate the data array
    fillData(data);
    // TestCases(data, 2);

    // int out  = Wavelet(data, dt, n, dj, s0, J, result, frequency);
    // int writeFlag = WriteFile(result, frequency, J, n, "DATA.log");

    free(data); free(result);
    fclose(out_file);
    return 0;
}
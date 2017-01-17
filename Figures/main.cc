#include "wavelet.h"
#include <math.h>
#include <omp.h>

#include <gsl/gsl_statistics.h>

int main(void)
{
    //Start the timer!
    double t = omp_get_wtime();

    //Initialize the necessary arrays.
    double *data;


    FILE * figure_out = fopen("data.log", "w");
    double start = 4.0;
    double time = -start;
    double dt = (double) 1.0/FS;

    int n = (int) (2 * start) * FS;
    printf("n = %d\n", n);
    // int n = DATA_SIZE;
    // const int J = (int) MAX_I - MIN_I;

    //Memory Allocations
    data    =  (double*) malloc(n *     sizeof(double));
    assert(data != NULL);

    for (int i = 0; i < n; ++i)
    {
        data[i] = CompleteRealMorlet(time, 1.0);
        time += dt;
        fprintf(figure_out, "%f\t%f\n", time, data[i]);
    }

    fclose(figure_out);



    // WriteFile(result, frequency, J, n, "DATA.log");
    // WriteDebug(data, n, FS, "data.log");

    //Free up Memory
    free(data);
    
    //Stop and print the timer. 
    t = omp_get_wtime() - t;
    printf("ERSP Execution Time: %f\n", t);
    return 0;
}